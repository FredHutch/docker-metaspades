#!/usr/bin/python
"""Align a set of reads against a reference database with DIAMOND, and save the results."""

import os
import sys
import time
import json
import uuid
import boto3
import shutil
import logging
import argparse
import traceback
import subprocess
from helpers.io_helpers import prokka_parser


def run_cmds(commands, retry=0, catchExcept=False):
    """Run commands and write out the log, combining STDOUT & STDERR."""
    logging.info("Commands:")
    logging.info(' '.join(commands))
    p = subprocess.Popen(commands,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()
    exitcode = p.wait()
    if stdout:
        logging.info("Standard output of subprocess:")
        for line in stdout.split('\n'):
            logging.info(line)
    if stderr:
        logging.info("Standard error of subprocess:")
        for line in stderr.split('\n'):
            logging.info(line)

    # Check the exit code
    if exitcode != 0 and retry > 0:
        msg = "Exit code {}, retrying {} more times".format(exitcode, retry)
        logging.info(msg)
        run_cmds(commands, retry=retry - 1)
    elif exitcode != 0 and catchExcept:
        msg = "Exit code was {}, but we will continue anyway"
        logging.info(msg.format(exitcode))
    else:
        assert exitcode == 0, "Exit code {}".format(exitcode)


def ena_url(accession):
    """See https://www.ebi.ac.uk/ena/browse/read-download for URL format."""
    url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
    folder1 = accession[:6]
    url = "{}/{}".format(url, folder1)
    if len(accession) > 9:
        if len(accession) == 10:
            folder2 = "00" + accession[-1]
        elif len(accession) == 11:
            folder2 = "0" + accession[-2:]
        elif len(accession) == 12:
            folder2 = accession[-3:]
        else:
            logging.info("This accession is too long: " + accession)
            assert len(accession) <= 12
        url = "{}/{}".format(url, folder2)

    # Add the accession to the URL
    url = "{}/{}/{}".format(url, accession, accession)

    return url


def get_sra(accession, temp_folder):
    """Get the FASTQ for an SRR accession via ENA (falling back to SRA)."""

    # Download from ENA via FTP
    url = ena_url(accession)

    logging.info("Base info for downloading from ENA: {}".format(url))
    # There are three possible file endings
    file_endings = ["_1.fastq.gz", "_2.fastq.gz", ".fastq.gz"]
    # Try to download each file
    for end in file_endings:
        run_cmds(["curl",
                  "-o", os.path.join(temp_folder, accession + end),
                  url + end], catchExcept=True)

    # Local paths for each of the three possible file endings
    r0_fp = "{}/{}{}".format(temp_folder, accession, ".fastq.gz")
    r1_fp = "{}/{}{}".format(temp_folder, accession, "_1.fastq.gz")
    r2_fp = "{}/{}{}".format(temp_folder, accession, "_2.fastq.gz")

    # If the forward and reverse reads were downloaded, return that pair
    if os.path.exists(r1_fp) and os.path.exists(r2_fp):
        # Return a tuple of filepaths, and a bool indicating paired-end reads
        logging.info("Both forward and reverse reads were found")
        return (r1_fp, r2_fp), True
    # If the file was downloaded with no _1/_2, return that
    elif os.path.exists(r0_fp):
        logging.info("Only a single set of unpaired reads were found")
        return r0_fp, False
    # Hedging against oddly incomplete data, return either R1 or R2, if alone
    elif os.path.exists(r1_fp):
        logging.info("Only a single set of unpaired reads were found")
        return r1_fp, False
    elif os.path.exists(r2_fp):
        logging.info("Only a single set of unpaired reads were found")
        return r2_fp, False

    # If none of those URLs downloaded, fall back to trying NCBI
    logging.info("No data was found on ENA, falling back to SRA")
    run_cmds([
        "fastq-dump",
        "--outdir",
        "--split-files",
        temp_folder, accession])

    # Local paths for each of the three possible file endings
    r0_fp = "{}/{}{}".format(temp_folder, accession, ".fastq")
    r1_fp = "{}/{}{}".format(temp_folder, accession, "_1.fastq")
    r2_fp = "{}/{}{}".format(temp_folder, accession, "_2.fastq")

    # If the forward and reverse reads were downloaded, return that pair
    if os.path.exists(r1_fp) and os.path.exists(r2_fp):
        # Return a tuple of filepaths, and a bool indicating paired-end reads
        logging.info("Both forward and reverse reads were found")
        return (r1_fp, r2_fp), True
    # If the file was downloaded with no _1/_2, return that
    elif os.path.exists(r0_fp):
        logging.info("Only a single set of unpaired reads were found")
        return r0_fp, False
    # Hedging against oddly incomplete data, return either R1 or R2, if alone
    elif os.path.exists(r1_fp):
        logging.info("Only a single set of unpaired reads were found")
        return r1_fp, False
    elif os.path.exists(r2_fp):
        logging.info("Only a single set of unpaired reads were found")
        return r2_fp, False

    # If no files were downloaded, throw an error
    msg = "File could not be downloaded from SRA: {}".format(accession)
    raise Exception(msg)


def get_reads_from_url(input_str, temp_folder, interleaved=False):
    """Get a set of reads from a URL -- return the downloaded filepath."""
    logging.info("Getting reads from {}".format(input_str))

    filename = input_str.split('/')[-1]
    local_path = os.path.join(temp_folder, filename)

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    if not input_str.startswith(('s3://', 'sra://', 'ftp://')):
        logging.info("Treating as local path")
        msg = "Input file does not exist ({})".format(input_str)
        assert os.path.exists(input_str), msg
        logging.info("Making symbolic link in temporary folder")
        os.symlink(input_str, local_path)
        return local_path, interleaved

    # Get files from AWS S3
    if input_str.startswith('s3://'):
        logging.info("Getting reads from S3")
        run_cmds([
            'aws', 's3', 'cp', '--quiet', '--sse',
            'AES256', input_str, temp_folder
            ])

    # Get files from an FTP server
    elif input_str.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', temp_folder, input_str])

    # Get files from SRA
    elif input_str.startswith('sra://'):
        accession = filename
        logging.info("Getting reads from SRA: " + accession)
        local_path, interleaved = get_sra(accession, temp_folder)

        # If a set of paired end data was found in SRA, the `local_path`
        # will be a tuple of two files, and `interleaved` will be true.
        # In all other circumstances, when `interleaved` is true, `local_path`
        # will represent the location of a single interleaved file

    else:
        raise Exception("Did not recognize prefix to fetch reads: " + input_str)

    return local_path, interleaved


def return_results(out, read_prefix, output_folder, temp_folder):
    """Write out the final results as a JSON object and write them to the output folder."""
    # Make a temporary file
    temp_fp = os.path.join(temp_folder, read_prefix + '.json')
    with open(temp_fp, 'wt') as fo:
        json.dump(out, fo)
    # Compress the output
    run_cmds(['gzip', temp_fp])
    temp_fp = temp_fp + '.gz'

    if output_folder.startswith('s3://'):
        # Copy to S3
        run_cmds([
            'aws', 's3', 'cp',
            '--quiet', '--sse', 'AES256',
            temp_fp, output_folder
            ])
        os.unlink(temp_fp)
    else:
        # Copy to local folder
        run_cmds(['mv', temp_fp, output_folder])


def run_metaspades(input_str,
                   output_folder,
                   threads=16,
                   max_mem=240,
                   temp_folder='/mnt/temp',
                   overwrite=False,
                   interleaved=False):
    """Assemble a set of reads using metaSPAdes."""

    # Record the start time
    start_time = time.time()

    # Use the read prefix to name the output and temporary files
    read_prefix = input_str.split('/')[-1]

    # Check to see if the output already exists, if so, skip this sample
    output_fp = output_folder.rstrip('/') + '/' + read_prefix + '.json.gz'
    if output_fp.startswith('s3://'):
        # Check S3
        logging.info("Ensure that the output path doesn't already exist on S3")
        bucket = output_fp[5:].split('/')[0]
        prefix = '/'.join(output_fp[5:].split('/')[1:])
        client = boto3.client('s3')
        results = client.list_objects(Bucket=bucket, Prefix=prefix)
        if 'Contents' in results:
            if overwrite:
                logging.info("Overwriting output ({})".format(output_fp))
            else:
                logging.info(
                    "Output already exists, skipping ({})".format(output_fp))
                return
    else:
        # Check local filesystem
        if os.path.exists(output_fp):
            if overwrite:
                logging.info("Overwriting output ({})".format(output_fp))
            else:
                logging.info(
                    "Output already exists, skipping ({})".format(output_fp))
                return

    # Get the reads
    read_fp, interleaved = get_reads_from_url(
        input_str,
        temp_folder,
        interleaved=interleaved
    )

    # If a set of paired end data was found in SRA, the `local_path`
    # will be a tuple of two files, and `interleaved` will be true.
    # In all other circumstances, when `interleaved` is true, `local_path`
    # will represent the location of a single interleaved file

    logging.info("\nRunning metaSPAdes\n")
    # Check for paired-end data
    if interleaved:
        logging.info("Paired-end information present")
        if isinstance(read_fp, tuple):
            logging.info("Two separate files for forward and reverse reads")
            # TWO PAIRED END READ FILES
            r1_fp, r2_fp = read_fp
            run_cmds([
                "metaspades.py",
                "-1", r1_fp,         # Forward read
                "-2", r2_fp,         # Reverse read
                "-o", temp_folder,   # Output folder
                "-t", str(threads),  # Threads
                "-m", str(max_mem)   # Maximum memory used
                ])
        else:
            logging.info("Single file, interleaved forward and reverse reads")
            # ONE INTERLEAVED READ FILE
            run_cmds([
                "metaspades.py",
                "--12", read_fp,     # Interleaved forward and reverse reads
                "-o", temp_folder,   # Output folder
                "-t", str(threads),  # Threads
                "-m", str(max_mem)   # Maximum memory used
                ])
    else:
        logging.info("Single file, unpaired reads")
        # ONE SINGLE-END READ FILE
        # NOTE: MetaSPAdes does not currently work with single-end data
        # Therefore falling back to SPAdes
        run_cmds([
            "spades.py",
            "-s", read_fp,       # Single-end reads file
            "-o", temp_folder,   # Output folder
            "-t", str(threads),  # Threads
            "-m", str(max_mem)   # Maximum memory used
            ])
    logging.info("Done with assembly")

    # Make sure the output file exists (assembled scaffolds)
    scaffold_fp = os.path.join(temp_folder, "scaffolds.fasta")
    assert os.path.exists(scaffold_fp)
    logging.info("Scaffolds written to {}".format(scaffold_fp))

    # Run Prokka on the assembled scaffold sequences
    prokka_folder = os.path.join(temp_folder, "prokka")
    logging.info("Running Prokka")
    run_cmds([
        "prokka",
        "--outdir",
        prokka_folder,
        "--prefix",
        read_prefix,
        "--cpus",
        str(threads),
        "--metagenome",
        scaffold_fp
        ])

    # Collect the results
    logging.info("Parsing the output")
    output = prokka_parser(prokka_folder, read_prefix)

    # Read in the logs
    logging.info("Reading in the logs")
    logs = open(log_fp, 'rt').readlines()
    spades_logs = {
        "logs": open(
            os.path.join(temp_folder, "spades.log"), "rt"
            ).readlines(),
        "params": open(
            os.path.join(temp_folder, "params.txt"), "rt"
            ).readlines(),
        "warnings": open(
            os.path.join(temp_folder, "warnings.log"), "rt"
            ).readlines()
        }

    # Make an object with all of the results
    out = {
        "input_path": input_str,
        "input": read_prefix,
        "output_folder": output_folder,
        "logs": logs,
        "spades_logs": spades_logs,
        "results": output,
        "time_elapsed": time.time() - start_time
    }

    # Write out the final results as a JSON object
    # and copy them to the output folder
    return_results(out, read_prefix, output_folder, temp_folder)
    logging.info("Done\n\n\n----------------\n\n\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Assemble a set of reads using metaSPAdes.
    """)

    parser.add_argument("--input",
                        type=str,
                        help="""Location for input file(s). Comma-separated.
                                (Supported: sra://, s3://, or ftp://).""")
    parser.add_argument("--output-folder",
                        type=str,
                        help="""Folder to place results.
                                (Supported: s3://, or local path).""")
    parser.add_argument("--interleaved",
                        action="store_true",
                        help="""Treat input as interleaved by default \
                                (ignored for SRA datasets).""")
    parser.add_argument("--overwrite",
                        action="store_true",
                        help="""Overwrite output files. Off by default.""")
    parser.add_argument("--threads",
                        type=int,
                        default=16,
                        help="Number of threads to use aligning.")
    parser.add_argument("--max-mem",
                        type=int,
                        default=240,
                        help="Maximum memory to use (in Gb).")
    parser.add_argument("--temp-folder",
                        type=str,
                        default='/share',
                        help="Folder used for temporary files.")

    args = parser.parse_args()

    # Check that the temporary folder exists
    assert os.path.exists(args.temp_folder)

    # Set a random string, which will be appended to all temporary files
    random_string = uuid.uuid4()

    # Make a temporary folder within the --temp-folder with the random string
    temp_folder = os.path.join(args.temp_folder, str(random_string))
    # Make sure it doesn't already exist
    msg = "Collision, {} already exists".format(temp_folder)
    assert os.path.exists(temp_folder) is False, msg
    # Make the directory
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = '{}/log.txt'.format(temp_folder)
    logFormatter = logging.Formatter('%(asctime)s %(levelname)-8s [run.py] %(message)s')
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to file
    fileHandler = logging.FileHandler(log_fp)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Also write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # Align each of the inputs and calculate the overall abundance
    for input_str in args.input.split(','):
        logging.info("Processing input argument: " + input_str)
        # Capture in a try statement
        try:
            run_metaspades(input_str,
                           args.output_folder,
                           threads=args.threads,
                           max_mem=args.max_mem,
                           temp_folder=temp_folder,
                           overwrite=args.overwrite,
                           interleaved=args.interleaved)
        except:
            # There was some error
            # Capture the traceback
            logging.info("There was an unexpected failure")
            exc_type, exc_value, exc_traceback = sys.exc_info()
            for line in traceback.format_tb(exc_traceback):
                logging.info(line)

            # Delete any files that were created in this process
            logging.info("Deleting temporary folder: {}".format(temp_folder))
            shutil.rmtree(temp_folder)

            # Exit
            logging.info("Exit type: {}".format(exc_type))
            logging.info("Exit code: {}".format(exc_value))
            sys.exit(exc_value)

    # Delete any files that were created in this process
    logging.info("Deleting temporary folder: {}".format(temp_folder))
    shutil.rmtree(temp_folder)

    # Stop logging
    logging.info("Done")
    logging.shutdown()

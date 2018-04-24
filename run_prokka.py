#!/usr/bin/python
"""Annotate an assembly with Prokka."""

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
from helpers.io_helpers import output_file_parser, truncate_fasta_headers


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


def get_file_from_url(input_str, temp_folder):
    """Get a file from a URL -- return the downloaded filepath."""
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
        return local_path

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

    return local_path


def return_results(out, sample_name, output_folder, temp_folder):
    """Write out the final results as a JSON object and write them to the output folder."""

    # Keep a list of files that need to be copied
    to_copy = []

    # Make a temporary file for the JSON format data
    json_temp_fp = os.path.join(temp_folder, sample_name + '.prokka.json')
    # Remember to copy this file
    to_copy.append(json_temp_fp)
    # Open up the file handle
    with open(json_temp_fp, 'wt') as fo:
        # Write out the JSON data
        json.dump(out, fo)

    # Write out the sequences and annotations as their own files
    for entry_name, file_ending in [
        ("gff", "gff"),
        ("proteins", "fastp"),
        ("genbank", "gbk")
    ]:
        # Make a file path
        temp_fp = os.path.join(temp_folder, sample_name + '.' + file_ending)
        # Remember to copy the file
        to_copy.append(temp_fp)
        # Open up the file path
        with open(temp_fp, "wt") as fo:
            # Write out each item in this list
            for item in out["results"][entry_name]:
                # Write out the contigs and proteins in FASTA format
                if entry_name in ["contigs", "proteins"]:
                    header, seq = item
                    fo.write(">{}\n{}\n".format(header, seq))
                # Everything else is written out straight
                else:
                    fo.write(item)

    # Compress each of the files with GZIP
    for fp in to_copy:
        run_cmds(['gzip', fp])
    to_copy = [s + ".gz" for s in to_copy]

    # Write to S3
    if output_folder.startswith('s3://'):
        for fp in to_copy:
            run_cmds([
                'aws', 's3', 'cp',
                '--quiet', '--sse', 'AES256',
                fp, output_folder
                ])
            os.unlink(fp)
    # Copy to a local folder
    else:
        # Copy to local folder
        for fp in to_copy:
            run_cmds(['mv', fp, output_folder])


def run_prokka(input_str,
               sample_name,
               output_folder,
               threads=16,
               temp_folder='/mnt/temp'):
    """Annotate an assembly with Prokka."""

    # Record the start time
    start_time = time.time()

    # Get the assembly
    scaffold_fp = get_file_from_url(
        input_str,
        temp_folder
    )

    # If the file is gzipped, decompress it
    if scaffold_fp.endswith(".gz"):
        logging.info("Decompressing contig FASTA")
        run_cmds(["gunzip", scaffold_fp])
        scaffold_fp = scaffold_fp[:-3]

    logging.info("")
    logging.info("")
    logging.info("Running Prokka")
    logging.info("")
    logging.info("")

    # Run Prokka on the assembled scaffold sequences
    prokka_folder = os.path.join(temp_folder, "prokka")
    logging.info("Running Prokka")
    run_cmds([
        "prokka",
        "--outdir",
        prokka_folder,
        "--prefix",
        sample_name,
        "--cpus",
        str(threads),
        "--metagenome",
        scaffold_fp
        ])

    # Collect the results
    logging.info("Parsing the output")
    output = output_file_parser(prokka_folder, sample_name)

    # Make an object with all of the results
    out = {
        "input_path": input_str,
        "input": sample_name,
        "output_folder": output_folder,
        "results": output,
        "time_elapsed": time.time() - start_time
    }

    # Write out the final results as a JSON object
    # and copy them to the output folder
    return_results(out, sample_name, output_folder, temp_folder)
    logging.info("Done")
    logging.info("")
    logging.info("")
    logging.info("---------------")
    logging.info("")
    logging.info("")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Annotate an assembly using Prokka.
    """)

    parser.add_argument("--input",
                        type=str,
                        required=True,
                        help="""Location for input file. Scaffolds in FASTA format.
                                (Supported: sra://, s3://, or ftp://).""")
    parser.add_argument("--sample-name",
                        type=str,
                        required=True,
                        help="""Prefix for output files.""")
    parser.add_argument("--output-folder",
                        type=str,
                        required=True,
                        help="""Folder to place results.
                                (Supported: s3://, or local path).""")
    parser.add_argument("--threads",
                        type=int,
                        default=16,
                        help="Number of threads to use.")
    parser.add_argument("--temp-folder",
                        type=str,
                        default='/share',
                        help="Folder used for temporary files.")

    args = parser.parse_args()

    # Check that the temporary folder exists
    assert os.path.exists(args.temp_folder)

    # Set a random string, which will be appended to all temporary files
    random_string = str(uuid.uuid4())[:8]

    # Make a temporary folder within the --temp-folder with the random string
    temp_folder = os.path.join(args.temp_folder, str(random_string))
    # Make sure it doesn't already exist
    msg = "Collision, {} already exists".format(temp_folder)
    assert os.path.exists(temp_folder) is False, msg
    # Make the directory
    os.mkdir(temp_folder)

    # Set up logging
    log_fp = '{}/log.txt'.format(temp_folder)
    logFormatter = logging.Formatter('%(asctime)s %(levelname)-8s [run_prokka.py] %(message)s')
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
        # Make a new temp folder for just this sample
        sample_temp_folder = os.path.join(temp_folder, str(uuid.uuid4())[:8])
        os.mkdir(sample_temp_folder)
        # Capture in a try statement
        try:
            run_prokka(input_str,
                       args.sample_name,
                       args.output_folder,
                       threads=args.threads,
                       temp_folder=sample_temp_folder)
            logging.info("Deleting temp folder {}".format(sample_temp_folder))
            shutil.rmtree(sample_temp_folder)
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

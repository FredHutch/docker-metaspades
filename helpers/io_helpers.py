#!/usr/bin/python

import os
import logging
from Bio.SeqIO.FastaIO import SimpleFastaParser


def truncate_fasta_headers(fasta_fp, max_len):
    """Truncate the headers in a FASTA file, modifying the file in place."""
    assert os.path.exists(fasta_fp)
    # Read in the FASTA information
    with open(fasta_fp, "rt") as f:
        fasta = [record for record in SimpleFastaParser(f)]

    # Write out to the same file location
    with open(fasta_fp, "wt") as fo:
        for record in fasta:
            header, seq = record
            # Check the header length
            if len(header) > max_len:
                # Truncate
                header = header[:max_len]
            # Write out in FASTA format
            fo.write(">{}\n{}\n".format(header, seq))


def read_tsv(fp):
    """Read a TSV (with header) and format it as a list of dicts."""
    with open(fp, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    output = []
    with open(fp, "rt") as f:
        for ix, line in enumerate(f):
            if ix > 0 and len(line) > 1:
                line = line.rstrip("\n").split("\t")
                if len(line) <= len(header):
                    output.append(dict(zip(header, line)))
                else:
                    logging.info("Line too long, skipping")
                    logging.info("\t".join(line))


def prokka_parser(folder, prefix):
    """Collect the set of results from prokka."""
    output = {}

    # Collect the FASTA records for contigs, transcripts, and proteins
    for tag, file_ending in [
        ("contigs", ".fna"),
        ("transcripts", ".ffn"),
        ("proteins", ".faa"),
    ]:
        filepath = os.path.join(folder, prefix + file_ending)
        assert os.path.exists(filepath)
        # Read in the FASTA
        logging.info("Reading in {}".format(filepath))
        records = [r for r in SimpleFastaParser(open(filepath, "rt"))]
        output[tag] = records

    # Record the features from the TSV
    features_fp = os.path.join(folder, prefix + ".tsv")
    assert os.path.exists(features_fp)
    logging.info("Reading in {}".format(features_fp))
    output["features"] = read_tsv(features_fp)

    # Also read in the Genbank file
    genbank_fp = os.path.join(folder, prefix + ".gbk")
    assert os.path.exists(genbank_fp)
    logging.info("Reading in {}".format(genbank_fp))
    with open(genbank_fp, "rt") as f:
        output["genbank"] = f.readlines()

    # Also read in the GFF file
    gff_fp = os.path.join(folder, prefix + ".gff")
    assert os.path.exists(gff_fp)
    logging.info("Reading in {}".format(gff_fp))
    with open(gff_fp, "rt") as f:
        output["gff"] = f.readlines()

    return output

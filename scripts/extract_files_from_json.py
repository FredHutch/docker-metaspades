#!/usr/bin/python
"""Extract the contigs, proteins, and genbank files from an output JSON."""

# The typical output for this workflow is a JSON containing all data
# This script will extract the contig, protein, and annotation information
# each in their own isolated file.

import os
import json
import gzip
import argparse


def extract_files_from_json(json_fp):
    """Extract the contig, protein, and annotation data from this JSON."""
    # Read in the data
    if json_fp.endswith(".json"):
        file_prefix = json_fp.replace(".json", "")
        dat = json.load(open(json_fp, "rt"))
    elif json_fp.endswith(".json.gz"):
        file_prefix = json_fp.replace(".json.gz", "")
        dat = json.load(gzip.open(json_fp, "rt"))

    # Make sure the output file paths do not exist
    for suffix in [".fasta", "fastp", ".gbk"]:
        fp = file_prefix + suffix
        msg = "File already exists, exiting ({})".format(fp)
        assert os.path.exists(fp) is False, msg

    # Write out the contigs
    with open(file_prefix + '.fasta', "wt") as fo:
        for header, seq in dat["results"]["contigs"]:
            fo.write(">{}\n{}\n".format(header, seq))
    print("Wrote out contigs to {}".format(file_prefix + '.fasta'))

    # Write out the proteins
    with open(file_prefix + '.fastp', "wt") as fo:
        for header, seq in dat["results"]["proteins"]:
            fo.write(">{}\n{}\n".format(header, seq))
    print("Wrote out proteins to {}".format(file_prefix + '.fastp'))

    # Write out the genbank annotations
    with open(file_prefix + '.gbk', "wt") as fo:
        for line in dat["results"]["genbank"]:
            fo.write(line)
    print("Wrote out GenBank annotations to {}".format(file_prefix + '.gbk'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Extract the contigs, proteins, and genbank files from an output JSON.
    """)

    parser.add_argument("--json",
                        type=str,
                        help="""JSON file containing assembly output.""")

    args = parser.parse_args()

    assert os.path.exists(args.json)
    msg = "Input file must end with either .json or .json.gz"
    assert args.json.endswith(".json") or args.json.endswith(".json.gz"), msg

    extract_files_from_json(args.json)

# docker-metaspades
Docker image running metaSPAdes

[![Docker Repository on Quay](https://quay.io/repository/fhcrc-microbiome/metaspades/status "Docker Repository on Quay")](https://quay.io/repository/fhcrc-microbiome/metaspades)

The purpose of this repository is to build a Docker image running metaSPAdes.
It will also contain a wrapper script that will be included within the Docker image
that will make it more convenient to run within an HPC (e.g. Slurm) or 'cloud'
(e.g. AWS Batch).


### Wrapper script options

#### --input

Specifies the set of FASTQ reads that will be aligned. Supports files from SRA, S3, or FTP. 
Use the file prefix to specify the source (`s3://`, `sra://`, or `ftp://`). Note that for 
SRA, users can just provide the accession (e.g. `sra://SRR123456`).

#### --output-folder

Folder to place the output in, supporting either `s3://` or a local path. 
Output files will take the form of `<prefix>.json.gz`, where `<prefix>` 
is the SRA accession (if specified), or otherwise the prefix of the input file from S3 or ftp. 

#### --interleaved

Treat the input data as interleaved by default. This is not considered when you provide
an SRA accession, as the paired-end nature of the data is guessed automatically.

#### --threads

Number of threads used by metaSPAdes during assembly, defaults to 16.

#### --max-mem

Maximum amount of memory to use, in Gb

#### --overwrite

Overwrite any output data, if present. Default behavior is to skip if a file is present
in the output filepath.

#### --temp-folder

The path to the folder used for temporary space. This should be somewhere with a bit of space
for all of the temporary files created by metaSPAdes.


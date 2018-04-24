#!/bin/bash

set -e

test_image(){
	
	img_tag=$1

	[[ "${#img_tag}" == "0" ]] && echo "Please specify image" && return

	[[ "$(docker run $img_tag echo True)" != "True" ]] && echo "Tag not found ($img_tag)" && return

	docker run \
		-v $PWD/output:/output \
		-v $PWD/temp:/share \
		-v $HOME:/root \
		--rm \
		$img_tag \
			run.py \
				--input /share/example.fastq.gz \
				--sample-name example \
				--output-folder /output \
				--temp-folder /share \
				--overwrite \
				--max-mem 8 \
				--threads 4 \
				--interleaved && \
	echo "Completed assembly + annotation" && \
	docker run \
		-v $PWD/output:/output \
		-v $PWD/temp:/share \
		-v $HOME:/root \
		-v ~/.aws/credentials:/root/.aws/credentials \
		--rm \
		$img_tag \
			run_metaspades.py \
				--input s3://fh-pi-fredricks-d/lab/Sam_Minot/data/temp/metaspades_tests/example.fastq.gz \
				--sample-name example \
				--output-folder s3://fh-pi-fredricks-d/lab/Sam_Minot/data/temp/metaspades_tests/ \
				--temp-folder /share \
				--overwrite \
				--max-mem 8 \
				--threads 4 \
				--interleaved && \
	echo "Completed assembly only" && \
	docker run \
		-v $PWD/output:/output \
		-v $PWD/temp:/share \
		-v $HOME:/root \
		-v ~/.aws/credentials:/root/.aws/credentials \
		--rm \
		$img_tag \
			run_prokka.py \
				--input s3://fh-pi-fredricks-d/lab/Sam_Minot/data/temp/metaspades_tests/example.fasta.gz \
				--sample-name example \
				--output-folder s3://fh-pi-fredricks-d/lab/Sam_Minot/data/temp/metaspades_tests/ \
				--temp-folder /share \
				--threads 4 && \
	echo "Completed annotation only"
}

test_image $1


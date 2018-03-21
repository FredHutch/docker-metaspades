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
				--input sra://SRR6174207 \
				--output-folder /output \
				--temp-folder /share \
				--overwrite \
				--max-mem 8 \
				--threads 4
}

test_image $1


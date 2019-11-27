#!/bin/bash
# convert a tsv that has one read per line into two files
# See https://gist.github.com/nathanhaigh/3521724 deinterleave.sh
# Usage:
#  cat fastq.tsv | tsv_to_fastq.sh R1.fastq.gz R2.fastq.gz compress
# Set up some defaults
GZIP_OUTPUT=0
PIGZ_COMPRESSION_THREADS=10

# If the third argument is the word "compress" then we'll compress the output using pigz
if [[ $3 == "compress" ]]; then
  GZIP_OUTPUT=1
fi

if [[ ${GZIP_OUTPUT} == 0 ]]; then
   tee >(cut -f 1-4 | tr "\t" "\n" > $1) | cut -f 5-8 | tr "\t" "\n" > $2
else
   tee >(cut -f 1-4 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $1) | cut -f 5-8 | tr "\t" "\n" | pigz --best --processes ${PIGZ_COMPRESSION_THREADS} > $2
fi
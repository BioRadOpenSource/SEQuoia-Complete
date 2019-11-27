#!/bin/bash
# Combine two fastq files so that there is one read per line
# https://gist.github.com/nathanhaigh/3521724 deinterleave_fastq.sh
# Usage:
#  fastq_to_tsv.sh biorad_v2_R1_test_100k.fastq.gz biorad_v2_R2_test_100k.fastq.gz > fastq.tsv
paste <(pigz -d -c $1 | paste - - - -) <(pigz -d -c $2 | paste - - - -)
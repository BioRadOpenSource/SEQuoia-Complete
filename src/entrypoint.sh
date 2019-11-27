#!/bin/bash
set -x

# enable debugging
[[ "$NXF_DEBUG_ENTRY" ]] && set -x

# wrap cli args with single quote to avoid wildcard expansion
cli=''; for x in "$@"; do cli+="'$x' "; done

nextflow -log /work/nextflow.log run /opt/2d-complete-all/main.nf $cli

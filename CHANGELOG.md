# v1.2.0 

## config

- Updated time out of based processing to 12 hours to match rest of pipeline

## DSL2

- Updated to match coding standards of DSL2 and not change main functions of the code outside of making compatible with current nextflow
- Moved the scripts to appropriate modules configuration 

### main

- Changed structure to make DSL2 using modules format and passing values to processes instead of global vars
- Import modules to match DSL2 standards

## conda

- Created a dev and locked conda env to keep constant the packages need to run the software with out depeneding on apt-get or version of ubuntu that might be retired

## Docker 

- Updated Dockerfile to match the current standards of using a locked conda yaml as the basis for the build to improve stability and reproducibility
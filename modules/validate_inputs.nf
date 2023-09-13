/*
Validation of the input reads
    - reads should be fastq file format
    - less that 500 million reads input to the toolkit
*/
process validateInputs {
    tag "Validation on $sample_id"
    publishDir "${params.outDir}/validation", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    
    script:
    """
    python3 /opt/biorad/src/validate.py $reads 2>&1 > validation.log
    """
}
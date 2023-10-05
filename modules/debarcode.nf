/*
Parse the reads and pull out barcodes to check the quality
The quality includes length, N's present 
Then creates a report to inform the user
*/
process debarcode{
    label 'mid_cpu'
    tag "debarcode on $sample_id"
    publishDir "${params.outDir}/debarcode", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*_R1.fastq.gz'), emit: reads 
    path('debarcode_stats.txt'), emit: report
    
    script:
    """
    bash /opt/biorad/src/fastq_to_tsv.sh $reads \
        | parallel --pipe python3 /opt/biorad/src/debarcode.py \
        | tee >(awk '/^MM/{bad=bad+1}/^@/{good=good+1}END{print "Total Reads: " good + bad; print "Good Reads: " good; print "Bad Reads: " bad}' > debarcode_stats.txt) \
        | grep -ve '^MM' \
        | bash /opt/biorad/src/tsv_to_fastq.sh ${sample_id}_debarcoded_R1.fastq.gz ${sample_id}_debarcoded_R2.fastq.gz compress
    """
}
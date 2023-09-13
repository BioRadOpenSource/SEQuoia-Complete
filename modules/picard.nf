/*
Calculate alignment statistics based on BAM files
Count how many reads fall on different features in the genome
*/
process picardAlignSummary {
    label 'low_memory'
    tag "picardAlignSummary on $sample_id"
    publishDir "${params.outDir}/picardAlignSummary", mode: 'copy'

    input:
    tuple val(sample_id), path(bams)
    path(refFlatFile)
    path(ribosomalIntervalFile)

    output:
    path('rna_metrics.txt'), emit: report_picard

    script:
    (bam, bai) = bams
    strand = params.reverseStrand ? "SECOND_READ_TRANSCRIPTION_STRAND" : "FIRST_READ_TRANSCRIPTION_STRAND"
    if(task.memory.toGiga() ==0){
	    picard_mem = task.memory.toBytes()
    }else{
	    picard_mem = task.memory.toGiga() 
    }
    """
    picard CollectRnaSeqMetrics -I $bam \
    -O rna_metrics.txt \
    -REF_FLAT $refFlatFile \
    -STRAND $strand \
    -RIBOSOMAL_INTERVALS $ribosomalIntervalFile \
    -Xmx${picard_mem}g
    """
}
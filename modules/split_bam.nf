/*
Use bedtools along with bed file to specify which regions are miRNA
to extracate the miRNA reads in to thier own BAM file for counting 
*/
process splitBamMi {
    label 'low_memory'
    tag "splitBamMi on $sample_id"
    publishDir "${params.outDir}/splitBamMi", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bams)
    path(miRNABedFile)

    output:
    tuple val(sample_id), path('out.miRNAs.bam'), emit: mirna_bam_ch

    script:
    (bam, bai) = bams
    """
    bedtools intersect -a $bam -b $miRNABedFile -f 1 -s > out.miRNAs.bam
    """
}
/*
Use bedtools along with bed file to specify which regions are not miRNA
to extracate the longer RNA reads in to thier own BAM file for counting 
The difference is the -v which does the opposite of the splitBamMi process
*/
process splitBamLong {
    label 'low_memory'
    tag "splitBamLong on $sample_id"
    publishDir "${params.outDir}/splitBamLong", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bams)
    path(miRNABedFile)

    output:
    tuple val(sample_id), path('out.longRNAs.bam'), emit: longrna_bam_ch

    script:
    (bam, bai) = bams
    """
    bedtools intersect -a $bam -b $miRNABedFile -f 1 -v -s > out.longRNAs.bam
    """
}
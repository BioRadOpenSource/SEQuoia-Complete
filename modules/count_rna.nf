/*
Using feature counts with provided gtf files to count which reads overlap
with known genomic feaures for genes tRNA etc 
*/
process countMicroRNA {
    label 'mid_cpu'
    label 'low_memory'
    tag "countMicroRNA on $sample_id"
    publishDir "${params.outDir}/microRNACounts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path(miRNAgtfFile)
    val(geneId)

    output:
    path("gene_counts_miRNA*"), emit: report_miRNACounts
    tuple val(sample_id), path("gene_counts_miRNA"), emit: counts_ch_mi 

    script:
    strand = params.reverseStrand ? "-s 2" : "-s 1"
    just_bam = bam[0]
    """
    featureCounts -T $task.cpus \
    --primary \
    -O \
    -M \
    -t exon \
    -g product $strand \
    -Q $params.minMapqToCount \
    -a $miRNAgtfFile \
    -o ./gene_counts_miRNA \
    -R BAM $just_bam
    """
}

process countLongRNA {
    label 'mid_cpu'
    label 'low_memory'
    tag "countLongRNA on $sample_id"
    publishDir "${params.outDir}/longRNACounts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam)
    path(longRNAgtfFile)
    val(geneId)
    
    output:
    tuple val(sample_id), path("gene_counts_longRNA"), emit: counts_ch_long
    path("gene_counts_longRNA*"), emit: report_longRNACounts

    script:
    strand = params.reverseStrand ? "-s 2" : "-s 1"
    just_bam = bam[0] 
    """
    featureCounts -T $task.cpus \
    --primary \
    -M \
    -t exon \
    -g $geneId $strand \
    -Q $params.minMapqToCount \
    -a $longRNAgtfFile \
    -o ./gene_counts_longRNA \
    -R BAM $just_bam
    """
}
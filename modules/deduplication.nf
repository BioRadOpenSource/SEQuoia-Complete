/*
Run UMI tools to find reads with the exact same UMI at the same
genomic coordinates in the genome to remove them prior to counting
This relies on the UMI tagging to put the UMI correctly in th XU tag
Then outputs a re-indexed BAM file for downstream operations 
*/
process deduplication {
    label 'high_memory'
    label 'mid_cpu'
    tag "deduplication for $sample_id"
    publishDir "${params.outDir}/dedup", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bams)
    
    output:
    tuple val(sample_id), path("Aligned.sortedByCoord.deduplicated.out.bam*"), emit: bams
    path("dedup.log"), emit: report_dedup
    
    script:
    (bam, bai) = bams
    """ 
    mkdir -p ./deduplicated
    umi_tools dedup -I $bam \
    --output-stats=./deduplicated \
    --method unique \
    --log ./dedup.log \
    --extract-umi-method=tag \
    --umi-tag=XU \
    > ./Aligned.sortedByCoord.deduplicated.out.bam
    sambamba index -t $task.cpus ./Aligned.sortedByCoord.deduplicated.out.bam
    printf "unique_input_reads: " >> ./dedup.log; samtools view $bam | cut -f1 | sort -u | wc -l >> ./dedup.log
    printf "unique_output_reads: " >> ./dedup.log; samtools view ./Aligned.sortedByCoord.deduplicated.out.bam | cut -f1 | sort -u | wc -l >> ./dedup.log                   
    """
}
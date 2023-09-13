/*
Get sequence coordinates from BAM file for chromosome + ERCC mapping if present
Move the UMI from the alignment to the tag to be used in deduplication 
Merge the tags with the BAM file and index it for deduplication 
*/
process umiTagging {
    label 'mid_cpu'
    label 'mid_memory'
    tag "umiTagging on $sample_id"
    publishDir "${params.outDir}/umiTagging", mode: 'copy'
    
    input:
    tuple val(sample_id) , path(bams)
    
    output:
    tuple val(sample_id), path("Aligned.sortedByCoord.tagged.bam*"), emit: dedup_in_ch
    
    script:
    (bam, bai) = bams
    """
    samtools idxstats $bam | cut -f 1 | grep -E 'chr|ERCC-*' > ./Aligned.sortedByCoord.idxstats.txt
    mkdir -p ./tmp/
    python3 /opt/biorad/src/tagBamFile.py $bam ./Aligned.sortedByCoord.idxstats.txt ./tmp/ $task.cpus
    sambamba merge -t $task.cpus ./Aligned.sortedByCoord.tagged.bam \$(find ./tmp/ | grep -E 'chr|ERCC-*')
    sambamba index -t $task.cpus ./Aligned.sortedByCoord.tagged.bam
    rm -r ./tmp/
    """
}
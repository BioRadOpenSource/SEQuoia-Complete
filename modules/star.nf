/*
Align the fastq reads to the specified refrence genome 
Index output BAM file to better down stream processing 
*/
process starAlign {
    label 'high_memory'
    label 'mid_cpu'
    tag "starAlign on $sample_id"
    publishDir "${params.outDir}/star", mode: 'copy'

    input:
    tuple val(sample_id), path(input)
    path(genomeDirPath)
    path(sjdbGTFFile)

    output:
    tuple val(sample_id), path("Aligned.sortedByCoord.out.bam*"), emit: starBam_ch
    path("Unmapped.out.mate*")
    path("Log.final.out"), emit: report_star //also for meta star

    script:
    """
    STAR \
        --readFilesIn $input \
        --readFilesCommand zcat \
        --genomeDir $genomeDirPath \
        --runThreadN $task.cpus \
        --sjdbGTFfile $sjdbGTFFile \
        --outFilterMismatchNoverLmax 0.05 \
        --outFilterMatchNmin 15 \
        --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 \
        --outReadsUnmapped Fastx \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ./ > star_log.txt 2>&1
    rm -rf _STARgenome
    sambamba index -t $task.cpus Aligned.sortedByCoord.out.bam
    """ 
}
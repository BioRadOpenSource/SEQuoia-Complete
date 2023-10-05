/*
Run fastqc to determine the quality of the reads
*/
process fastQc {
    tag "FASTQC on $sample_id"
    publishDir "${params.outDir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.{zip,html}"), emit: fastqc_results

    script:
    """
    fastqc $reads
    """
}
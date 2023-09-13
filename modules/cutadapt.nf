/*
Using CutAdapt to remove the trailing A's that might be present
in addition to trimming based on the quality score the user provided
If no trim parameter is true then only basic quality trimming is performed
*/
process cutAdapt {
    label 'mid_cpu'
    tag "cutAdapt on $sample_id"
    publishDir "${params.outDir}/cutAdapt", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path( "trimmed_R1.fastq.gz"), emit: trimmed_ch
    path("trimlog.log"), emit: report_trim

    script:
    read1 = reads[0]
    if (!params.noTrim) {
    """
    cutadapt -u 1 \
             -a "A{10}" \
             -m ${params.minBp} \
             -j $task.cpus \
             -q $params.fivePrimeQualCutoff,$params.threePrimeQualCutoff \
             -o trimmed_R1.fastq.gz \
            $read1 1> trimlog.log
    """
    } else { // No trimming is done
    """
    cutadapt -m ${params.minBp} \
             -j $task.cpus \
             -q $params.fivePrimeQualCutoff,$params.threePrimeQualCutoff \
             -o trimmed_R1.fastq.gz \
            $read1 1> trimlog.log
    """
    }
}
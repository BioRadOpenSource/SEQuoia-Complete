/*
Combine the counts files in to one file
For miRNA skip the header information in the first three lines
Then calculate the nofrmalized counts in both TPM and RPKM
TMP = Transcripts Per Kilobase Million
RPKM = Reads Per Kilobase Million
*/
process calcRPKMTPM {
    label 'low_memory'
    label 'low_cpu'
    tag "calcRPKMTPM on $sample_id"
    publishDir "${params.outDir}/calcRPKMTPM", mode: 'copy'
    
    input:
    tuple val(sample_id), path(mi_counts), path(long_counts)

    output:
    path("gene_counts_rpkmtpm.txt"), emit: rpkm_tpm_ch  //normalize_xls
    

    script:
    """
    cat gene_counts_longRNA <(tail -n +3 gene_counts_miRNA) > ./tmp_counts.txt 
    python3 /opt/biorad/src/calc_rpkm_tpm.py ./tmp_counts.txt ./gene_counts_rpkmtpm.txt
    """

}
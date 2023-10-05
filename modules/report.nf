/*
Generate the HTML and PDF reports for the pipleine 
Takes in the outputs and stages them in a way that R 
can access them with the markdown
*/
process assembleReport {
    label 'low_memory'
    tag "assembleReport"
    publishDir "${params.outDir}/report", mode: 'copy' // TODO: Filter down the outputs since so much stuff will be in this dir

    input:
    
    path(fastqc, stageAs: 'out/fastqc/*')
    path(debarcode , stageAs: 'out/debarcode/*') // optional
    path(trim, stageAs: 'out/cutAdapt/*')
    path(star, stageAs: 'out/star/*')
    path(picard, stageAs: 'out/star/*') // Goes into star for reasons
    path(dedup, stageAs: 'out/umitools/*') // optional
    path(miRNACounts, stageAs: 'out/counts/*')
    path(longRNACounts, stageAs: 'out/counts/*')
    path(annoDirPath)

    output:
    path('htmlReport.html')
    path('pdfReport.pdf')

    script:
    """
    mkdir -p ./tmp
    cp /opt/biorad/src/htmlReport.R ./tmp/htmlReport.R
    cp /opt/biorad/src/pdfReport.R ./tmp/pdfReport.R
    Rscript /opt/biorad/src/generateRmdReport.R \$(readlink -f ./out) \$(readlink -f ./tmp)  \$(readlink -f $annoDirPath)
    cp ./tmp/htmlReport.html ./
    cp ./tmp/pdfReport.pdf ./
    """
    
    
    
}
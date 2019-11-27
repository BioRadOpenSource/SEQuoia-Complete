
def paramsWithUsage = readParamsFromJsonSettings()

// Constants
acceptableGenomes = ["rnor6", "hg38", "mm10"]
allowedSpikes     = ["ercc", "miltenyi"]

// Show help emssage
if (params.help){
    helpMessage(paramsWithUsage)
    exit 0
}

// Validate that genome is in correct set
if ( !acceptableGenomes.contains(params.genome) ) {
    log.error "$params.genome not in acceptable genomes $acceptableGenomes"
    exit 1
}

if (params.spikeType != "NONE" && !allowedSpikes.contains(params.spikeType)) {
    log.error "$params.spikeType is not in $allowedSpikes"
    exit 1
}

// Make all the genome related things file resources
genomeDirPath        = file(params.genomes[params.genome][params.spikeType].genomeDir)
annoDirPath          = file(params.genomes[params.genome][params.spikeType].annoDir)
geneId               = params.genomes[params.genome][params.spikeType].geneId
sjdbGTFFile          = file(params.genomes[params.genome][params.spikeType].sjdbGTFFile)
refFlatFile          = file(params.genomes[params.genome][params.spikeType].refFlatFile)
ribosomalIntervalFile = file(params.genomes[params.genome][params.spikeType].ribosomalIntervalFile)
miRNABedFile         = file(params.genomes[params.genome][params.spikeType].miRNABedFile)
miRNAgtfFile         = file(params.genomes[params.genome][params.spikeType].miRNAgtfFile)
longRNAgtfFile       = file(params.genomes[params.genome][params.spikeType].longRNAgtfFile)
sizesFile            = file(params.genomes[params.genome][params.spikeType].sizesFile)

// Check that inputs are viable / R1 only vs both
if (params.reads == "NOINPUT") {
    log.error "No input reads were supplied"
    exit 1
}


// Create Summary
def summary = [:]
summary['Run Name'] = workflow.runName
summary['Reads'] = params.reads
summary['Genome'] = params.genome
summary['Spike Type'] = params.spikeType
summary['Reference Dir'] = params.genomes_base
summary['Annotations Dir'] = annoDirPath
summary['Skip UMI?'] = params.skipUmi
summary['Min MAPQ To Count'] = params.minMapqToCount
summary['Output Dir'] = params.outDir
summary['Trace Dir'] = params.tracedir
/*summary['Max Cores'] = task.cpus*/
summary['geneId'] = geneId
summary['sjdb GTF File'] = sjdbGTFFile
summary['ref Flat File'] = refFlatFile
summary['Ribosoma Interval File'] = ribosomalIntervalFile
summary['miRNA Bed File'] = miRNABedFile
summary['miRNA GTF File'] = miRNAgtfFile
summary['Long RNA GTF File'] = longRNAgtfFile
summary['Sizes File'] = sizesFile
if(workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
log.info bioradHeader()
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "----------------------------------------------------"

// Create Channels
Channel
    //Todo, make this more robust
    .fromFilePairs( params.reads, size: params.skipUmi ? 1 : 2) //, size: params.skipUmi ? 1 : 2 ) // Assume we always pass in R1 and R2, but if skipumi, only use R1
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf not using R2, pass --skipUmi on CLI" }
    .into { raw_reads_fastqc; raw_reads; raw_reads_validation }


// Begin Processing

if (params.validateInputs) {
    process validateInputs {
        tag "Validation on $sample_id"
        publishDir "${params.outDir}/validation", mode: 'copy'

        input:
        set val(sample_id), file(reads) from raw_reads_validation

        script:
        """
        python3 /opt/biorad/src/validate.py $reads 2>&1 > validation.log
        """
    }
}
process fastQc {
    tag "FASTQC on $sample_id"
    publishDir "${params.outDir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(sample_id), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results, report_fastqc

    script:
    """
    fastqc $reads
    """
}

// Only extract barcodes if umiAware
if (!params.skipUmi) {
    process debarcode{
        label 'mid_cpu'
        tag "debarcode on $sample_id"
        publishDir "${params.outDir}/debarcode", mode: 'copy'

        input:
        set val(sample_id), file(reads) from raw_reads

        output:
        set val(sample_id), file('*_R1.fastq.gz') into debarcoded_ch
        file 'debarcode_stats.txt' into report_debarcode

        script:
        """
        bash /opt/biorad/src/fastq_to_tsv.sh $reads \
            | parallel --pipe python3 /opt/biorad/src/debarcode.py \
            | tee >(awk '/^MM/{bad=bad+1}/^@/{good=good+1}END{print "Total Reads: " good + bad; print "Good Reads: " good; print "Bad Reads: " bad}' > debarcode_stats.txt) \
            | grep -ve '^MM' \
            | bash /opt/biorad/src/tsv_to_fastq.sh ${sample_id}_debarcoded_R1.fastq.gz ${sample_id}_debarcoded_R2.fastq.gz compress
        """
    }
} else {
    /*debarcoded_ch = raw_reads.map{ [it[0], it[1][0]]} // Reduce inputs to just read one in the event that there were 2*/
    debarcoded_ch = raw_reads
    report_debarcode = Channel.empty()
}

// There should only be a R1 at 
process cutAdapt {
    label 'mid_cpu'
    tag "cutAdapt on $name"
    publishDir "${params.outDir}/cutAdapt", mode: 'copy'
    
    input:
    set val(name), file(read1) from debarcoded_ch

    output:
    set val(name), file( "trimmed_R1.fastq.gz") into trimmed_ch
    file "trimlog.log" into report_trim

    script:
    if (!params.noTrim) {
    """
    cutadapt -u 1 -a "A{10}" -m ${params.minBp} -j $task.cpus \
             -q $params.fivePrimeQualCutoff,$params.threePrimeQualCutoff \
             -o trimmed_R1.fastq.gz $read1 1> trimlog.log
    """
    } else { // No trimming is done
    """
    cutadapt -m ${params.minBp} -j $task.cpus \
             -q $params.fivePrimeQualCutoff,$params.threePrimeQualCutoff -o trimmed_R1.fastq.gz $read1 1> trimlog.log
    """
    }
}

process starAlign {
    label 'high_memory'
    label 'mid_cpu'
    tag "starAlign on $name"
    publishDir "${params.outDir}/star", mode: 'copy'


    input:
    set val(name), file(input) from trimmed_ch
    file genomeDirPath
    file sjdbGTFFile

    output:
    set val(name), file("Aligned.sortedByCoord.out.bam*") into umiTagging_ch, picardBam_ch
    file "Log.final.out" into report_star

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
        # --outSAMmultNmax 1 \
        # ^ do not print mms, just report in NH tag
        # --outMultimapperOrder Random \
        # --runRNGseed 1234 \
        --outFileNamePrefix ./ > star_log.txt 2>&1
    rm -rf _STARgenome
    sambamba index -t $task.cpus Aligned.sortedByCoord.out.bam
    """ 
}

process picardAlignSummary {
    label 'low_memory'
    tag "picardAlignSummary on $name"
    publishDir "${params.outDir}/picardAlignSummary", mode: 'copy'

    input:
    set val(name), file(bams) from picardBam_ch
    file refFlatFile
    file ribosomalIntervalFile

    output:
    file 'rna_metrics.txt' into report_picard

    script:
    (bam, bai) = bams
    strand = params.reverseStrand ? "SECOND_READ_TRANSCRIPTION_STRAND" : "FIRST_READ_TRANSCRIPTION_STRAND"
    """
    java -jar /opt/picard/picard.jar CollectRnaSeqMetrics I=$bam \
    O=rna_metrics.txt \
    REF_FLAT=$refFlatFile \
    STRAND=$strand \
    RIBOSOMAL_INTERVALS=$ribosomalIntervalFile
    """
}

if (!params.skipUmi) {
    process umiTagging {
        label 'mid_cpu'
        label 'mid_memory'
        tag "umiTagging on $name"
        publishDir "${params.outDir}/umiTagging", mode: 'copy'
        input:
        set val(name) , file(bams) from umiTagging_ch

        output:
        set val(name), file("Aligned.sortedByCoord.tagged.bam*") into dedup_in_ch
        
        script:
        (bam, bai) = bams
        """
        samtools idxstats $bam | cut -f 1 | grep chr > ./Aligned.sortedByCoord.idxstats.txt
        # samtools idxstats $bam | cut -f 1 | uniq > ./Aligned.sortedByCoord.idxstats.txt
        mkdir -p ./tmp/
        python3 /opt/biorad/src/tagBamFile.py $bam ./Aligned.sortedByCoord.idxstats.txt ./tmp/ $task.cpus
        sambamba merge -t $task.cpus ./Aligned.sortedByCoord.tagged.bam \$(find ./tmp/ | grep chr)
        #sambamba merge -t $task.cpus ./Aligned.sortedByCoord.tagged.bam \$(find ./tmp/ | grep .bam)
        sambamba index -t $task.cpus ./Aligned.sortedByCoord.tagged.bam
        rm -r ./tmp/
        """

    }
    process deduplication {
        label 'high_memory'
        label 'mid_cpu'
        tag "dedup on $name"
        publishDir "${params.outDir}/dedup", mode: 'copy'

        input:
        set val(name), file(bams) from dedup_in_ch

        output:
        set val(name), file("Aligned.sortedByCoord.deduplicated.out.bam*") into splitBamMi_ch, splitBamLong_ch
        file 'dedup.log' into report_dedup

        script:
        (bam, bai) = bams
        """
        mkdir -p ./deduplicated
        umi_tools dedup -I $bam --output-stats=./deduplicated \
        --method unique --log ./dedup.log \
        --extract-umi-method=tag --umi-tag=XU \
        > ./Aligned.sortedByCoord.deduplicated.out.bam
        sambamba index -t $task.cpus ./Aligned.sortedByCoord.deduplicated.out.bam
        printf "unique_input_reads: " >> ./dedup.log; samtools view $bam | cut -f1 | sort -u | wc -l >> ./dedup.log
        printf "unique_output_reads: " >> ./dedup.log; samtools view ./Aligned.sortedByCoord.deduplicated.out.bam | cut -f1 | sort -u | wc -l >> ./dedup.log
        """
    }
} else {
    umiTagging_ch.into { splitBamMi_ch; splitBamLong_ch }
    report_dedup = Channel.empty()
}

process splitBamMi {
    label 'low_memory'
    tag "splitBamMi on $name"
    publishDir "${params.outDir}/splitBamMi", mode: 'copy'
    input:
    set val(name), file(bams) from splitBamMi_ch
    file miRNABedFile 

    output:
    set val(name), file('out.miRNAs.bam') into mirna_bam_ch, bigwig_mirna_bam_ch

    script:
    (bam, bai) = bams
    """
    # Split bam_to_split into long and short RNAs
    bedtools intersect -a $bam -b $miRNABedFile -f 1 -s > out.miRNAs.bam
    """
}

process splitBamLong {
    label 'low_memory'
    tag "splitBamLong on $name"
    publishDir "${params.outDir}/splitBamLong", mode: 'copy'
    input:
    set val(name), file(bams) from splitBamLong_ch
    file miRNABedFile

    output:
    set val(name), file('out.longRNAs.bam') into longrna_bam_ch, bigwig_long_rna_bam_ch

    script:
    (bam, bai) = bams
    """
    bedtools intersect -a $bam -b $miRNABedFile -f 1 -v -s > out.longRNAs.bam
    """
}

process countMicroRNA {
    label 'mid_cpu'
    label 'low_memory'
    tag "countMicroRNA on $name"
    publishDir "${params.outDir}/microRNACounts", mode: 'copy'

    input:
    set val(name), file(bam) from mirna_bam_ch
    file miRNAgtfFile

    output:
    set val(name), file('gene_counts_miRNA') into mi_counts_ch
    file 'gene_counts_miRNA*' into report_miRNACounts

    script:
    strand = params.reverseStrand ? "-s 2" : "-s 1" 
    """
    featureCounts -T $task.cpus --primary -O -M -t exon -g product $strand -Q $params.minMapqToCount \
    -a $miRNAgtfFile \
    -o ./gene_counts_miRNA \
    -R BAM ./out.miRNAs.bam
    """
}

process countLongRNA {
    label 'mid_cpu'
    label 'low_memory'
    tag "countLongRNA on $name"
    publishDir "${params.outDir}/longRNACounts", mode: 'copy'

    input:
    set val(name), file(bam) from longrna_bam_ch
    file longRNAgtfFile
    output:
    set val(name), file('gene_counts_longRNA') into long_counts_ch
    file "gene_counts_longRNA*" into report_longRNACounts

    script:
    strand = params.reverseStrand ? "-s 2" : "-s 1" 
    """
    featureCounts -T $task.cpus --primary -M -t exon -g $geneId $strand -Q $params.minMapqToCount \
    -a $longRNAgtfFile \
    -o ./gene_counts_longRNA \
    -R BAM ./out.longRNAs.bam
    """
}

process calcRPMKTPM {
    label 'low_memory'
    label 'low_cpu'
    tag "calcRPMKTPM on $l_name"
    publishDir "${params.outDir}/calcRPMKTPM", mode: 'copy'
    input:
    set val(l_name), file(long_counts) from long_counts_ch
    set val(m_name), file(mi_counts) from mi_counts_ch

    output:
    file 'gene_counts_rpkmtpm.txt' into rpkm_tpm_ch

    script:
    """
    cat $long_counts <(tail -n +3 $mi_counts) > ./tmp_counts.txt
    python3 /opt/biorad/src/calc_rpkm_tpm.py ./tmp_counts.txt ./gene_counts_rpkmtpm.txt
    """

}

/*process genBigWig {*/
    /*label 'mid_memory'*/
    /*label 'low_cpu'*/
    /*tag "genBigWig from $l_name"*/
    /*publishDir "${params.outDir}/bigwig", mode: 'copy'*/
    /*input:*/
    /*set val(l_name), file(long_bam) from bigwig_long_rna_bam_ch*/
    /*set val(m_name), file(mi_bam) from bigwig_mirna_bam_ch*/
    /*file sizesFile */

    /*script:*/
    /*//TODO: split this into two processes so long and mi can run in parallel*/
    /*//      or read from one channel and pass in the output name or something*/
    /*//      use a .filter { it[1].size() > 0 } to prevent 0 size stuff running*/
    /*"""*/
    /*if [ \$(samtools view $long_bam | head -n1 | wc -l) -eq 1 ]; then*/
        /*bedtools genomecov -bg -ibam $long_bam | sort -k1,1 -k2,2n > ./longRNAs.bedGraph*/
        /*bedGraphToBigWig ./longRNAs.bedGraph $sizesFile ./longRNAs.bigWig*/
    /*else */
        /*echo $long_bam " is empty"*/
    /*fi*/

    /*if [ \$(samtools view $mi_bam | head -n1 | wc -l) -eq 1 ]; then */
        /*bedtools genomecov -bg -ibam $mi_bam | sort -k1,1 -k2,2n > ./miRNAs.bedGraph*/
        /*bedGraphToBigWig ./miRNAs.bedGraph $sizesFile ./miRNAs.bigWig*/
    /*else*/
        /*echo $mi_bam " is empty"*/
    /*fi*/
    /*"""*/
/*}*/

process assembleReport {
    label 'low_memory'
    tag "assembleReport"
    publishDir "${params.outDir}/report", mode: 'copy' // TODO: Filter down the outputs since so much stuff will be in this dir

    input:
    file annoDirPath
    file(fastqc: 'out/fastqc/*') from report_fastqc.collect()
    file('out/debarcode/*') from report_debarcode.collect().ifEmpty([]) // optional
    file('out/cutAdapt/*') from report_trim.collect()
    file('out/star/*') from report_star.collect() 
    file('out/star/*') from report_picard.collect() // Goes into star for reasons
    file('out/umitools/*') from report_dedup.collect().ifEmpty([]) // optional
    file('out/counts/*') from report_miRNACounts.collect()
    file('out/counts/*') from report_longRNACounts.collect()

    output:
    file 'htmlReport.html'
    file 'pdfReport.pdf'

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

/* Helper Functions */
def readParamsFromJsonSettings() {
    List paramsWithUsage
    try {
        paramsWithUsage = tryReadParamsFromJsonSettings()
    } catch (Exception e) {
        println "Could not rea parameters settings from json. $e"
        pramsWithUsage = Collections.emptyMap()
    }
    return paramsWithUsage
}

def tryReadParamsFromJsonSettings() throws Exception {
    def paramsContent = new File(config.params_description.path).text
    def paramsWithUsage = new groovy.json.JsonSlurper().parseText(paramsContent)
    return paramsWithUsage.get('parameters')
}

String prettyFormatParamsWithPaddingAndIndent(List paramGroup, Integer padding=2, Integer indent=4) {
    def fields = ["name", "usage", "type", "default_value", "pattern", "choices"]
    def maxFields = fields.collectEntries { String field ->
        [(field): paramGroup.collect {
            def val = it.get(field)
            val ? val.toString().size() : 1
        }.max()]
    }
    def formatter = {param -> 
        sprintf("%${indent}s%-${maxFields.name}s (%-${maxFields.type}s) %-${maxFields.default_value}s %-${maxFields.usage}s %-${maxFields.choices}s %-${maxFields.pattern}s\n", "",
                                param.name ?: "", param.type ?: "",  param.default_value ?: "", param.usage ?: "", param.choices ?: "", param.pattern ?: "")
    }
    def requiredParamsFormattedList = paramGroup.sort { it.name }.findAll { it.required }.collect { Map param -> formatter(param) }
    def optionalParamsFormattedList = paramGroup.sort { it.name }.findAll { !it.required }.collect { Map param -> formatter(param) }
    return String.format("REQUIRED:\n%s\nOPTIONAL:\n%s\n", requiredParamsFormattedList.join(), optionalParamsFormattedList.join())
}

def helpMessage(paramsWithUsage) {
    def helpMessage = String.format(
    """\
    %s
    
    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run . --reads './tests/*_R{1,2}.fastq.gz' --genome hg38 --outDir /data/out --skipUmi --genomes_base /mnt/genome-annotations

    Args:

    %s
    """.stripIndent(), bioradHeader(), prettyFormatParamsWithPaddingAndIndent(paramsWithUsage, 8, 4))
    log.info helpMessage
}

def bioradHeader() {
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    return """
    ${c_reset}
    ${c_green}__________.__                __________             .___
    ${c_green}\\______   \\__| ____          \\______   \\_____     __| _/
    ${c_green}|    |  _/  |/  _ \\   ______ |       _/\\__  \\   / __ | 
    ${c_green}|    |   \\  (  <_> ) /_____/ |    |   \\ / __ \\_/ /_/ | 
    ${c_green}|______  /__|\\____/          |____|_  /(____  /\\____ | 
    ${c_green}       \\/                           \\/      \\/      \\/ 
    ${c_reset}
    """.stripIndent()
}

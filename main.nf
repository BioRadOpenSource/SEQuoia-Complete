nextflow.enable.dsl = 2

def paramsWithUsage = readParamsFromJsonSettings()

// Constants
acceptableGenomes = ['rnor6', 'hg38', 'mm10']
allowedSpikes     = ['ercc', 'miltenyi']

// Show help emssage
if (params.help) {
    helpMessage(paramsWithUsage)
    exit 0
}

// Validate that genome is in correct set
if (!acceptableGenomes.contains(params.genome)) {
    log.error "$params.genome not in acceptable genomes $acceptableGenomes"
    exit 1
}

if (params.spikeType != 'NONE' && !allowedSpikes.contains(params.spikeType)) {
    log.error "$params.spikeType is not in $allowedSpikes"
    exit 1
}

//include modules containing processes
include { validateInputs }                    from './modules/validate_inputs'
include { fastQc }                            from './modules/fastqc'
include { debarcode }                         from './modules/debarcode'
include { cutAdapt }                          from './modules/cutadapt'
include { starAlign }                         from './modules/star'
include { picardAlignSummary }                from './modules/picard'
include { umiTagging }                        from './modules/umi_tagging'
include { deduplication }                     from './modules/deduplication'
include { splitBamLong; splitBamMi }          from './modules/split_bam'
include { countLongRNA; countMicroRNA }       from './modules/count_rna'
include { calcRPKMTPM }                       from './modules/rpkm'
include { assembleReport }                    from './modules/report'

// Make all the genome related things file resources
genomeDirPath         = file(params.genomes[params.genome][params.spikeType].genomeDir)
annoDirPath           = file(params.genomes[params.genome][params.spikeType].annoDir)
geneId                = params.genomes[params.genome][params.spikeType].geneId
sjdbGTFFile           = file(params.genomes[params.genome][params.spikeType].sjdbGTFFile)
refFlatFile           = file(params.genomes[params.genome][params.spikeType].refFlatFile)
ribosomalIntervalFile = file(params.genomes[params.genome][params.spikeType].ribosomalIntervalFile)
miRNABedFile          = file(params.genomes[params.genome][params.spikeType].miRNABedFile)
miRNAgtfFile          = file(params.genomes[params.genome][params.spikeType].miRNAgtfFile)
longRNAgtfFile        = file(params.genomes[params.genome][params.spikeType].longRNAgtfFile)
sizesFile             = file(params.genomes[params.genome][params.spikeType].sizesFile)

// Check that inputs are viable / R1 only vs both
if (params.reads == 'NOINPUT') {
    log.error 'No input reads were supplied'
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
summary['Max Cores'] = params.max_cpus
summary['Max Ram'] = params.max_memory
summary['geneId'] = geneId
summary['sjdb GTF File'] = sjdbGTFFile
summary['ref Flat File'] = refFlatFile
summary['Ribosoma Interval File'] = ribosomalIntervalFile
summary['miRNA Bed File'] = miRNABedFile
summary['miRNA GTF File'] = miRNAgtfFile
summary['Long RNA GTF File'] = longRNAgtfFile
summary['Sizes File'] = sizesFile
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile == 'awsbatch') {
    summary['AWS Region']    = params.awsregion
    summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
log.info bioradHeader()
log.info summary.collect { k, v -> "${k.padRight(22)}: $v" }.join('\n')
log.info '----------------------------------------------------'

workflow {
    // Create Channels
    raw_reads = Channel
        .fromFilePairs(params.reads, size: params.skipUmi ? 1 : 2) //, size: params.skipUmi ? 1 : 2) // Assume we always pass in R1 and R2, but if skipumi, only use R1
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf not using R2, pass --skipUmi on CLI" }

    // Begin Processing

    if (params.validateInputs) {
        validateInputs(raw_reads)
    }
    fastQc(raw_reads)

    // Only extract barcodes if umiAware
    if (!params.skipUmi) {
        debarcode(raw_reads)
        debarcoded_ch = debarcode.out.reads
        report_debarcode = debarcode.out.report
    } else {
        debarcoded_ch = raw_reads
        report_debarcode = Channel.empty()
    }

    // There should only be a R1 at this point
    cutAdapt(debarcoded_ch)
    starAlign(cutAdapt.out.trimmed_ch, genomeDirPath, sjdbGTFFile)
    picardAlignSummary(starAlign.out.starBam_ch, refFlatFile, ribosomalIntervalFile)

    if (!params.skipUmi) {
        umiTagging(starAlign.out.starBam_ch)
        deduplication(umiTagging.out.dedup_in_ch)
        splitBamMi_ch = deduplication.out.bams
        splitBamLong_ch = deduplication.out.bams
    } else {
        splitBamMi_ch = starAlign.out.starBam_ch
        splitBamLong_ch = starAlign.out.starBam_ch
        report_dedup = Channel.empty()
    }

    splitBamMi(splitBamMi_ch, miRNABedFile)
    splitBamLong(splitBamLong_ch, miRNABedFile)

    //count RNAs
    countMicroRNA(splitBamMi.out.mirna_bam_ch, miRNAgtfFile, geneId)
    countLongRNA(splitBamLong.out.longrna_bam_ch, longRNAgtfFile, geneId)

    //RPKM and TMP
    calcRPKMTPM(countLongRNA.out.counts_ch_long.join(countMicroRNA.out.counts_ch_mi))

    //generate report
    assembleReport(
        fastQc.out.fastqc_results,
        debarcode.out.report.ifEmpty([]),
        cutAdapt.out.report_trim,
        starAlign.out.report_star,
        picardAlignSummary.out.report_picard,
        deduplication.out.report_dedup.ifEmpty([]),
        countMicroRNA.out.report_miRNACounts,
        countLongRNA.out.report_longRNACounts,
        annoDirPath
        )
}//end workflow

/* Helper Functions */
def readParamsFromJsonSettings() {
    List paramsWithUsage
    try {
        paramsWithUsage = tryReadParamsFromJsonSettings()
    } catch (Exception e) {
        println "Could not read parameters settings from json. $e"
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
    def fields = ['name', 'usage', 'type', 'default_value', 'pattern', 'choices']
    def maxFields = fields.collectEntries { String field ->
        [(field): paramGroup.collect {
            def val = it.get(field)
            val ? val.toString().size() : 1
        }.max()]
    }
    def formatter = { param ->
        sprintf("%${indent}s%-${maxFields.name}s (%-${maxFields.type}s) %-${maxFields.default_value}s %-${maxFields.usage}s %-${maxFields.choices}s %-${maxFields.pattern}s\n", '',
                                param.name ?: '', param.type ?: '',  param.default_value ?: '', param.usage ?: '', param.choices ?: '', param.pattern ?: '')
    }
    def requiredParamsFormattedList = paramGroup.sort { it.name }.findAll { it.required }.collect { Map param -> formatter(param) }
    def optionalParamsFormattedList = paramGroup.sort { it.name }.findAll { !it.required }.collect { Map param -> formatter(param) }
    return String.format('REQUIRED:\n%s\nOPTIONAL:\n%s\n', requiredParamsFormattedList.join(), optionalParamsFormattedList.join())
}

def helpMessage(paramsWithUsage) {
    def helpMessage = String.format(
    """\
    %s

    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run . --reads './tests/*_R{1,2}.fastq.gz' --genome hg38 --outDir ~/data/out --skipUmi --genomes_base ~/genome-annotations

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

#!/usr/bin/env nextflow
/*

========================================================================================
               nf-core/ E X O S E Q    B E S T    P R A C T I C E
========================================================================================
 #### Homepage / Documentation
 https://github.com/scilifelab/NGI-ExoSeq
 #### Authors
 Senthilkumar Panneerselvam @senthil10 <senthilkumar.panneerselvam@scilifelab.se>
 Phil Ewels @ewels <phil.ewels@scilifelab.se>
 Alex Peltzer @alex_peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
 Marie Gauder <marie.gauder@student.uni-tuebingen.de>
 Marc Hoeppner <m.hoeppner@ikmb.uni-kiel.de>

 Some code parts were lent from other NGI-Pipelines (e.g. CAW), specifically the error
 handling, logging messages. Thanks for that @CAW guys.
----------------------------------------------------------------------------------------
Developed based on GATK's best practise, takes set of FASTQ files and performs:
 - alignment (BWA)
 - recalibration (GATK)
 - variant calling (GATK)
 - variant evaluation (SnpEff)
*/

// Help message
helpMessage = """
===============================================================================
nf-core/ExoSeq : Exome/Targeted sequence capture best practice analysis v${params.version}
===============================================================================

Usage: nextflow nf-core/ExoSeq --reads '*_R{1,2}.fastq.gz' --genome GRCh37

This is a typical usage where the required parameters (with no defaults) were
given. The available paramaters are listed below based on category

Required parameters:
--reads	                       Absolute path to project directory (excludes the --samples option)
--samples		       A CSV formatted file containing sample information (excludes the --reads option)
--genome                       Name of iGenomes reference

Output:
--outdir                       Path where the results to be saved [Default: './results']
--bam 			       Generate the recalibrated alignment file in BAM format (standard: CRAM)

Kit files:
--kit                          Kit used to prep samples [Default: 'agilent_v5']
--bait                         Absolute path to bait file
--target                       Absolute path to target file

Genome/Variation files:
--dbsnp                        Absolute path to dbsnp file
--thousandg                    Absolute path to 1000G file
--mills                        Absolute path to Mills file
--omni                         Absolute path to Omni file
--gfasta                       Absolute path to genome fasta file

Other options:
--project                      Project name, can be passed to scheduler for accounting purposes (e.g. Uppmax)

For more detailed information regarding the parameters and usage refer to package
documentation at https://github.com/nf-core/ExoSeq
"""

// Variables and defaults
params.name = false
params.help = false
params.reads = false
params.samples = false
params.genome = false
params.run_id = false
params.aligner = 'bwa' //Default, but stay tuned for later ;-) 
params.saveReference = true

// Output configuration
params.outdir = './results'
params.saveAlignedIntermediates = false
params.saveIntermediateVariants = false
params.bam = false // default alignment format is CRAM

// Clipping options
params.notrim = false
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Kit options
params.kit = 'agilent_v5'
params.bait = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].bait ?: false : false
params.target = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].target ?: false : false
params.dbsnp = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].dbsnp ?: false : false
params.thousandg = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].thousandg ?: false : false
params.mills = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].mills ?: false : false
params.omni = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].omni ?: false : false
params.gfasta = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].gfasta ?: false : false

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// Check blocks for certain required parameters, to see they are given and exist
//if (!params.reads || !params.genome){
//    exit 1, "Parameters '--reads' and '--genome' are required to run the pipeline"
//
//}

if (!params.kitFiles[ params.kit ] && ['bait', 'target'].count{ params[it] } != 2){
    exit 1, "Kit '${params.kit}' is not available in pre-defined config, so " +
            "provide all kit specific files with option '--bait' and '--target' "
}
 if (!params.metaFiles[ params.genome ] && ['gfasta', 'dbsnp', 'thousandg', 'mills', 'omni'].count{ params[it] } != 5){
     exit 1, "Genome '${params.genome}' is not available in pre-defined config, so you need to provide all genome specific " +
             "files with options '--gfasta', '--dbsnp', '--thousandg', '--mills' and '--omni' "
 }

// Create a channel for input files
// Can be either a CSV file or a folder/regexp pointing to fastq files

if (params.samples) {
	Channel.from(file(params.samples))
       .splitCsv(sep: ';', header: true)
       .set {  read_files_trimming }
} else if (params.reads) {
        extractFastqFromDir(params.reads)
        .set { read_files_trimming }
} else {
	exit 1, "No samples were defined, see --help"
}

// Create a summary for the logfile
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Genome']       = params.genome
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
summary['Aligner'] = "BWA Mem"
summary['Fasta Ref']    = params.gfasta
summary['Save Intermediate Aligned Files'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Save Intermediate Variant Files'] = params.saveIntermediateVariants ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


try {
    if( ! workflow.nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
        }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * STEP 1 - trim with trim galore
*/

if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = []
    trimgalore_logs = []
} else {
    process trim_galore {
        tag "${patientID}|${sampleID}|${libraryID}"
        publishDir "${params.outdir}/trim_galore", mode: 'copy'

        input:
        set val(patientID), val(sampleID), val(libraryID), val(rgID),R1,R2 from read_files_trimming

        output:
        set val(patientID), val(sampleID), val(libraryID), val(rgID),file("*_val_1*.fq.gz"),file("*_val_2*.fq.gz") into trimmed_reads
        file '*trimming_report.txt' into trimgalore_results, trimgalore_logs
   	file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
        
        script:
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''

        """
            trim_galore --fastqc --paired --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $R1 $R2
        """
    }
}

/*
 * STEP 2 - Map with BWA Mem
*/

process bwamem {
    tag "${patientID}|${sampleID}"

    input:
    set val(patientID),val(sampleID),val(libraryID),val(rgID), file(R1),file(R2) from trimmed_reads

    output:
    set val(patientID),val(sampleID), file(bam) into samples_sorted_bam
    file '.command.log' into bwa_stdout

    script:
    bam = sampleID + "_" + libraryID + "_" + rgID  + "_bwa.bam"
    avail_mem = task.memory ? "-m ${task.memory.toMega().intdiv(task.cpus)}M" : ''
    rg="\'@RG\\tID:${rgID}\\tLB:${libraryID}\\tSM:${sampleID}\\tDS:${params.gfasta}\\tPL:illumina\'"
        
    """
    bwa mem \\
    -M \\
    -R $rg \\
    -t ${task.cpus} \\
    ${params.gfasta} \\
    $R1 $R2 | samtools ${avail_mem} sort -O bam - > ${bam}
    """
}

/*
* Merge BAM files per sample
*/

bwa_grouped_by_sample = samples_sorted_bam.groupTuple(by: [0,1])

process merge_bam_by_sample {

        tag "${patientID}|${sampleID}"
	
	input:
        set val(patientID),val(sampleID), file(aligned_bam_list) from bwa_grouped_by_sample

	output:
	set val(patientID),val(sampleID),file(merged_bam) into merged_bam_by_sample

	script:
	merged_bam = sampleID + "_merged.bam"

	"""
		picard MergeSamFiles \
			INPUT=${aligned_bam_list.join(' INPUT=')} \
			OUTPUT=${merged_bam} \
			CREATE_INDEX=false \
			CREATE_MD5_FILE=false \
			SORT_ORDER=coordinate
	"""
}

/*
*  STEP 3 - Mark PCR duplicates in sorted BAM file
*/ 

process markDuplicates {
    tag "${patientID}|${sampleID}"
    publishDir "${params.outdir}/Picard_Markduplicates/metrics", mode: 'copy',
        saveAs: { filename -> filename.indexOf("_dedup_metrics.txt") > 0 ? filename : null }

    input:
    set val(patientID), val(sampleID), file(merged_bam) from merged_bam_by_sample

    output:
    set val(patientID),val(sampleID),file(dedup_bam),file(dedup_bai) into samples_markdup_bam
    file(metrics) into markdup_results
    file '.command.log' into markDuplicates_stdout

    script:
    dedup_bam = patientID + "_" + sampleID + "_dedup.bam"
    dedup_bai = dedup_bam + ".bai"
    metrics = patientID + "_" + sampleID + "_dedup_metrics.txt"

    """
        mkdir `pwd`/tmp
        gatk MarkDuplicates \\
        --INPUT ${sorted_bam} \\
        --OUTPUT ${dedup_bam} \\
        --METRICS_FILE ${metrics} \\
        --REMOVE_DUPLICATES false \\
        --CREATE_INDEX true \\
        --TMP_DIR=`pwd`/tmp \\
        --java-options -Xmx${task.memory.toGiga()}G
    """
}

/*
 * Step 5 - Recalibrate BAM file with known variants and BaseRecalibrator
 * 
*/ 

process recal_bam_files {
    tag "${patientID}|${sampleID}"
    publishDir "${params.outdir}/GATK_Recalibration", mode: 'copy'

    input:
    set val(patientID),val(sampleID), file(markdup_bam), file(markdup_bam_ind) from samples_markdup_bam

    output:
    set val(patientID),val(sampleID), val(markdup_bam), file(recal_table) into samples_recal_reports
    file '.command.log' into gatk_stdout
    file '.command.log' into gatk_base_recalibration_results

    script:
    recal_table = patientID + "_" + sampleID + "_table.recal"
    """
    gatk-launch BaseRecalibrator \\
        -R ${params.gfasta} \\
        -I ${markdup_bam} \\
        -O ${recal_table} \\
        -L ${params.target} \\
        --known-sites ${params.dbsnp} \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}

process applyBQSR {
    tag "${patientID}|${sampleID}"
    publishDir "${params.outdir}/GATK_ApplyBQSR", mode: 'copy'

    input:
    set val(patientID),val(sampleID),file(markdup_bam), file(recal_table) from samples_recal_reports

    output:
    set val(patientID), val(sampleID), file(clean_bam), file(clean_bai) into bam_vcall, bam_metrics, bam_for_multiple_metrics, bam_for_hs_metrics

    script:
    clean_bam = patientID + "_" + sampleID + "clean.bam"
    clean_bai = clean_bam + ".bai"
    """
    gatk-launch ApplyBQSR \\
        -R ${params.gfasta} \\
        -I ${markdup_bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${clean_bam} \\
        -L ${params.target} \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """

}

/*
 * Generate relevant statistics
 *
*/

process picard_multiple_metrics {
	tag "${patientID}|${sampleID}"
	publishDir "${params.outdir}/Processing/Picard_Metrics", mode: 'copy'
 
	input:
	set val(patientID), val(sampleID), file(bam), file(bai) from bam_for_multiple_metrics

	output:
	file("${prefix}*") into CollectMultipleMetricsOutput mode flatten

	script:       
	prefix = patientID + "_" + sampleID

	"""
		picard CollectMultipleMetrics \
		PROGRAM=MeanQualityByCycle \
		PROGRAM=QualityScoreDistribution \
		PROGRAM=CollectAlignmentSummaryMetrics \
		PROGRAM=CollectInsertSizeMetrics\
       	        PROGRAM=CollectSequencingArtifactMetrics \
                PROGRAM=CollectQualityYieldMetrics \
	        PROGRAM=CollectGcBiasMetrics \
		PROGRAM=CollectBaseDistributionByCycle \
		INPUT=${bam} \
		REFERENCE_SEQUENCE=${params.gfasta} \
		DB_SNP=${params.dbsnp} \
		INTERVALS=${params.bait} \
		ASSUME_SORTED=true \
		QUIET=true \
		OUTPUT=${prefix} \
		TMP_DIR=tmp
	"""
}	

process picard_hc_metrics {

    tag "${patientID}|${sampleID}"
    publishDir "${params.outdir}/Picard_Metrics", mode: 'copy'

    input:
    set val(patientID),val(sampleID), file(bam), file(bai) from bam_for_hs_metrics

    output:
    file(outfile) into HybridCaptureMetricsOutput mode flatten

    script:
    outfile = patientID + "_" + sampleID  + ".hybrid_selection_metrics.txt"

    """
        picard CollectHsMetrics \
               INPUT=${bam} \
               OUTPUT=${outfile} \
               TARGET_INTERVALS=${params.target} \
               BAIT_INTERVALS=${params.bait} \
               REFERENCE_SEQUENCE=${params.gfasta} \
               TMP_DIR=tmp
        """
}

/*
 * Step 8 - Call Variants with HaplotypeCaller in GVCF mode (differentiate between exome and whole genome data here)
 * 
*/ 

process variantCall {

    tag "${patientID}|${sampleID}"
    publishDir "${params.outdir}/GATK_VariantCalling/", mode: 'copy'

    input:
    set val(patientID),val(sampleID), file(realign_bam), file(realign_bam_ind) from bam_vcall

    output:
    set val(patientID),val(sampleID), file(gvcf), file(gvcf_tbi) into raw_variants

    script:
    gvcf = patientID + "_" + sampleID + ".g.vcf.gz"
    gvcf_tbi = gvcf + ".tbi"

    """
    gatk-launch HaplotypeCaller \\
        -I ${realign_bam} \\
        -R ${params.gfasta} \\
        -O ${gvcf} \\
        -ERC GVCF \\
        -L ${params.target} \\
        --create-output-variant-index \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp ${params.dbsnp} \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}

/*
* Step 9 - Generate a YAML file for software versions in the pipeline
* This is then parsed by MultiQC and the report feature to produce a final report with the software Versions in the pipeline.
*/ 

process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo "$params.version" &> v_nfcore_exoseq.txt
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    samtools --version &> v_samtools.txt
    bwa &> v_bwa.txt 2>&1 || true
    gatk-launch BaseRecalibrator --version &> v_gatk.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
* Step 10 - Generate MultiQC config file
*
*/ 

process GenerateMultiQCconfig {

  publishDir "${params.outdir}/MultiQC/", mode: 'copy'

  output:
  file("multiqc_config.yaml") into ( multiQCconfig_fastq , multiQCconfig_sample )

  script:
  """
  touch multiqc_config.yaml
  echo "custom_logo_title: 'Exome Analysis Workflow'" >> multiqc_config.yaml
  echo "extra_fn_clean_exts:" >> multiqc_config.yaml
  echo "- _R1" >> multiqc_config.yaml
  echo "- _R2" >> multiqc_config.yaml
  echo "report_header_info:" >> multiqc_config.yaml
  echo "- Exoseq version: ${params.version}" >> multiqc_config.yaml
  echo "- Command Line: ${workflow.commandLine}" >> multiqc_config.yaml
  echo "- Directory: ${workflow.launchDir}" >> multiqc_config.yaml
  echo "- Genome: "${params.gfasta} >> multiqc_config.yaml
  echo "  dbSNP : ${params.dbsnp}" >> multiqc_config.yaml
  echo "  Omni: ${params.omni}" >> multiqc_config.yaml
  echo "  Mills: ${params.mills}" >> multiqc_config.yaml
  echo "top_modules:" >> multiqc_config.yaml
  echo "- 'fastqc'" >> multiqc_config.yaml
  echo "- 'cutadapt'" >> multiqc_config.yaml
  echo "- 'bwa'" >> multiqc_config.yaml
  echo "- 'samtools'" >> multiqc_config.yaml
  echo "- 'gatk'" >> multiqc_config.yaml
  """

}

/*
* Step 12 - Collect metrics, stats and other resources with MultiQC 
*/ 

process multiqc_fastq {

    publishDir "${params.outdir}/MultiQC/Library", mode: 'copy'

    input:
    file multiQCconfig from multiQCconfig_fastq
    file (fastqc:'fastqc/*') from trimgalore_fastqc_reports.flatten().toList()
    file ('trimgalore/*') from trimgalore_results.flatten().toList()
    file ('software_versions/*') from software_versions_yaml.toList()

    output:
    file '*multiqc_report_fastq.html' into multiqc_report_fastq

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report_fastq" : ''
    
    """
    multiqc -f $rtitle $rfilename --config $multiQCconfig . 
    """
}

process multiqc_sample {

    publishDir "${params.outdir}/MultiQC/Sample", mode: 'copy'

    input:
    file multiQCconfig from multiQCconfig_sample
    file('*')  from markdup_results.flatten().toList()
    file('*') from CollectMultipleMetricsOutput.flatten().toList()
    file('*') from HybridCaptureMetricsOutput.flatten().toList()
    file ('software_versions/*') from software_versions_yaml.toList()

    output:
    file '*multiqc_report_sample.html' into multiqc_report_sample

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report_sample" : ''
    
    """
    multiqc -f $rtitle $rfilename --config $multiQCconfig .
    """

}


/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def exoMessage() {
  // Display ExoSeq message
  log.info "nf-core/ExoSeq ANALYSIS WORKFLOW ~ ${params.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outdir
  log.info "Genome      : " + params.gfasta
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def versionMessage() {
  // Display version message
  log.info "nf-core/ExoSeq ANALYSIS WORKFLOW"
  log.info "  version   : " + version
  log.info workflow.commitId ? "Git info    : ${workflow.repository} - ${workflow.revision} [${workflow.commitId}]" : "  revision  : " + this.grabRevision()
}

workflow.onComplete {
  // Display complete message
  this.nextflowMessage()
  this.exoMessage()
  this.minimalInformationMessage()
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError {
  // Display error message
  this.nextflowMessage()
  this.exoMessage()
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}

def extractFastqFromDir(pattern) {

 // Curtesy of: Sarek developers
  // create a channel of FASTQs from a directory
  // All samples are considered 'normal'.
  // All FASTQ files are collected and emitted;
  // they must have _R1_ and _R2_ in their names.

  def fastq = Channel.create()

  // a temporary channel does all the work
  Channel.fromFilePairs(pattern)
    .ifEmpty { error "No files found  matching pattern '${pattern}'" }
    .subscribe onNext: { name,reads ->
      // the name of the library minus the lane number is assumed to be a unique sample id
      libraryId = name.toString().split('_L0')[0]
      (flowcell, lane) = flowcellLaneFromFastq(reads[0])
      patient = libraryId
      sampleId = libraryId
      libraryId = libraryId
      rgId = "${flowcell}.${sampleId}.${lane}"
      result = [patient, sampleId,libraryId, rgId, reads[0], reads[1] ]
      fastq.bind(result)
    }, onComplete: { fastq.close() }

  fastq
}

def flowcellLaneFromFastq(path) {
  // parse first line of a FASTQ file (optionally gzip-compressed)
  // and return the flowcell id and lane number.
  // expected format:
  // xx:yy:FLOWCELLID:LANE:... (seven fields)
  // or
  // FLOWCELLID:LANE:xx:... (five fields)
  InputStream fileStream = new FileInputStream(path.toFile())
  InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
  Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
  BufferedReader buffered = new BufferedReader(decoder)
  def line = buffered.readLine()
  assert line.startsWith('@')
  line = line.substring(1)
  def fields = line.split(' ')[0].split(':')
  String fcid
  int lane
  if (fields.size() == 7) {
    // CASAVA 1.8+ format
    fcid = fields[2]
    lane = fields[3].toInteger()
  }
  else if (fields.size() == 5) {
    fcid = fields[0]
    lane = fields[1].toInteger()
  }
  [fcid,lane]
}



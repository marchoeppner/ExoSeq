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

 Some code parts were lent from other NGI-Pipelines (e.g. CAW), specifically the error
 handling, logging messages. Thanks for that @CAW guys.
----------------------------------------------------------------------------------------
Developed based on GATK's best practise, takes set of FASTQ files and performs:
 - alignment (BWA)
 - recalibration (GATK)
 - realignment (GATK)
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
--reads                        Absolute path to project directory
--genome                       Name of iGenomes reference


Output:
--outdir                       Path where the results to be saved [Default: './results']

Kit files:
--kit                          Kit used to prep samples [Default: 'agilent_v5']
--bait                         Absolute path to bait file
--target                       Absolute path to target file
--target_bed                   Absolute path to target bed file (snpEff compatible format)

Genome/Variation files:
--dbsnp                        Absolute path to dbsnp file
--thousandg                    Absolute path to 1000G file
--mills                        Absolute path to Mills file
--omni                         Absolute path to Omni file
--gfasta                       Absolute path to genome fasta file
--bwa_index                    Absolute path to bwa genome index

Other options:
--exome                        Exome data, if this is not set, run as genome data
--project                      Uppnex project to user for SLURM executor

For more detailed information regarding the parameters and usage refer to package
documentation at https://github.com/nf-core/ExoSeq
"""

// Variables and defaults
params.name = false
params.help = false
params.reads = false
params.singleEnd = false
params.genome = false
params.run_id = false
params.exome = false //default genome, set to true to run restricting to exome positions
params.aligner = 'bwa' //Default, but stay tuned for later ;-) 
params.saveReference = true

// Output configuration
params.outdir = './results'
params.saveAlignedIntermediates = false
params.saveIntermediateVariants = false

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
params.target_bed = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].target_bed ?: false : false
params.dbsnp = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].dbsnp ?: false : false
params.thousandg = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].thousandg ?: false : false
params.mills = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].mills ?: false : false
params.omni = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].omni ?: false : false
params.gfasta = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].gfasta ?: false : false
params.bwa_index = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].bwa_index ?: false : false

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
if (!params.reads || !params.genome){
    exit 1, "Parameters '--reads' and '--genome' are required to run the pipeline"
}
if (!params.kitFiles[ params.kit ] && ['bait', 'target'].count{ params[it] } != 2){
    exit 1, "Kit '${params.kit}' is not available in pre-defined config, so " +
            "provide all kit specific files with option '--bait' and '--target'"
}
 if (!params.metaFiles[ params.genome ] && ['gfasta', 'bwa_index', 'dbsnp', 'thousandg', 'mills', 'omni'].count{ params[it] } != 6){
     exit 1, "Genome '${params.genome}' is not available in pre-defined config, so you need to provide all genome specific " +
             "files with options '--gfasta', '--bwa_index', '--dbsnp', '--thousandg', '--mills' and '--omni'"
 }

// Create a channel for input files

Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming }


// Validate Input indices for BWA Mem and GATK

if(params.aligner == 'bwa' ){
    bwaId = Channel
        .fromPath("${params.gfasta}.bwt")
        .ifEmpty { exit 1, "BWA index not found: ${params.gfasta}.bwt" }
}


// Create a summary for the logfile
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']       = params.genome
summary['WES/WGS']      = params.exome ? 'WES' : 'WGS'
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
if(params.aligner == 'bwa'){
    summary['Aligner'] = "BWA Mem"
    if(params.bwa_index)          summary['BWA Index']   = params.bwa_index
    else if(params.gfasta)          summary['Fasta Ref']    = params.gfasta
}
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

// Build BWA Index if this is required

if(params.aligner == 'bwa' && !params.bwa_index){
    // Create Channels
    fasta_for_bwa_index = Channel
        .fromPath("${params.gfasta}")
    fasta_for_samtools_index = Channel
        .fromPath("${params.gfasta}")
    // Create a BWA index for non-indexed genomes
    process makeBWAIndex {
        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_for_bwa_index

        output:
        file "*.{amb,ann,bwt,pac,sa}" into bwa_index

        script:
        """
        bwa index $fasta
        """
    }
    // Create a FastA index for non-indexed genomes
    process makeFastaIndex {
        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
        
        input:
        file fasta from fasta_for_samtools_index

        output:
        file "*.fai" into samtools_index

        script:
        """
        samtools faidx $fasta
        """
    }
} else {
    bwa_index = file("${params.bwa_index}")
}


/*
* STEP 0 - FastQC
*
*/

process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    file '.command.out' into fastqc_stdout

    script:
    """
    fastqc -q $reads
    """
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
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy'

        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file(reads) into trimmed_reads
        file '*trimming_report.txt' into trimgalore_results, trimgalore_logs

        script:
        single = reads instanceof Path
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (params.singleEnd) {
            """
            trim_galore --gzip $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}


/*
 * STEP 2 - Map with BWA Mem
 */

process bwamem {
    tag "$name"

    input:
    set val(name), file(reads) from trimmed_reads
    file(bwa_index) from bwa_index

    output:
    set val(name), file("${name}_bwa.${bam_suffix}") samples_sorted_bam
    file '.command.log' into bwa_stdout

    script:
    
    def avail_mem = task.memory ? "-m ${task.memory.toMega().intdiv(task.cpus)}M" : ''
    rg="\'@RG\\tID:${params.run_id}\\tSM:${params.run_id}\\tPL:illumina\'"
    
    """
    bwa mem \\
    -R $rg \\
    -t ${task.cpus} \\
    $params.gfasta \\
    $reads \\
    || samtools ${avail_mem} sort -O bam - > $outfile ${name}_bwa.bam
    """
}

/*
*  STEP 3 - Mark PCR duplicates in sorted BAM file
*/ 

process markDuplicates {
    tag "${name}"
    publishDir "${params.outdir}/Picard_Markduplicates/metrics", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".dup_metrics") > 0 ? filename : null }

    input:
    set val(name), file(sorted_bam) from samples_sorted_bam

    output:
    set val(name), file("${name}_markdup.bam"), file("${name}_markdup.bai") into samples_markdup_bam, samples_for_applyBQSR
    file("${name}.dup_metrics") into markdup_results
    file '.command.log' into markDuplicates_stdout

    script:
    """
        mkdir `pwd`/tmp
        gatk-launch MarkDuplicates \\
        --INPUT $sorted_bam \\
        --OUTPUT ${name}_markdup.bam \\
        --METRICS_FILE ${name}.dup_metrics \\
        --REMOVE_DUPLICATES false \\
        --CREATE_INDEX true \\
        --TMP_DIR=`pwd`/tmp \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}




/*
 * Step 5 - Recalibrate BAM file with known variants and BaseRecalibrator
 * 
*/ 
process recal_bam_files {
    tag "${name}"
    publishDir "${params.outdir}/GATK_Recalibration", mode: 'copy'


    input:
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_markdup_bam

    output:
    set val(name), file("${name}_table.recal") into samples_recal_reports
    file '.command.log' into gatk_stdout
    file '.command.log' into gatk_base_recalibration_results

    script:
    if(params.exome){
    """
    gatk-launch BaseRecalibrator \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        -O ${name}_table.recal \\
        -L $params.target \\
        --known-sites $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    } else {
    """
    gatk-launch BaseRecalibrator \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        -O ${name}_table.recal \\
        --known-sites $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """    
    }
}

process applyBQSR {
    tag "${name}"
    publishDir "${params.outdir}/GATK_ApplyBQSR", mode: 'copy'

    input:
    set val(name), file("${name}_table.recal") from samples_recal_reports
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_for_applyBQSR

    output:
    set val(name), file("${name}.bam"), file("${name}.bai") into bam_vcall, bam_metrics

    script:
    if(params.exome){
    """
    gatk-launch ApplyBQSR \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        --bqsr-recal-file ${name}_table.recal \\
        -O ${name}.bam \\
        -L $params.target \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    } else {
    """
    gatk-launch ApplyBQSR \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        --bqsr-recal-file ${name}_table.recal \\
        -O ${name}.bam \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """    
    }

}

/*
 * Step 7 - Determine quality metrics of mapped BAM files using QualiMap 2
 * 
*/ 
process qualiMap {
    tag "${name}"
    publishDir "${params.outdir}/Qualimap", mode: 'copy'
    
    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_metrics

    output:
    file "${name}" into qualimap_results
    file '.command.log' into qualimap_stdout

    script:
    gcref = ''
    gff = ''
    if(params.genome == 'GRCh37') gcref = '-gd HUMAN'
    if(params.genome == 'GRCm38') gcref = '-gd MOUSE'
    if(params.exome) gff ="-gff ${params.target_bed}"
    """
    qualimap bamqc $gcref \\
    -bam $realign_bam \\
    -outdir ${name} \\
    --skip-duplicated \\
    --collect-overlap-pairs \\
    -nt ${task.cpus} \\
    $gff \\
    --java-mem-size=${task.memory.toGiga()}G \\
    """
}

/*
 * Step 8 - Call Variants with HaplotypeCaller in GVCF mode (differentiate between exome and whole genome data here)
 * 
*/ 
process variantCall {
    tag "${name}"
    publishDir "${params.outdir}/GATK_VariantCalling/", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/variants/, "raw_variants")}

    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_vcall

    output:
    set val(name), file("${name}_variants.vcf"), file("${name}_variants.vcf.idx") into raw_variants

    script:
    if(params.exome){
    """
    gatk-launch HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -O ${name}_variants.vcf \\
        -ERC GVCF \\
        -L $params.target \\
        --create-output-variant-index \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    } else { //We have a winner (genome)
    """
    gatk-launch HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -O ${name}_variants.vcf \\
        -ERC GVCF \\
        --create-output-variant-index \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """    
    }
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
    qualimap --version &> v_qualimap.txt
    gatk-launch BaseRecalibrator --version &> v_gatk.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/**
* Step 10 - Generate MultiQC config file
*
*/ 

process GenerateMultiQCconfig {
  publishDir "${params.outdir}/MultiQC/", mode: 'copy'

  input:

  output:
  file("multiqc_config.yaml") into multiQCconfig

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
  echo "- 'qualimap'" >> multiqc_config.yaml
  echo "- 'gatk'" >> multiqc_config.yaml
  """
}


/*
* Step 12 - Collect metrics, stats and other resources with MultiQC in a single call
*/ 

process multiqc {
    tag "$name"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiQCconfig
    file (fastqc:'fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('gatk_base_recalibration/T*') from gatk_base_recalibration_results.toList()
    file ('gatk_picard_duplicates/*') from markdup_results.toList()
    file ('qualimap/*') from qualimap_results.toList()
    file ('software_versions/*') from software_versions_yaml.toList()


    output:
    file '*multiqc_report.html' into multiqc_report
    file '*_data' into multiqc_data
    file '.command.err' into multiqc_stderr
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
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


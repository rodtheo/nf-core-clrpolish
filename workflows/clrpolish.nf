/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowClrpolish.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.preview.recursion=true

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MERYL_COUNT as MERYL_COUNT_READS_01;
          MERYL_COUNT as MERYL_COUNT_GENOME_01  } from '../modules/nf-core/meryl/count/main'
include { CAT_FASTQ } from '../modules/nf-core/cat/fastq/main'
include { MERYL_HISTOGRAM as MERYL_HISTOGRAM_READS_PRE;
          MERYL_HISTOGRAM as MERYL_HISTOGRAM_GENOME_PRE } from '../modules/nf-core/meryl/histogram/main'

include { MERQURY as MERQURY_PRE } from '../modules/nf-core/merqury/main'
include { GENOMESCOPE2 as GENOMESCOPE2_PRE } from '../modules/nf-core/genomescope2/main'
include { MERFIN_COMPLETENESS } from '../modules/local/merfin_completeness'
include { MERFIN_HIST } from '../modules/local/merfin_hist'
include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'
include { BWA_INDEX } from '../modules/nf-core/bwa/index/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { FREEBAYES } from '../modules/nf-core/freebayes/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { TRIMMOMATIC } from '../modules/nf-core/trimmomatic/main'
include { MOSDEPTH } from '../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTWGSMETRICS } from '../modules/nf-core/picard/collectwgsmetrics/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTER;
          BCFTOOLS_VIEW as BCFTOOLS_VIEW_COMPRESS } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_BEFORE;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_POLISHED;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_CONCAT } from '../modules/nf-core/bcftools/index/main'
include { MERFIN_POLISH } from '../modules/local/merfin_polish'
include { FREEBAYES_FASTAGENERATEREGIONS } from '../modules/local/freebayes/fastagenerateregions'
include { BCFTOOLS_CONCAT } from '../modules/nf-core/bcftools/concat/main'
include { VCFLIB_VCFUNIQ } from '../modules/nf-core/vcflib/vcfuniq/main'
include { BCFTOOLS_SORT } from '../modules/nf-core/bcftools/sort/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

process GET_PEAK {
    input:
    tuple val(meta), path(genomescope2model)

    output:
    path("*_peak.txt"), emit: peak
// cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))' > file

    // peak="\$(cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))')"
    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat ${genomescope2model} | grep "^kmercov" | python -c 'import sys; print(float([x.split() for x in sys.stdin.readlines()][0][1]))' > ${prefix}_peak.txt
    """
}


process tick {
  input:
    path 'input.txt'
  output:
    path 'result.txt'
  script:
    """
    cat input.txt > result.txt
    echo "Task ${task.index} : tick" >> result.txt
    """
}

process tock {
  input:
    path 'input.txt'
  output:
    path 'result.txt'
  script:
    """
    cat input.txt > result.txt
    echo "Task ${task.index} : tock" >> result.txt
    """
}

workflow clock {
  take: infile
  main:
    infile | tick | tock
  emit:
    tock.out
}

workflow CLRPOLISH {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )

    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    reads_ch = INPUT_CHECK.out.reads
    genome_path_ch = Channel.fromPath("$params.genome")
    

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    TRIMMOMATIC (
        INPUT_CHECK.out.reads
    )

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Create a channel where paired-end data is mentioned as single-end
    // This is necessary in order to concatenate paired-end in a single file by cat/fastq nf-core module
    reads_ch_single = reads_ch.map { it -> [['id': it[0]['id'], 'single_end': true], it[1]] }
            

    CAT_FASTQ (
        reads_ch_single
    )

    MERYL_COUNT_READS_01 (
        CAT_FASTQ.out.reads
    )

    // Prepare genome channel to count kmers with meryl
    ch_genome_description = Channel.of(['id': 'genome', 'single_end': true])
    genome_ch = ch_genome_description.concat( genome_path_ch ).toList()
    genome_ch.view()

    MERYL_COUNT_GENOME_01 (
        genome_ch
    )

    MERYL_HISTOGRAM_READS_PRE (
        // MERYL_COUNT_READS_01.out.meryl_db.mix(MERYL_COUNT_GENOME_01.out.meryl_db)
        MERYL_COUNT_READS_01.out.meryl_db
    )

    MERYL_HISTOGRAM_GENOME_PRE (
        MERYL_COUNT_GENOME_01.out.meryl_db
    )

    MERYL_COUNT_READS_01.out.meryl_db.combine(genome_path_ch).first().view()

    MERQURY_PRE (
        MERYL_COUNT_READS_01.out.meryl_db.combine(genome_path_ch).first()
    )

    GENOMESCOPE2_PRE (
        MERYL_HISTOGRAM_READS_PRE.out.hist
    )

    peak_out = GET_PEAK (
        GENOMESCOPE2_PRE.out.model
    )

    peak_ch_val = peak_out.map{ it.text.trim() }.first()
    peak_ch_val.view()

    MERFIN_COMPLETENESS (
        MERYL_COUNT_GENOME_01.out.meryl_db,
        MERYL_COUNT_READS_01.out.meryl_db,
        GENOMESCOPE2_PRE.out.lookup_table,
        peak_ch_val
    )

    MERFIN_HIST (
        genome_ch,
        MERYL_COUNT_READS_01.out.meryl_db,
        GENOMESCOPE2_PRE.out.lookup_table,
        peak_ch_val
    )

    // 
    // MODULE: ILLUMINA READ MAPPING
    // 
    genome_index_ch = BWA_INDEX (
        genome_ch
    )

    
    BWA_MEM (
        TRIMMOMATIC.out.trimmed_reads,
        genome_index_ch.index,
        true
    )

    SAMTOOLS_INDEX (
        BWA_MEM.out.bam
    )


    genome_ch_fai = genome_ch.map{ it[0] }.concat( genome_ch.map{ it[1]+".fai" } ).toList()
    genome_ch_fai.view()
    
    SAMTOOLS_FAIDX (
        genome_ch,
        genome_ch_fai
    )

    ch_bam = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)       // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_bam.map { meta, bam, bai ->
                        return [meta, bam, bai, []]
            }
            .set { ch_mosdepth_in }

    MOSDEPTH ( ch_mosdepth_in, [[:],[]])
    PICARD_COLLECTWGSMETRICS ( BWA_MEM.out.bam
                                        .join(SAMTOOLS_INDEX.out.bai),
				genome_ch,
				genome_ch_fai,
				[]
			)

    
    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowClrpolish.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowClrpolish.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTWGSMETRICS.out.metrics.collect{it[1]}.ifEmpty([]))

    // println("A partir daqui")
    // ch_multiqc_files.view()

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

    // MODULE: VariantCalling
        
    genome_ch_fai_file = genome_ch_fai.map { it[1] }

    freebayes_chr = Channel.value("bsubtilis1 bsubtilis2")

    FREEBAYES_FASTAGENERATEREGIONS (
        genome_ch_fai,
        freebayes_chr
    )

    freebayes_input_ch = BWA_MEM.out.bam.flatten().concat(SAMTOOLS_INDEX.out.bai.map { it[1] }).toList()

    freebayes_input_ch_el = freebayes_input_ch.map( { [it] } )
    // freebayes_input_ch_el = Channel.of(1, 2)
    freebayes_input_ch_el.view{ "Freebayes: " + it }
    freebayes_beds = FREEBAYES_FASTAGENERATEREGIONS.out.bed.map { it[1] }
    freebayes_input_lists = freebayes_input_ch_el.combine(freebayes_beds.flatten())
    freebayes_input_lists_ok = freebayes_input_lists.map{ [ [ 'id': it[0][0]['id']+'_'+it[1].name.split('\\.')[1]+'_'+it[1].name.split('\\.')[3], 'single_end': false ],
                                 it[0][1], it[0][2], [], [], it[1]] }
    freebayes_input_lists_ok.view{ "LISTS: " + it }
    freebayes_input_lists_ok.count().view()                          

    FREEBAYES (
        freebayes_input_lists_ok,
        genome_path_ch.first(),
        genome_ch_fai_file.first(),
        [],
        [],
        []
    )

    BCFTOOLS_INDEX_BEFORE (
        FREEBAYES.out.vcf
    )

    out_vcfs = FREEBAYES.out.vcf.map{ it[1] }.collect()
    out_vcfs_tbis = BCFTOOLS_INDEX_BEFORE.out.tbi.map{ it[1] }.collect()

    // // out_vcfs.view {"AAAAH: " + it  }

    sample_name_ch = BWA_MEM.out.bam.map{ it[0] }
    // // out_vcfs.map{ it[1] }.toList().view{ "AAAAH: " + it }
    vcfs_ch_list = sample_name_ch.concat(out_vcfs).concat(out_vcfs_tbis).toList()
    // vcfs_ch_list.view { "VCFS: " + it}
    

    // BCFTOOLS_INDEX_BEFORE.out.tbi.view { "INDEX: " + it }

    // vcfs_ch_index_list = vcfs_ch_list.concat(BCFTOOLS_INDEX_BEFORE.out.tbi)
    // vcfs_ch_index_list.view { "VCFS_TBI: " + it}
    

    BCFTOOLS_CONCAT (
        vcfs_ch_list
    )

    BCFTOOLS_SORT (
        BCFTOOLS_CONCAT.out.vcf
    )

    BCFTOOLS_INDEX_CONCAT (
        BCFTOOLS_SORT.out.vcf
    )

    uniqvcf_ch = BCFTOOLS_SORT.out.vcf.concat(BCFTOOLS_INDEX_CONCAT.out.tbi)
    uniqvcf_ch.view{ "UNIQ: " + it}

    // VCFLIB_VCFUNIQ (
    //     BCFTOOLS_CONCAT.out.vcf.concat(BCFTOOLS_INDEX_CONCAT.out.tbi)
    // )

    // FREEBAYES_FASTAGENERATEREGIONS.out.bed.view()

    // BCFTOOLS_INDEX_BEFORE (
    //     FREEBAYES.out.vcf
    // )

    // BCFTOOLS_VIEW_FILTER (
    //     FREEBAYES.out.vcf,
    //     BCFTOOLS_INDEX.out.tbi,
    //     [],
    //     [],
    // )

    // HERE

    // MERFIN_POLISH (
    //     genome_ch,
    //     MERYL_COUNT_READS_01.out.meryl_db,
    //     GENOMESCOPE2_PRE.out.lookup_table,
    //     peak_ch_val,
    //     FREEBAYES.out.vcf 
    // )

    // BCFTOOLS_INDEX_POLISHED (
    //     MERFIN_POLISH.out.vcf
    // )


    // BCFTOOLS_VIEW_COMPRESS (
    //     MERFIN_POLISH.out.vcf.join(BCFTOOLS_INDEX_POLISHED.out.tbi),
    //     [],
    //     [],
    //     [],
    // )

    // TO HERE



    // FREEBAYES.out.view()
//   clock
//     .recurse(file(params.input))
//     .until { it -> it.size() > 100 }

//   clock
//     .out
//     .view(it -> it.text)


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

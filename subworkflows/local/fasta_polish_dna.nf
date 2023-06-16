// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'


include { MERYL_COUNT as MERYL_COUNT_READS_01;
          MERYL_COUNT as MERYL_COUNT_GENOME  } from '../../modules/nf-core/meryl/count/main'
include { CAT_FASTQ } from '../../modules/nf-core/cat/fastq/main'
include { MERYL_HISTOGRAM as MERYL_HISTOGRAM_READS_PRE;
          MERYL_HISTOGRAM as MERYL_HISTOGRAM_GENOME_PRE } from '../../modules/nf-core/meryl/histogram/main'

include { MERQURY as MERQURY_PRE } from '../../modules/nf-core/merqury/main'
include { GENOMESCOPE2 as GENOMESCOPE2_PRE } from '../../modules/nf-core/genomescope2/main'
include { MERFIN_COMPLETENESS } from '../../modules/local/merfin_completeness'
include { MERFIN_HIST as MERFIN_HIST_EVALUATE_POLISH;
          MERFIN_HIST as MERFIN_HIST_EVALUATE_PRIMARY_ASSEMBLY } from '../../modules/local/merfin_hist'
include { BWA_MEM } from '../../modules/nf-core/bwa/mem/main'
include { BWA_INDEX } from '../../modules/nf-core/bwa/index/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/samtools/index/main'
include { FREEBAYES } from '../../modules/nf-core/freebayes/main'
include { SAMTOOLS_FAIDX } from '../../modules/nf-core/samtools/faidx/main'
include { TRIMMOMATIC } from '../../modules/nf-core/trimmomatic/main'
include { MOSDEPTH } from '../../modules/nf-core/mosdepth/main'
include { PICARD_COLLECTWGSMETRICS } from '../../modules/nf-core/picard/collectwgsmetrics/main'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTER;
          BCFTOOLS_VIEW as BCFTOOLS_VIEW_COMPRESS } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_BEFORE;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_POLISHED;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_CONCAT } from '../../modules/nf-core/bcftools/index/main'
include { MERFIN_POLISH } from '../../modules/local/merfin_polish'
include { FREEBAYES_FASTAGENERATEREGIONS } from '../../modules/local/freebayes/fastagenerateregions'
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { VCFLIB_VCFUNIQ } from '../../modules/nf-core/vcflib/vcfuniq/main'
include { BCFTOOLS_SORT } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP as TABIX_BGZIP_VCF_POLISHED } from '../../modules/nf-core/tabix/bgzip/main'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/bcftools/consensus/main'

process SUBSET_VCF {
    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

  input:
    tuple val(meta), path(vcf_in)
  output:
    tuple val(meta), path("*_filtered.vcf.gz"), emit: vcf_out
  script:
    """
    #!/bin/bash

    zcat $vcf_in | head -n 10000 > result_filtered.vcf
    bgzip -c result_filtered.vcf > result_filtered.vcf.gz
    tabix -p vcf result_filtered.vcf.gz
    """
}

process GET_CHR_NAMES {
    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

  input:
    tuple path(fai)
  output:
    env str_chr_names, emit: chr_names
  script:
    """
    #!/bin/bash

    str_chr_names=`cut -f1 ${fai} | tr '\\n' ' ' | awk '{print \$0}'`
    """
}


workflow FASTA_POLISH_DNA {

    take:
    // TODO nf-core: edit input (take) channels
    ch_iteration
    ch_reads // channel: [ val(meta), [ TRIMMOMATIC.out.trimmed_reads ] ]
    ch_genome // channel: [ val(meta), [ assembly ] ]
    ch_read_meryl_db            // MERYL_COUNT_READS_01.out.meryl_db,
    ch_lookup_table            // GENOMESCOPE2_PRE.out.lookup_table,
    ch_peak_val            // peak_ch_val,
                

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    genome_index_ch = BWA_INDEX (
        ch_genome
    )

    
    BWA_MEM (
        ch_reads,
        genome_index_ch.index,
        true
    )

    SAMTOOLS_INDEX (
        BWA_MEM.out.bam
    )


    // genome_ch_fai = ch_genome.map{ it[0] }.concat( ch_genome.map{ it[1]+".fai" } ).toList()
    // genome_ch_fai.view()
    
    SAMTOOLS_FAIDX (
        ch_genome
        // genome_ch_fai
    )
    genome_ch_fai = SAMTOOLS_FAIDX.out.fai

    // MODULE: VariantCalling
        
    genome_ch_fai_file = genome_ch_fai.map { it[1] }

    GET_CHR_NAMES (genome_ch_fai_file)

    freebayes_chr = GET_CHR_NAMES.out.chr_names

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
    // freebayes_input_lists_ok.count().view()                          

    FREEBAYES (
        freebayes_input_lists_ok,
        ch_genome.map { it[1] }.first(),
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

    uniqvcf_ch = BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_INDEX_CONCAT.out.tbi)
    uniqvcf_ch.view{ "UNIQ: " + it}

    VCFLIB_VCFUNIQ (
        uniqvcf_ch
    )

    
    if (params.config_profile_name == 'Test profile') {

        SUBSET_VCF (
            VCFLIB_VCFUNIQ.out.vcf
        )

        vcf_to_polish_ch = SUBSET_VCF.out.vcf_out
        // vcf_to_polish_ch = VCFLIB_VCFUNIQ.out.vcf
    } else {
        vcf_to_polish_ch = VCFLIB_VCFUNIQ.out.vcf 
    }

    MERFIN_POLISH (
            ch_genome,
            ch_read_meryl_db,
            ch_lookup_table,
            ch_peak_val,
            vcf_to_polish_ch
        )

    TABIX_BGZIP_VCF_POLISHED (
        MERFIN_POLISH.out.vcf
    )

     BCFTOOLS_INDEX_POLISHED (
        TABIX_BGZIP_VCF_POLISHED.out.output
    )

    // TO DO: Check validity of using this task
    BCFTOOLS_VIEW_COMPRESS (
        MERFIN_POLISH.out.vcf.join(BCFTOOLS_INDEX_POLISHED.out.tbi),
        [],
        [],
        [],
    )

    ch_n = ch_iteration.map { ['id': 'round_' + it[0], 'single_end': true] }
    consensus_fasta_ch = TABIX_BGZIP_VCF_POLISHED.out.output.join(BCFTOOLS_INDEX_POLISHED.out.tbi).join(ch_genome)
    consensus_fasta_ch = consensus_fasta_ch.map {[ ['id': it[0]['id'], 'single_end': true ], it[1], it[2], it[3] ]}
    consensus_fasta_ch = ch_n.concat(consensus_fasta_ch).toList().map{[ ['id': it[0]['id']+'_'+it[1][0]['id'], 'single_end': true], it[1][1], it[1][2], it[1][3] ]}
    consensus_fasta_ch.view { "CONSENSUS: " + it }

    BCFTOOLS_CONSENSUS (
        consensus_fasta_ch
    )
    
    ch_genome_polished = BCFTOOLS_CONSENSUS.out.fasta


    MERFIN_HIST_EVALUATE_POLISH (
        ch_genome_polished,
        ch_read_meryl_db,
        ch_lookup_table,
        ch_peak_val
    )

    MERYL_COUNT_GENOME (
        ch_genome
    )

    MERFIN_COMPLETENESS (
        MERYL_COUNT_GENOME.out.meryl_db,
        ch_read_meryl_db,
        ch_lookup_table,
        ch_peak_val
    )

    ch_merfin_hist = MERFIN_COMPLETENESS.out.log

    emit:
    ch_iteration
    ch_reads // channel: [ val(meta), [ TRIMMOMATIC.out.trimmed_reads ] ]
    ch_genome_polished // channel: [ val(meta), [ assembly ] ]
    ch_read_meryl_db            // MERYL_COUNT_READS_01.out.meryl_db,
    ch_lookup_table            // GENOMESCOPE2_PRE.out.lookup_table,
    ch_peak_val            // peak_ch_val,
    ch_merfin_hist
    // TODO nf-core: edit emitted channels
    
    // versions = ch_versions                     // channel: [ versions.yml ]
}


// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules
nextflow.enable.dsl = 2
nextflow.preview.recursion=true

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
          BCFTOOLS_VIEW as BCFTOOLS_VIEW_COMPRESS;
          BCFTOOLS_VIEW as BCFTOOLS_VIEW_COMPRESS_NORM } from '../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_BEFORE;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_POLISHED;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_CONCAT;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_NORM;
          BCFTOOLS_INDEX as BCFTOOLS_INDEX_COMPRESS_NORM } from '../../modules/nf-core/bcftools/index/main'
include { MERFIN_POLISH } from '../../modules/local/merfin_polish'
include { FREEBAYES_FASTAGENERATEREGIONS } from '../../modules/local/freebayes/fastagenerateregions'
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { VCFLIB_VCFUNIQ } from '../../modules/nf-core/vcflib/vcfuniq/main'
include { BCFTOOLS_SORT } from '../../modules/nf-core/bcftools/sort/main'
include { TABIX_BGZIP as TABIX_BGZIP_VCF_POLISHED;
            TABIX_BGZIP as TABIX_BGZIP_VCF_NORM } from '../../modules/nf-core/tabix/bgzip/main'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/bcftools/consensus/main'
include { GET_COMPLETNESS_VAL } from '../../modules/local/merfin_get_qv_value'
include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main'

process SUBSET_VCF {
    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

  input:
    tuple val(meta), path(vcf_in)
  output:
    tuple val(meta), path("*_filtered.vcf.gz"), emit: vcf_out

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    zcat $vcf_in | head -n 10000 > result_filtered.vcf
    bgzip -c result_filtered.vcf > result_filtered.vcf.gz
    tabix -p vcf result_filtered.vcf.gz
    """
}

process GET_QV_VALUE {
    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

  input:
    path(merfin_hist_log)
  output:
    env(qv_star_score)

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    qv_star_score=`grep "^Merfin QV*" ${merfin_hist_log} | awk '{split(\$0,a," "); print a[3]}'`
    
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

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    str_chr_names=`cut -f1 ${fai} | tr '\\n' ' ' | awk '{print \$0}'`
    """
}

params.n_chromosome = 2
params.n_splits_per_chr = 2

workflow FASTA_POLISH_DNA {

    take:
    // TODO nf-core: edit input (take) channels
    ch_tuple_interation
    // ch_iteration
    // ch_reads // channel: [ val(meta), [ TRIMMOMATIC.out.trimmed_reads ] ]
    // ch_genome // channel: [ val(meta), [ assembly ] ]
    // ch_read_meryl_db            // MERYL_COUNT_READS_01.out.meryl_db,
    // ch_lookup_table            // GENOMESCOPE2_PRE.out.lookup_table,
    // ch_peak_val            // peak_ch_val,
                

    main:

    ch_iteration = Channel.empty()
    ch_genome = Channel.empty()
    ch_qv_value = Channel.empty()
    ch_completness = Channel.empty()

    ch_iteration = ch_tuple_interation.map{ it[0] }
    ch_reads = ch_tuple_interation.map{ it[1] }
    ch_genome = ch_tuple_interation.map{ it[2] }
    ch_read_meryl_db = ch_tuple_interation.map{ it[3] }
    ch_lookup_table = ch_tuple_interation.map{ it[4] }
    ch_peak_val = ch_tuple_interation.map{ it[5] }
    ch_qv_value = ch_tuple_interation.map{ it[6] }
    ch_completness = ch_tuple_interation.map{ it[8] }

    // ch_tuple_interation.view{ "ITERATION CH: " + it}

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow
    // ch_genome.view{ "CH_GENOME: "+it }
    genome_index_ch = BWA_INDEX (
        ch_genome
    )

    ch_genome.view{ "BWA MEM EXECUTADO AQUI: "+it }
    
    BWA_MEM (
        ch_reads,
        genome_index_ch.index,
        true
    )

    genome_index_ch.index.view{ "BWA MEM OUTPUT AQUI: "+it }

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

    // FREEBAYES_FASTAGENERATEREGIONS.out.bed.view{ "FREEBAYES_FASTAGENERATEREGIONS: " + it }

    // BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai).view()
    // FREEBAYES_FASTAGENERATEREGIONS.out.bed.flatten().view()
    // SAMTOOLS_INDEX.out.bai.map { it[1] }.view()

    // freebayes_input_ch = BWA_MEM.out.bam.flatten().concat(SAMTOOLS_INDEX.out.bai.map { it[1] }).toList()
    freebayes_input_ch = BWA_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)
    // freebayes_input_ch.view{ "freebayes_input_ch: " + it }
    ch_genome_meta = ch_genome.map { it[0] }
    // ch_genome_meta.view{ "PQP: "+it }
    
    
    // freebayes_input_ch_el = freebayes_input_ch.map( { [it] } )
    // freebayes_input_ch_el = Channel.of(1, 2)
    // freebayes_input_ch_el.view{ "freebayes_input_ch_el: " + it }
    freebayes_beds = FREEBAYES_FASTAGENERATEREGIONS.out.bed.map { it[1] }
    // freebayes_input_ch_el = ch_genome.first().mix(freebayes_input_ch.map { [it[1], it[2]]}.first()).buffer(size: 2)
    // freebayes_input_ch_el = ch_genome.join(genome_ch_fai).mix(freebayes_input_ch.map { [it[1], it[2]]}.first()).buffer(size: 2)
    freebayes_input_ch_el = ch_genome.join(genome_ch_fai).merge(freebayes_input_ch.map { [it[1], it[2]]})
    
    freebayes_input_ch_el_fb = freebayes_input_ch_el.map{ [it[0][0], it[1][0], it[1][1]] }
    freebayes_input_lists = freebayes_input_ch_el.combine(freebayes_beds.flatten())
    // freebayes_input_lists.join(ch_genome).join(genome_ch_fai).view{ "******* INPUT CH: "+it}
    // FREEBAYES_FASTAGENERATEREGIONS.out.bed.view{ "PQP: "+it }
    freebayes_input_lists.filter{it -> { def name=it[0]['id']; it[5].name =~ /^$name/;} }.view{"XXXXXXXXXXXXXXXXXXX VIEW HERE:"+it}
    // freebayes_input_lists_ok = freebayes_input_lists.map{ [ [ 'id': it[0][0]['id']+'_'+it[1].name.split('\\.')[1]+'_'+it[1].name.split('\\.')[3], 'single_end': false ],
    //                              it[0][1], it[0][2], [], [], it[1]] }
    // freebayes_input_lists_ok = freebayes_input_lists.map{it[3]}.view {it.name.split('\\.')}
    freebayes_input_lists = freebayes_input_lists.filter{it -> { def name=it[0]['id']; it[5].name =~ /^$name/;} }
    freebayes_input_lists_ok = freebayes_input_lists.map{ 
        [ [ 'id': it[0]['id']+'_'+it[5].name.split('\\.')[1]+'_'+it[5].name.split('\\.')[3], 'single_end': false ],
                                 it[3], it[4], [], [], it[5], it[1], it[2]] }
    freebayes_input_lists_ok.first().view{ "LISTS: " + it }
    // freebayes_input_lists_ok.count().view()
    ch_genome.view{ "=== ESSE EH O GENOMA: "+it }  
    genome_ch_fai_file.collate(3, false).view{ "=== ESSE EH O GENOMA (FAI): "+it }


    FREEBAYES (
        freebayes_input_lists_ok.map{ [it[0], it[1], it[2], it[3], it[4], it[5]] },
        freebayes_input_lists_ok.map{ it[6] },
        freebayes_input_lists_ok.map{ it[7] },
        // ch_genome.map { it[1] }.collate(1).first(),
        // genome_ch_fai_file.collate(1).first(),
        [],
        [],
        []
    )

    out_vcfs = FREEBAYES.out.vcf.map{ [it[0], it[1]] }
    // out_vcfs.view{ "OOOUT: "+it }
    
    BCFTOOLS_INDEX_BEFORE (
        FREEBAYES.out.vcf
        // out_vcfs
    )

    out_vcfs = FREEBAYES.out.vcf.map{ it[1] }
    out_vcfs_tbis = BCFTOOLS_INDEX_BEFORE.out.tbi.map{ it[1] }

    // sample_name_ch = BWA_MEM.out.bam.map{ it[0] }
    // sample_name_ch = ch_genome.map{ it[0] }
    // sample_name_ch.view{ "sample name: "+it }
    // // out_vcfs.map{ it[1] }.toList().view{ "AAAAH: " + it }
    // ch_vcfs_list = out_vcfs.collect()
    // def size_to_buffer = freebayes_chr.map{ it.text.split(' ') }
    

    // ch_buffer_size = FREEBAYES_FASTAGENERATEREGIONS.out.bed.map{ it[1].size }
    
    buffer_size = params.n_chromosome * params.n_splits_per_chr
    println "BUFFER SIZE: $buffer_size"
    ch_vcfs_list = out_vcfs.collate( buffer_size )
    ch_vcfs_list = ch_genome_meta.merge(ch_vcfs_list.map{ [it] })
    

    // ch_vcfs_list.view{ "VCFS LIST:"+it }

    ch_vcfs_tbis_list = out_vcfs_tbis.collate( buffer_size )
    ch_vcfs_tbis_list = ch_genome_meta.merge(ch_vcfs_tbis_list.map{ [it] })
    // ch_vcfs_tbis_list.view{ "TBIS LIST:"+it }

    
    // ch_genome_meta.view{ "PQPQQQQQQ: " + it}
    vcfs_ch_list = ch_genome_meta.join(ch_vcfs_list).join(ch_vcfs_tbis_list)
    // vcfs_ch_list.view{ "PQPQQQQQQ: " + it}
    // vcfs_ch_list = vcfs_ch_list.merge(ch_vcfs_tbis_list).buffer(size: 2).map{ [it[0][0], it[0][1], it[1]] }
    // vcfs_ch_list.view{ "PQP: "+it }
    // vcfs_ch_list = vcfs_ch_list.mix(ch_vcfs_tbis_list).buffer(size: 3)
    
    // vcfs_ch_list.view{ "VLISSSSST: " + it }
    // vcfs_ch_list = vcfs_ch_list.mix(ch_vcfs_tbis_list.first())
    // vcfs_ch_list = vcfs_ch_list.collect()
    // vcfs_ch_list = vcfs_ch_list.mix(ch_vcfs_list)
    // vcfs_ch_list = vcfs_ch_list.collect()
    // vcfs_ch_list = vcfs_ch_list.mix(ch_vcfs_list)
    // vcfs_ch_list = vcfs_ch_list.mix(ch_vcfs_tbis_list)
    // // vcfs_ch_list = vcfs_ch_list.collect()
    // // vcfs_ch_list = sample_name_ch.concat(ch_vcfs_list)
    // // .mix(out_vcfs_tbis.toList())
    // ch_vcfs_list = Channel.empty()
    // ch_vcfs_tbis_list = Channel.empty()
    // vcfs_ch_list.view { "VCFS: " + it}
    

    // BCFTOOLS_INDEX_BEFORE.out.tbi.view { "INDEX: " + it }

    // vcfs_ch_index_list = vcfs_ch_list.concat(BCFTOOLS_INDEX_BEFORE.out.tbi)
    vcfs_ch_list.view { "VCFS_TBI: " + it}
    

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
    // uniqvcf_ch.view{ "UNIQ: " + it}

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

    ch_vcf_norm = TABIX_BGZIP_VCF_POLISHED.out.output.join(BCFTOOLS_INDEX_POLISHED.out.tbi)
    
    ch_vcf_norm.view{"VCF NORM: "+it}

    BCFTOOLS_NORM (
        ch_vcf_norm,
        ch_genome
    )

    BCFTOOLS_NORM.out.vcf.view{"UUUUUUULA LA: "+it}

     BCFTOOLS_INDEX_NORM (
        BCFTOOLS_NORM.out.vcf
    )

    // TO DO: Check validity of using this task
    BCFTOOLS_VIEW_COMPRESS_NORM (
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_INDEX_NORM.out.tbi),
        [],
        [],
        [],
    )

    BCFTOOLS_INDEX_COMPRESS_NORM (
        BCFTOOLS_VIEW_COMPRESS_NORM.out.vcf
    )

    ch_n = ch_iteration.map { ['id': 'round_' + it[0], 'single_end': true] }
    
    consensus_fasta_ch = BCFTOOLS_VIEW_COMPRESS_NORM.out.vcf.join(BCFTOOLS_INDEX_COMPRESS_NORM.out.tbi).join(ch_genome)
    
    consensus_fasta_ch = consensus_fasta_ch.map {[ ['id': it[0]['id'], 'single_end': true ], it[1], it[2], it[3] ]}
    
    // consensus_fasta_ch = ch_n.combine(consensus_fasta_ch)
    consensus_fasta_ch = ch_n.merge(consensus_fasta_ch)
    consensus_fasta_ch.view { "CONSENSUS B: " + it }
    consensus_fasta_ch.filter{ it[0] }.view { "CONSENSUS B: " + it }
    consensus_fasta_ch = consensus_fasta_ch.map{[ ['id': it[0]['id']+'_'+it[1]['id'], 'single_end': true], it[2], it[3], it[4] ]}
    // consensus_fasta_ch.view { "CONSENSUS C: " + it }

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

    // GET_QV_VALUE ( MERFIN_HIST_EVALUATE_POLISH.out.log ).view{"GET QV VALUE: " + it}

    MERYL_COUNT_GENOME (
        ch_genome_polished
    )

    MERFIN_COMPLETENESS (
        MERYL_COUNT_GENOME.out.meryl_db,
        ch_read_meryl_db,
        ch_lookup_table,
        ch_peak_val
    )

    ch_merfin_hist = MERFIN_HIST_EVALUATE_POLISH.out.log
    ch_merfin_hist.view{ "MERFIN HIST: "+it }
    ch_qv_value_polished = GET_QV_VALUE (ch_merfin_hist).view{ "QV VALUE: " + it }
    ch_completness_polished = GET_COMPLETNESS_VAL(MERFIN_COMPLETENESS.out.completeness)

    // ch_out_iteration = ch_iteration.map{ [[it[0] + 1]] }.merge(ch_reads.map{ [it] }).merge(ch_genome_polished.map{ [it] }).merge(ch_read_meryl_db.map{ [it] }).merge(ch_lookup_table.map{ [it] }).merge(ch_peak_val.first()).map{ [it] }
    ch_diff_qv_value = ch_qv_value_polished.toFloat().merge(ch_qv_value.toFloat()).map { it[0] - it[1] }.toFloat()
    
    
    
    // ch_out_iteration_ = ch_iteration.map{ [[it[0] + 1]] }.first().concat(ch_reads.map{ [it] }.first()).concat(ch_genome_polished.map{ [it] }.first()).concat(ch_read_meryl_db.map{ [it] }.first()).concat(ch_lookup_table.map{ [it] }.first()).concat(ch_peak_val.first()).buffer(size: 6)
    ch_out_iteration_ = ch_iteration.map{ [[it[0] + 1]] }.merge(ch_reads.map{ [it] }, ch_genome_polished.map{ [it] }, ch_read_meryl_db.map{ [it] }, ch_lookup_table.map{ [it] }, ch_peak_val, ch_qv_value_polished, ch_diff_qv_value, ch_completness_polished)
    // ch_out_iteration_ = ch_out_iteration_.map{ [ it[0][0], it[1][0], it[2][0], it[3][0], it[4][0], it[5] ] }

    ch_out_iteration_.view { 'OUT INTERACTION: ' + it }

    emit:
    ch_out_iteration_
    // ch_iteration.map{ [it[0] + 1] }.first()
    // ch_reads.first() // channel: [ val(meta), [ TRIMMOMATIC.out.trimmed_reads ] ]
    // ch_genome_polished.first() // channel: [ val(meta), [ assembly ] ]
    // ch_read_meryl_db.first()            // MERYL_COUNT_READS_01.out.meryl_db,
    // ch_lookup_table.first()            // GENOMESCOPE2_PRE.out.lookup_table,
    // ch_peak_val.first()            // peak_ch_val,
    // ch_merfin_hist
    // TODO nf-core: edit emitted channels
    
    // versions = ch_versions                     // channel: [ versions.yml ]
}


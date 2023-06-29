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
    missing_kmers_env=`grep "^K-mers not found in reads (missing)" ${merfin_hist_log} | awk '{split(\$0,a,": "); print a[2]}'`
    
    """
}

process GET_KMERS_VALS {
    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

  input:
    path(merfin_hist_log)
  output:
    tuple env(missing_kmers_env), env(overly_kmers_env), env(found_kmers_env)

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    missing_kmers_env=`grep "^K-mers not found in reads (missing)" ${merfin_hist_log} | awk '{split(\$0,a,": "); print a[2]}'`
    overly_kmers_env=`grep "^K-mers overly represented in assembly:" ${merfin_hist_log} | awk '{split(\$0,a,": "); print a[2]}'`
    found_kmers_env=`grep "^K-mers found in the assembly:" ${merfin_hist_log} | awk '{split(\$0,a,": "); print a[2]}'`
    """
}

process GET_COMPLETNESS_VAL {
    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

  input:
    tuple val(meta), path(merfin_completeness_log)
  output:
    env(complete_env)

  when:
    task.ext.when == null || task.ext.when

  script:
    """
    #!/bin/bash

    complete_env=`grep -E "^COMPLETENESS:" ${merfin_completeness_log} | awk '{split(\$0,a,": "); print a[2]}' | sed 's/^[ \t]*//;s/[ \t]*\$//'`
    """
}
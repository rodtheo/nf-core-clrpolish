process GET_FRAMESHIFTS {
    tag "$meta.id"
    label 'process_low'

    publishDir "$params.outdir/frameshifts", mode: 'copy'
    conda "bioconda::gawk=5.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/gawk:5.3.1' }"

    input:
        tuple val(meta), path(gff)

    output:
        tuple val(meta), path("*_frameshifts.txt"), emit: frameshifts
        path "versions.yml", emit: versions

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        grep -v "#" ${gff} | grep 'Frameshift' | gawk -F'\t' '\$3=="mRNA" && match(\$9, /Frameshift=([^;]+)/, arr) {sum += arr[1]} END {print "Total:", sum}' > ${prefix}_frameshifts.txt

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            grep: \$(grep --version | head -n 1)
            awk: \$(awk --version | head -n 1)
        END_VERSIONS
        """
    stub:
        """
        echo "Total: 0" > frameshifts.txt
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            grep: stub
            awk: stub
        END_VERSIONS
        """
}
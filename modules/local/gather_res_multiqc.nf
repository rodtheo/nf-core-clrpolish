process GATHER_RES_MULTIQC {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    val(list_infos)
	tuple val(meta_a), path(genomescope_linear_plot_png)
	tuple val(meta_c), path(merqury_spectra_cn_fl_png)
	tuple val(meta_d), path(merqury_spectra_asm_fl_png)

    output:
    path '*_mqc.json'       , emit: mqc
	path 'table_results_mqc.txt'       , emit: general_stats
	path 'genomescope2_plot_mqc.png'       , emit: gs2_plot
	path 'merqury_spectra_cn_fl_plot_mqc.png', emit: merqry_spectra_cn
	path 'merqury_spectra_asm_fl_plot_mqc.png', emit: merqry_spectra_asm
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/clrpolish/bin/
    """
    create_merfin_multiqc_out.py \\
        \'$list_infos\' \\
        qvalues_gathered_linegraph_mqc.json

	ln -s -f ${genomescope_linear_plot_png} genomescope2_plot_mqc.png
	ln -s -f ${merqury_spectra_cn_fl_png} merqury_spectra_cn_fl_plot_mqc.png
	ln -s -f ${merqury_spectra_asm_fl_png} merqury_spectra_asm_fl_plot_mqc.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

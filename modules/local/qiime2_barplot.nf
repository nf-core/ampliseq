// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_BARPLOT {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

	input:
	path(metadata)
	path(table)
	path(taxonomy)

	output:
	path("barplot/*")   , emit: folder
    path "*.version.txt", emit: version

    script:
    def software     = getSoftwareName(task.process)
	"""
    export XDG_CONFIG_HOME="\${PWD}/HOME"

	qiime taxa barplot  \
		--i-table ${table}  \
		--i-taxonomy ${taxonomy}  \
		--m-metadata-file ${metadata}  \
		--o-visualization taxa-bar-plots.qzv  \
		--verbose
	qiime tools export --input-path taxa-bar-plots.qzv  \
		--output-path barplot

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
	"""
}
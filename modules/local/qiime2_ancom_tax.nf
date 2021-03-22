// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_ANCOM_TAX {
    tag "${table.baseName} - taxonomic level: ${taxlevel}"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

	input:
    tuple path(metadata), path(table), path(taxonomy) ,val(taxlevel)

	output:
	path("ancom/*")     , emit: ancom
    path "*.version.txt", emit: version

    script:
    def software     = getSoftwareName(task.process)
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

	qiime taxa collapse \
		--i-table ${table} \
		--i-taxonomy ${taxonomy} \
		--p-level ${taxlevel} \
		--o-collapsed-table lvl${taxlevel}-${table}
	qiime composition add-pseudocount \
		--i-table lvl${taxlevel}-${table} \
		--o-composition-table comp-lvl${taxlevel}-${table}
	qiime composition ancom \
		--i-table comp-lvl${taxlevel}-${table} \
		--m-metadata-file ${metadata} \
		--m-metadata-column ${table.baseName} \
		--o-visualization comp-lvl${taxlevel}-${table.baseName}.qzv
	qiime tools export --input-path comp-lvl${taxlevel}-${table.baseName}.qzv \
		--output-path ancom/Category-${table.baseName}-level-${taxlevel}

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_EXPORT_RELTAX {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

	input:
	path(table)
	path(taxonomy)

	output:
	path("*.tsv")        , emit: tsv
    path "*.version.txt" , emit: version

    script:
    def software     = getSoftwareName(task.process)
	"""
    export XDG_CONFIG_HOME="\${PWD}/HOME"

	##on several taxa level
	array=( 2 3 4 5 6 7 )

	for i in \${array[@]}
	do
		#collapse taxa
		qiime taxa collapse \
			--i-table ${table} \
			--i-taxonomy ${taxonomy} \
			--p-level \$i \
			--o-collapsed-table table-\$i.qza
		#convert to relative abundances
		qiime feature-table relative-frequency \
			--i-table table-\$i.qza \
			--o-relative-frequency-table relative-table-\$i.qza
		#export to biom
		qiime tools export --input-path relative-table-\$i.qza \
			--output-path relative-table-\$i
		#convert to tab separated text file
		biom convert \
			-i relative-table-\$i/feature-table.biom \
			-o rel-table-\$i.tsv --to-tsv
	done

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
	"""
}
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process QIIME2_INASV {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(asv)
    
    output:
    path("table.qza")    , emit: qza
    path "*.version.txt" , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    echo -n "#OTU Table" | cat - "$asv" > biom-table.txt
    biom convert -i biom-table.txt -o table.biom --table-type="OTU table" --to-hdf5
    qiime tools import \
        --input-path table.biom \
        --type 'FeatureTable[Frequency]' \
        --input-format BIOMV210Format \
        --output-path table.qza
    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_INSEQ {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(seq)
    
    output:
    path("rep-seqs.qza"), emit: qza
    path "*.version.txt", emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    qiime tools import \
        --input-path "$seq" \
        --type 'FeatureData[Sequence]' \
        --output-path rep-seqs.qza
    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_INTAX {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "quay.io/qiime2/core:2021.2"
    } else {
        container "quay.io/qiime2/core:2021.2"
    }

    input:
    path(tax)
    
    output:
    path "*.version.txt" , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    #requires: Feature_ID\tTaxon, see https://forum.qiime2.org/t/dada2-taxonomy-into-qiime2/15369/2
    #qiime tools import \
    #    --type 'FeatureData[Taxonomy]' \
    #    --input-path $tax \
    #    --output-path taxonomy.qza
    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_EXTRACT {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_low'
    label 'single_cpu'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    tuple val(meta), path(database)
    
    output:
    tuple val(meta), path("*.qza"), emit: qza
    path "*.version.txt"          , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    ### Import
    qiime tools import --type \'FeatureData[Sequence]\' \
        --input-path ${database[0]} \
        --output-path ref-seq.qza
    qiime tools import --type \'FeatureData[Taxonomy]\' \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path ${database[1]} \
        --output-path ref-taxonomy.qza
    #Extract sequences based on primers
    qiime feature-classifier extract-reads \
        --i-sequences ref-seq.qza \
        --p-f-primer ${meta.FW_primer} \
        --p-r-primer ${meta.RV_primer} \
        --o-reads ${meta.FW_primer}-${meta.RV_primer}-ref-seq.qza \
        --quiet

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_TRAIN {
    tag "${meta.FW_primer}-${meta.RV_primer}"
    label 'process_high'
    label 'single_cpu'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    tuple val(meta), path(qza)
    
    output:
    path("*-classifier.qza"), emit: qza
    path "*.version.txt"    , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    #Train classifier
    qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads ${meta.FW_primer}-${meta.RV_primer}-ref-seq.qza \
        --i-reference-taxonomy ref-taxonomy.qza \
        --o-classifier ${meta.FW_primer}-${meta.RV_primer}-classifier.qza \
        --quiet

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_CLASSIFY {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(trained_classifier)
    path(repseq)
    
    output:
    path("taxonomy.qza"), emit: qza
    path("taxonomy.tsv"), emit: tsv
    path "*.version.txt", emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    qiime feature-classifier classify-sklearn  \
        --i-classifier ${trained_classifier}  \
        --p-n-jobs ${task.cpus}  \
        --i-reads ${repseq}  \
        --o-classification taxonomy.qza  \
        --verbose
    qiime metadata tabulate  \
        --m-input-file taxonomy.qza  \
        --o-visualization taxonomy.qzv  \
        --verbose
    #produce "taxonomy/taxonomy.tsv"
    qiime tools export --input-path taxonomy.qza  \
        --output-path taxonomy
    qiime tools export --input-path taxonomy.qzv  \
        --output-path taxonomy
    cp taxonomy/taxonomy.tsv .

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

process QIIME2_FILTERTAXA {
    tag "taxa:${exclude_taxa};min-freq:${min_frequency};min-samples:${min_samples}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(table)
    path(repseq)
    path(taxonomy)
    val(min_frequency)
    val(min_samples)
    val(exclude_taxa)

    output:
    path("filtered-table.qza"), emit: asv
    path("filtered-table.tsv"), emit: tsv
    path("filtered-sequences.qza"), emit: seq
    path "*.version.txt"       , emit: version

    script:
    def minfrequency = "${min_frequency}" == "false" ? 1 : "${min_frequency}"
    def minsamples   = "${params.min_samples}" == "false" ? 1 : "${params.min_samples}"
    def software     = getSoftwareName(task.process)
    """
    if ! [ \"${exclude_taxa}\" = \"none\" ]; then
        #filter sequences
        qiime taxa filter-seqs \
            --i-sequences ${repseq} \
            --i-taxonomy ${taxonomy} \
            --p-exclude ${exclude_taxa} --p-mode contains \
            --o-filtered-sequences tax_filtered-sequences.qza
        #filter abundance table
        qiime taxa filter-table \
            --i-table ${table} \
            --i-taxonomy ${taxonomy} \
            --p-exclude ${exclude_taxa} --p-mode contains \
            --o-filtered-table tax_filtered-table.qza
        filtered_table="tax_filtered-table.qza"
        filtered_sequences="tax_filtered-sequences.qza"
    else
        filtered_table=${table}
        filtered_sequences=${repseq}
    fi
    qiime feature-table filter-features \
        --i-table \$filtered_table \
        --p-min-frequency ${minfrequency} \
        --p-min-samples ${minsamples} \
        --o-filtered-table filtered-table.qza
    
    qiime feature-table filter-seqs \
        --i-data \$filtered_sequences \
        --i-table filtered-table.qza \
        --o-filtered-data filtered-sequences.qza

    #produce raw count table in biom format "table/feature-table.biom"
    qiime tools export --input-path filtered-table.qza  \
        --output-path table
    #produce raw count table
    biom convert -i table/feature-table.biom \
        -o table/feature-table.tsv  \
        --to-tsv
    cp table/feature-table.tsv filtered-table.tsv

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}

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

process QIIME2_EXPORT_ABSOLUTE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

	input:
	path(table)
	path(repseq)
	path(taxonomy)

	output:
	path("rep-seq.fasta")            , emit: fasta
	path("feature-table.tsv")        , emit: tsv
	path("feature-table.biom")       , emit: biom
	path("seven_number_summary.tsv") , emit: summary
    path("descriptive_stats.tsv")    , emit: descr
	path("abs-abund-table-*.tsv")    , emit: abundtable
    path "*.version.txt"             , emit: version

    script:
    def software     = getSoftwareName(task.process)
	"""
	#produce raw count table in biom format "table/feature-table.biom"
	qiime tools export --input-path ${table}  \
		--output-path table
    cp table/feature-table.biom .
    
	#produce raw count table "table/feature-table.tsv"
	biom convert -i table/feature-table.biom \
		-o feature-table.tsv  \
		--to-tsv
    
	#produce representative sequence fasta file "sequences.fasta"
	qiime feature-table tabulate-seqs  \
		--i-data ${repseq}  \
		--o-visualization rep-seqs.qzv
	qiime tools export --input-path rep-seqs.qzv  \
		--output-path representative_sequences
    cp representative_sequences/sequences.fasta rep-seq.fasta
    cp representative_sequences/*.tsv .

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
		#export to biom
		qiime tools export --input-path table-\$i.qza \
			--output-path table-\$i
		#convert to tab separated text file
		biom convert \
			-i table-\$i/feature-table.biom \
			-o abs-abund-table-\$i.tsv --to-tsv
	done

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
	"""
}

process QIIME2_EXPORT_RELASV {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

	input:
	path(table)

	output:
	path("rel-table-ASV.tsv"), emit: tsv
    path "*.version.txt"     , emit: version

    script:
    def software     = getSoftwareName(task.process)
	"""
	#convert to relative abundances
	qiime feature-table relative-frequency \
		--i-table ${table} \
		--o-relative-frequency-table relative-table-ASV.qza
    
	#export to biom
	qiime tools export --input-path relative-table-ASV.qza --output-path relative-table-ASV

	#convert to tab separated text file "rel-table-ASV.tsv"
	biom convert -i relative-table-ASV/feature-table.biom \
		-o rel-table-ASV.tsv --to-tsv

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
	"""
}

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
	##on several taxa level
	array=( 2 3 4 5 6 )

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
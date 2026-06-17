process QIIME2_ANCOMBC2_ASV {
    tag "${table.baseName}-${formula_in}"
    label 'process_medium'
    label 'process_long'
    label 'error_ignore'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    tuple path(metadata), path(table), val(formula_in)

    output:
    path("${outfolder}/*.svg")          , emit: plot
    path("${outfolder}/*.tsv")          , emit: table
    path("${outfolder}/*.qzv")          , emit: barplot_qzv
    path("${outfolder}/differentials/*"), emit: differentials
    path("*.qza")                       , emit: qza
    path("*.qzv")                       , emit: qzv
    path "versions.yml"                 , emit: versions_qiime2_ancombc2_asv, topic: versions

    script:
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
    def formula     = formula_in ?: "${table.baseName}"
    outfolder       = "Category-${formula}-ASV"
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime composition ancombc2 \\
        --i-table "${table}" \\
        --m-metadata-file "${metadata}" \\
        $args \\
        --p-fixed-effects-formula '${formula}' \\
        --o-ancombc2-output "ASV-${formula}.differentials.qza" \\
        --p-num-processes ${task.cpus}  \\
        --verbose
    qiime tools export \\
        --input-path "ASV-${formula}.differentials.qza" \\
        --output-path "${outfolder}/differentials"

    # Generate tabular view of ANCOMBC2 output
    qiime composition tabulate \\
        --i-data "ASV-${formula}.differentials.qza" \\
        --o-visualization "ASV-${formula}.differentials.qzv"
    qiime tools export \\
        --input-path "ASV-${formula}.differentials.qzv" \\
        --output-path "${outfolder}/differentials"

    # Generate bar plot views of ANCOMBC2 output
    qiime composition ancombc2-visualizer \\
        --i-data "ASV-${formula}.differentials.qza" \\
        --o-visualization "ASV-${formula}.visualizer.qzv"
    # 'qiime tools export' does not produce a valid html
    mv "ASV-${formula}.visualizer.qzv" "${outfolder}/"

    # Generate volcano plot
    ancombc_volcanoplot.r \\
        "${outfolder}/differentials/lfc.jsonl" \\
        "${outfolder}/differentials/q.jsonl" \\
        $args2 \\
        "${outfolder}/differentials/p.jsonl" \\
        "${outfolder}/differentials/se.jsonl" \\
        "${outfolder}/differentials/passed_ss.jsonl" \\
        "ASV-"
    mv *.volcano_plot.* "${outfolder}/"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

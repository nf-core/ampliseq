process QIIME2_ANCOMBC2_TAX {
    tag "${table.baseName}-${formula_in}-${taxlevel}"
    label 'process_medium'

    conda "${moduleDir}/envs/rachis-qiime2-linux-64-conda.yml"
    container "qiime2/qiime2:2026.4"

    input:
    tuple path(metadata), path(table), path(taxonomy), val(taxlevel), val(formula_in)

    output:
    path("visualizer/*")         , emit: plot
    path("differentials/*")      , emit: differentials
    path("*.qza")                , emit: qza, optional: true
    path("*.qzv")                , emit: qzv, optional: true
    path "versions.yml"          , emit: versions_qiime2_ancombc2_tax, topic: versions

    script:
    def args        = task.ext.args ?: ''
    def formula     = formula_in ?: "${table.baseName}"
    def prefix      = "lvl${taxlevel}-${formula}"
    def outfolder   = "Category-${formula}-level-${taxlevel}"
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # Sum data at the specified level
    qiime taxa collapse \\
        --i-table "${table}" \\
        --i-taxonomy "${taxonomy}" \\
        --p-level ${taxlevel} \\
        --o-collapsed-table "${prefix}.qza"

    # Extract summarised table and output a file with the number of taxa
    qiime tools export \\
        --input-path "${prefix}.qza" \\
        --output-path exported/
    biom convert \\
        -i exported/feature-table.biom \\
        -o "${prefix}.feature-table.tsv" \\
        --to-tsv

    if [ \$(grep -v '^#' -c "${prefix}.feature-table.tsv") -lt 2 ]; then
        mkdir differentials
        echo ${taxlevel} > differentials/\"WARNING Summing your data at taxonomic level ${taxlevel} produced less than two rows (taxa), ANCOMBC2 can't proceed -- did you specify a bad reference taxonomy?\".txt
        mkdir visualizer
        echo ${taxlevel} > visualizer/\"WARNING Summing your data at taxonomic level ${taxlevel} produced less than two rows (taxa), ANCOMBC2 can't proceed -- did you specify a bad reference taxonomy?\".txt
    else
        qiime composition ancombc2 \\
            --i-table "${prefix}.qza" \\
            --m-metadata-file "${metadata}" \\
            $args \\
            --p-fixed-effects-formula '${formula}' \\
            --o-ancombc2-output "${prefix}.differentials.qza" \\
            --p-num-processes ${task.cpus}  \\
            --verbose
        qiime tools export \\
            --input-path "${prefix}.differentials.qza" \\
            --output-path "differentials/${outfolder}"

        # Generate tabular view of ANCOMBC2 output
        qiime composition tabulate \\
            --i-data "${prefix}.differentials.qza" \\
            --o-visualization "${prefix}.differentials.qzv"
        qiime tools export \\
            --input-path "${prefix}.differentials.qzv" \\
            --output-path "differentials/${outfolder}"

        # Generate bar plot views of ANCOMBC2 output
        qiime composition ancombc2-visualizer \\
            --i-data "${prefix}.differentials.qza" \\
            --o-visualization "${prefix}.visualizer.qzv"
        # 'qiime tools export' does not produce a valid html
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

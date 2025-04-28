process QIIME2_ANCOMBC_TAX {
    tag "${table.baseName}-${formula_in}-${taxlevel}"
    label 'process_medium'
    label 'single_cpu'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple path(metadata), path(table), path(taxonomy), val(taxlevel), val(formula_in)

    output:
    path("da_barplot/*")   , emit: da_barplot
    path("differentials/*"), emit: differentials
    path("*.qza")          , emit: qza, optional: true
    path("*.qzv")          , emit: qzv, optional: true
    path "versions.yml"    , emit: versions

    script:
    def args        = task.ext.args ?: ''
    def args2       = task.ext.args2 ?: ''
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
        echo ${taxlevel} > differentials/\"WARNING Summing your data at taxonomic level ${taxlevel} produced less than two rows (taxa), ANCOMBC can't proceed -- did you specify a bad reference taxonomy?\".txt
        mkdir da_barplot
        echo ${taxlevel} > da_barplot/\"WARNING Summing your data at taxonomic level ${taxlevel} produced less than two rows (taxa), ANCOMBC can't proceed -- did you specify a bad reference taxonomy?\".txt
    else
        qiime composition ancombc \\
            --i-table "${prefix}.qza" \\
            --m-metadata-file "${metadata}" \\
            $args \\
            --p-formula '${formula}' \\
            --o-differentials "${prefix}.differentials.qza" \\
            --verbose
        qiime tools export \\
            --input-path "${prefix}.differentials.qza" \\
            --output-path "differentials/${outfolder}"

        # Generate tabular view of ANCOM-BC output
        qiime composition tabulate \\
            --i-data "${prefix}.differentials.qza" \\
            --o-visualization "${prefix}.differentials.qzv"
        qiime tools export \\
            --input-path "${prefix}.differentials.qzv" \\
            --output-path "differentials/${outfolder}"

        # Generate bar plot views of ANCOM-BC output
        qiime composition da-barplot \\
            --i-data "${prefix}.differentials.qza" \\
            $args2 \\
            --o-visualization "${prefix}.da_barplot.qzv"
        qiime tools export --input-path "${prefix}.da_barplot.qzv" \\
            --output-path "da_barplot/${outfolder}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

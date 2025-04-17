process QIIME2_ANCOM_TAX {
    tag "${table.baseName}-taxonomic_level:${taxlevel}"
    label 'process_medium'
    label 'single_cpu'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    tuple path(metadata), path(table), path(taxonomy) ,val(taxlevel)

    output:
    path "ancom/*"      , emit: ancom
    path "versions.yml" , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"
    mkdir ancom

    # Sum data at the specified level
    qiime taxa collapse \\
        --i-table ${table} \\
        --i-taxonomy ${taxonomy} \\
        --p-level ${taxlevel} \\
        --o-collapsed-table lvl${taxlevel}-${table}

    # Extract summarised table and output a file with the number of taxa
    qiime tools export \\
        --input-path lvl${taxlevel}-${table} \\
        --output-path exported/
    biom convert \\
        -i exported/feature-table.biom \\
        -o ${table.baseName}-level-${taxlevel}.feature-table.tsv \\
        --to-tsv

    if [ \$(grep -v '^#' -c ${table.baseName}-level-${taxlevel}.feature-table.tsv) -lt 2 ]; then
        echo ${taxlevel} > ancom/\"WARNING ${table.baseName} Summing your data at taxonomic level ${taxlevel} produced less than two rows (taxa), ANCOM can't proceed -- did you specify a bad reference taxonomy?\".txt
    else
        qiime composition add-pseudocount \\
                --i-table lvl${taxlevel}-${table} \\
                --o-composition-table comp-lvl${taxlevel}-${table}
        qiime composition ancom \\
                --i-table comp-lvl${taxlevel}-${table} \\
                --m-metadata-file ${metadata} \\
                --m-metadata-column ${table.baseName} \\
                --o-visualization comp-lvl${taxlevel}-${table.baseName}.qzv
        qiime tools export --input-path comp-lvl${taxlevel}-${table.baseName}.qzv \\
                --output-path ancom/Category-${table.baseName}-level-${taxlevel}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

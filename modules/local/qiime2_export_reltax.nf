process QIIME2_EXPORT_RELTAX {
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(table)
    path(taxonomy)
    val(tax_agglom_min)
    val(tax_agglom_max)

    output:
    path("*.tsv")        , emit: tsv
    path "versions.yml"  , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    ##on several taxa level
    array=(\$(seq ${tax_agglom_min} 1 ${tax_agglom_max}))

    for i in \${array[@]}
    do
        #collapse taxa
        qiime taxa collapse \\
            --i-table ${table} \\
            --i-taxonomy ${taxonomy} \\
            --p-level \$i \\
            --o-collapsed-table table-\$i.qza
        #convert to relative abundances
        qiime feature-table relative-frequency \\
            --i-table table-\$i.qza \\
            --o-relative-frequency-table relative-table-\$i.qza
        #export to biom
        qiime tools export \\
            --input-path relative-table-\$i.qza \\
            --output-path relative-table-\$i
        #convert to tab separated text file
        biom convert \\
            -i relative-table-\$i/feature-table.biom \\
            -o rel-table-\$i.tsv --to-tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

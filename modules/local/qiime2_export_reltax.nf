process QIIME2_EXPORT_RELTAX {
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

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
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    ##on several taxa level
    array=(\$(seq ${tax_agglom_min} 1 ${tax_agglom_max}))

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
    END_VERSIONS
    """
}

process QIIME2_EXPORT_RELTAX {
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(table)
    path(taxonomy)
    val(tax_agglom_min)
    val(tax_agglom_max)

    output:
    path("*.tsv")        , emit: tsv
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
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

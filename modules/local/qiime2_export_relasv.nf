process QIIME2_EXPORT_RELASV {
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(table)

    output:
    path("rel-table-ASV.tsv"), emit: tsv
    path "versions.yml"      , emit: versions

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

    #convert to relative abundances
    qiime feature-table relative-frequency \\
        --i-table ${table} \\
        --o-relative-frequency-table relative-table-ASV.qza

    #export to biom
    qiime tools export \\
        --input-path relative-table-ASV.qza \\
        --output-path relative-table-ASV

    #convert to tab separated text file "rel-table-ASV.tsv"
    biom convert \\
        -i relative-table-ASV/feature-table.biom \\
        -o rel-table-ASV.tsv --to-tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

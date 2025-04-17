process QIIME2_EXPORT_RELASV {
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(table)

    output:
    path("rel-table-ASV.tsv"), emit: tsv
    path "versions.yml"      , emit: versions

    script:
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

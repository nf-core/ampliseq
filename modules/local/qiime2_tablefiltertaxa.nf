process QIIME2_TABLEFILTERTAXA {
    tag "taxa:${exclude_taxa};min-freq:${min_frequency};min-samples:${min_samples}"
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(table)
    path(taxonomy)
    val(min_frequency)
    val(min_samples)
    val(exclude_taxa)

    output:
    path("filtered-table.qza"), emit: qza
    path("filtered-table.tsv"), emit: tsv
    path "versions.yml"       , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    if ! [ \"${exclude_taxa}\" = \"none\" ]; then
        qiime taxa filter-table \\
            --i-table ${table} \\
            --i-taxonomy ${taxonomy} \\
            --p-exclude "${exclude_taxa}" \\
            --p-mode contains \\
            --o-filtered-table tax_filtered-table.qza
        filtered_table="tax_filtered-table.qza"
    else
        filtered_table=${table}
    fi

    qiime feature-table filter-features \\
        --i-table \$filtered_table \\
        --p-min-frequency ${min_frequency} \\
        --p-min-samples ${min_samples} \\
        --o-filtered-table filtered-table.qza

    #produce raw count table in biom format "table/feature-table.biom"
    qiime tools export \\
        --input-path filtered-table.qza  \\
        --output-path table
    #produce raw count table
    biom convert \\
        -i table/feature-table.biom \\
        -o table/feature-table.tsv  \\
        --to-tsv
    cp table/feature-table.tsv filtered-table.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

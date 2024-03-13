process QIIME2_TABLEFILTERTAXA {
    tag "taxa:${exclude_taxa};min-freq:${min_frequency};min-samples:${min_samples}"
    label 'process_low'

    container "qiime2/core:2023.7"

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

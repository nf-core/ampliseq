process QIIME2_DIVERSITY_CORE {
    label 'process_low'

    container "quay.io/qiime2/core:2022.11"

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }

    input:
    path(metadata)
    path(table)
    path(tree)
    path(stats)
    val(mindepth)

    output:
    path("diversity_core/*_pcoa_results.qza")   , emit: pcoa
    path("diversity_core/*_vector.qza")         , emit: vector
    path("diversity_core/*_distance_matrix.qza"), emit: distance
    path "versions.yml"                         , emit: versions
    path("*rarefaction.txt")                    , emit: depth

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    mindepth=\$(count_table_minmax_reads.py $stats minimum 2>&1)
    if [ \"\$mindepth\" -lt \"$mindepth\" ]; then mindepth=$mindepth; fi
    if [ \"\$mindepth\" -gt \"10000\" ]; then echo \$mindepth >\"Use the sampling depth of \$mindepth for rarefaction.txt\" ; fi
    if [ \"\$mindepth\" -lt \"10000\" -a \"\$mindepth\" -gt \"5000\" ]; then echo \$mindepth >\"WARNING The sampling depth of \$mindepth is quite small for rarefaction.txt\" ; fi
    if [ \"\$mindepth\" -lt \"5000\" -a \"\$mindepth\" -gt \"1000\" ]; then echo \$mindepth >\"WARNING The sampling depth of \$mindepth is very small for rarefaction.txt\" ; fi
    if [ \"\$mindepth\" -lt \"1000\" ]; then echo \$mindepth >\"WARNING The sampling depth of \$mindepth seems too small for rarefaction.txt\" ; fi

    qiime diversity core-metrics-phylogenetic \\
        --m-metadata-file ${metadata} \\
        --i-phylogeny ${tree} \\
        --i-table ${table} \\
        --p-sampling-depth \$mindepth \\
        --output-dir diversity_core \\
        --p-n-jobs-or-threads ${task.cpus} \\
        --verbose

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

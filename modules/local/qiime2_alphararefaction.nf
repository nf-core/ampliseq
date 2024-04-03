process QIIME2_ALPHARAREFACTION {
    label 'process_low'

    container "qiime2/core:2023.7"

    input:
    path(metadata)
    path(table)
    path(tree)
    path(stats)

    output:
    path("alpha-rarefaction/*"), emit: rarefaction
    path("*.qzv")              , emit: qzv
    path "versions.yml"        , emit: versions

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

    maxdepth=\$(count_table_minmax_reads.py $stats maximum 2>&1)

    #check values
    if [ \"\$maxdepth\" -gt \"75000\" ]; then maxdepth=\"75000\"; fi
    if [ \"\$maxdepth\" -gt \"5000\" ]; then maxsteps=\"250\"; else maxsteps=\$((maxdepth/20)); fi
    qiime diversity alpha-rarefaction  \\
        --i-table ${table}  \\
        --i-phylogeny ${tree}  \\
        --p-max-depth \$maxdepth  \\
        --m-metadata-file ${metadata}  \\
        --p-steps \$maxsteps  \\
        --p-iterations 10  \\
        --o-visualization alpha-rarefaction.qzv
    qiime tools export --input-path alpha-rarefaction.qzv  \\
        --output-path alpha-rarefaction

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

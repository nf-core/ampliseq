process QIIME2_ALPHARAREFACTION {
    label 'process_low'

    conda "${projectDir}/modules/local/envs/qiime2-amplicon-2024.10-py310-linux-conda.yml"
    container "qiime2/amplicon:2024.10"

    input:
    path(metadata)
    path(table)
    path(tree)
    path(stats)

    output:
    path("alpha-rarefaction/*"), emit: rarefaction
    path("*.qzv")              , emit: qzv
    path "versions.yml"        , emit: versions

    script:
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

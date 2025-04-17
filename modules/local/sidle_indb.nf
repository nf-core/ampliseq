process SIDLE_INDB {
    label 'process_single'

    conda "${projectDir}/modules/local/envs/pipesidle-0-1-0-beta.yml"
    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    path(seq)
    path(tax)

    output:
    path("db_sequences.qza"), emit: seq
    path("db_taxonomy.qza") , emit: tax
    path "versions.yml"     , emit: versions

    script:
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # db_seq
    qiime tools import \\
        --input-path $seq \\
        --output-path db_sequences.qza \\
        --type 'FeatureData[Sequence]'

    # db_tax
    qiime tools import \\
        --input-path $tax \\
        --output-path db_taxonomy.qza \\
        --type 'FeatureData[Taxonomy]' \\
        --input-format HeaderlessTSVTaxonomyFormat

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

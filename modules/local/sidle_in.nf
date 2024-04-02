process SIDLE_IN {
    tag "$meta.region"
    label 'process_single'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    tuple val(meta), path(table), path(seq)

    output:
    tuple val(meta), path("*_table.qza"), path("*_rep-seqs.qza"), emit: table_seq
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sidle in QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.region}"
    """
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    # seq
    qiime tools import \\
        --input-path "$seq" \\
        --type 'FeatureData[Sequence]' \\
        --output-path ${prefix}_rep-seqs.qza

    # table
    biom convert -i "$table" -o table.biom --table-type="OTU table" --to-hdf5
    qiime tools import \\
        --input-path table.biom \\
        --type 'FeatureTable[Frequency]' \\
        --input-format BIOMV210Format \\
        --output-path ${prefix}_table.qza

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
    END_VERSIONS
    """
}

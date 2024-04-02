process SIDLE_TREERECON {
    label 'process_medium'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    path(reconstruction_fragments)
    path(ref_db_tree)

    output:
    path("reconstructed_tree.qza")       , emit: qza
    path("reconstruction_placements.qza"), emit: qza_placements
    path("reconstructed_tree.nwk")       , emit: nwk
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sidle in QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    # https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#reconstructing-the-phylogenetic-tree
    # required: SEPP file https://forum.qiime2.org/t/sidle-tutorial-missing-aligned-sequence-file/20604/8
    # SEPP file only available for Greengenes 13_8 or SILVE 128 (not 138!): https://forum.qiime2.org/t/error-in-reconstructing-the-phylogenetic-tree/23757/8
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime fragment-insertion sepp \\
        --p-threads $task.cpus \\
        --i-representative-sequences $reconstruction_fragments \\
        --i-reference-database $ref_db_tree \\
        --o-tree reconstructed_tree.qza \\
        --o-placements reconstruction_placements.qza

    #export tree file
    qiime tools export \\
        --input-path reconstructed_tree.qza \\
        --output-path exported
    cp exported/tree.nwk reconstructed_tree.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        q2-fragment-insertion: \$( qiime fragment-insertion --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

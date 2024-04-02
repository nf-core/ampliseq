process SIDLE_SEQRECON {
    label 'process_medium'
    label 'single_cpu'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    path(reconstruction_map)
    path(reconstruction_summary)
    path(db_aligned_sequences)

    output:
    path("reconstruction_fragments.qza") , emit: qza
    path("reconstruction_fragments/*")   , emit: visualisation
    path("reconstructed_fragments.fasta"), emit: fasta
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sidle in QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#reconstructing-the-phylogenetic-tree
    #https://forum.qiime2.org/t/sidle-tutorial-missing-aligned-sequence-file/20604/4 for db_aligned_sequences
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    #CPU=1
    qiime sidle reconstruct-fragment-rep-seqs \\
        --i-reconstruction-map ${reconstruction_map} \\
        --i-reconstruction-summary ${reconstruction_summary} \\
        --i-aligned-sequences ${db_aligned_sequences} \\
        --o-representative-fragments reconstruction_fragments.qza

    #export visualisation
    qiime metadata tabulate \\
        --m-input-file reconstruction_fragments.qza \\
        --o-visualization reconstruction_fragments.qzv
    qiime tools export \\
        --input-path reconstruction_fragments.qzv \\
        --output-path "reconstruction_fragments"

    #export fasta file
    qiime tools export \\
        --input-path reconstruction_fragments.qza \\
        --output-path exported
    cp exported/dna-sequences.fasta reconstructed_fragments.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle: \$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

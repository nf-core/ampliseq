process SIDLE_TAXRECON {
    label 'process_single'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    path(reconstruction_map)
    path(tax)

    output:
    path("reconstruction_taxonomy.qza"), emit: qza
    path("reconstruction_taxonomy/*")  , emit: visualisation
    path("reconstruction_taxonomy.tsv"), emit: tsv
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sidle in QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#taxonomic-reconstruction
    #https://forum.qiime2.org/t/sidle-reconstruct-database/25439
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime sidle reconstruct-taxonomy \\
        --i-reconstruction-map ${reconstruction_map} \\
        --i-taxonomy ${tax} \\
        $args \\
        --o-reconstructed-taxonomy reconstruction_taxonomy.qza

    #export visualisation
    qiime metadata tabulate \\
        --m-input-file reconstruction_taxonomy.qza \\
        --o-visualization reconstruction_taxonomy.qzv
    qiime tools export \\
        --input-path reconstruction_taxonomy.qzv \\
        --output-path "reconstruction_taxonomy"

    #export taxonomic tsv
    qiime tools export \\
        --input-path reconstruction_taxonomy.qza \\
        --output-path exported
    cp exported/taxonomy.tsv reconstruction_taxonomy.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle: \$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

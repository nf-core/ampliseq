process SIDLE_DBRECON {
    label 'process_medium'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    val(metaid)
    path(map)
    path(aligned_map)

    output:
    path("reconstruction_map.qza")    , emit: reconstruction_map
    path("reconstruction_summary.qza"), emit: reconstruction_summary
    path("reconstruction_summary/*")  , emit: visualisation
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sidle in QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def db_input = ""
    // sort the input so that the regions are sorted by sequence
    def df = [metaid, map, aligned_map].transpose().sort{ it[0] }
    for (i in df) {
        db_input += " --p-region "+i[0]+" --i-kmer-map "+i[1]+" --i-regional-alignment "+i[2]
    }
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#database-reconstruction
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime sidle reconstruct-database \\
        --p-n-workers $task.cpus \\
        $db_input \\
        $args \\
        --o-database-map reconstruction_map.qza \\
        --o-database-summary reconstruction_summary.qza

    #database summary can be used to evaluate the quality of the reconstruction; see Fuks, C; Elgart, M; Amir, A; et al (2018) “Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling.” Microbiome. 6:17. doi: 10.1186/s40168-017-0396-x
    qiime metadata tabulate \\
        --m-input-file reconstruction_summary.qza \\
        --o-visualization reconstruction_summary.qzv
    qiime tools export \\
        --input-path reconstruction_summary.qzv \\
        --output-path "reconstruction_summary"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle: \$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

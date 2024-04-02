process SIDLE_TABLERECON {
    label 'process_medium'

    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    val(metaid)
    path(table)
    path(aligned_map)
    path(reconstruction_map)
    path(reconstruction_summary)

    output:
    path("reconstruction_table.qza")        , emit: qza
    path("reconstruction_table/*")          , emit: exported
    path("reconstructed_feature-table.biom"), emit: biom
    path("reconstructed_feature-table.tsv") , emit: tsv
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "Sidle in QIIME2 does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def region_input = ""
    // sort the input so that the regions are sorted by sequence
    def df = [metaid, aligned_map, table].transpose().sort{ it[0] }
    for (i in df) {
        region_input += " --p-region "+i[0]+" --i-regional-alignment "+i[1]+" --i-regional-table "+i[2]
    }
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#table-reconstruction
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    qiime sidle reconstruct-counts \\
        --p-n-workers $task.cpus \\
        $region_input \\
        --i-database-map $reconstruction_map \\
        --i-database-summary $reconstruction_summary \\
        $args \\
        --o-reconstructed-table reconstruction_table.qza

    #export visualisation
    qiime feature-table summarize \\
        --i-table reconstruction_table.qza \\
        --o-visualization reconstruction_table.qzv
    qiime tools export \\
        --input-path reconstruction_table.qzv \\
        --output-path "reconstruction_table"

    #export feature table in biom and tsv format
    qiime tools export \\
        --input-path reconstruction_table.qza \\
        --output-path exported
    biom convert \\
        -i exported/feature-table.biom \\
        -o reconstructed_feature-table.tsv \\
        --to-tsv
    cp exported/feature-table.biom reconstructed_feature-table.biom

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qiime2: \$( qiime --version | sed '1!d;s/.* //' )
        qiime2 plugin sidle: \$( qiime sidle --version | sed 's/ (.*//' | sed 's/.*version //' )
        q2-sidle: \$( qiime sidle --version | sed 's/.*version //' | sed 's/)//' )
    END_VERSIONS
    """
}

process SIDLE_TABLERECON {
    label 'process_medium'

    conda "${moduleDir}/envs/pipesidle-0-1-0-beta.yml"
    container 'nf-core/pipesidle:0.1.0-beta'

    input:
    val(region)
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

    script:
    def args = task.ext.args ?: ''
    def region_input = ""
    // input must be sorted already by regions
    def df = [region, aligned_map, table].transpose()
    df.each { i ->
        region_input += " --p-region "+i[0]+" --i-regional-alignment "+i[1]+" --i-regional-table "+i[2]
    }
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#table-reconstruction
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"

    #understand input args
    echo "Args:"
    echo "$args"
    echo "$aligned_map"
    echo "___"

    #visualize arguments
    echo "$args" > args.txt
    cat args.txt
    cat "$region_input"
    # Print original and filtered tables for all region tables in $aligned_map for investigation
    for region_table in ${aligned_map}; do
        # Print original table
        echo "Original table for \$region_table:"
        region_table_base=\$(basename "\$region_table" .qza)
        qiime tools export \
            --input-path "\$region_table" \
            --output-path "\${region_table_base}_exported"
        cat "\${region_table_base}_exported/feature-table.tsv"

        # Filter zero-count samples and print filtered table
        filtered_table="\${region_table_base}.filtered.qza"
        filtered_table_exported_folder="\${region_table_base}_filtered"
        qiime feature-table filter-samples \
            --i-table "\$region_table" \
            --p-min-frequency 1 \
            --o-filtered-table "\$filtered_table"
        echo "Filtered table for \$region_table:"
        qiime tools export \
            --input-path "\$filtered_table" \
            --output-path "\${filtered_table_exported_folder}"
        cat "\${filtered_table_exported_folder}/feature-table.tsv"
    done



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

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

    path("*_feature-table.tsv")             , emit: regional_tsv
    path("total_counts.tsv")                , emit: total_counts
    path("kept_samples.tsv")                , emit: kept_samples

    script:
    def args = task.ext.args ?: ''
    def region_input = ""
    // input must be sorted already by regions
    def df = [region, aligned_map, table].transpose()

    def minCountsMatcher = args =~ /(^|\s)--p-min-counts(?:=|\s+)(\d+)(?=\s|$)/
    def minCounts = minCountsMatcher.find() ? minCountsMatcher.group(2) : '0'

    df.each { i ->
        def table_base = i[2].toString().replaceAll(/\.qza$/, '')
        region_input += " --p-region ${i[0]} --i-regional-alignment ${i[1]} --i-regional-table ${table_base}.filtered.qza"
    }
    """
    #https://q2-sidle.readthedocs.io/en/latest/reconstruction.html#table-reconstruction
    export XDG_CONFIG_HOME="./xdgconfig"
    export MPLCONFIGDIR="./mplconfigdir"
    export NUMBA_CACHE_DIR="./numbacache"


    # Export original and filtered regional tables as TSVs in the task work dir.
    for region_table in ${table}; do
        region_table_base=\$(basename "\$region_table" .qza)
        original_table_exported_folder="\${region_table_base}_exported"
        original_table_tsv="\${region_table_base}_feature-table.tsv"

        qiime tools export \
            --input-path "\$region_table" \
            --output-path "\${original_table_exported_folder}"
        biom convert \
            -i "\${original_table_exported_folder}/feature-table.biom" \
            -o "\${original_table_tsv}" \
            --to-tsv
    done

    # Determine total counts per sample across all regional tables, to be used as input for reconstruction
    {
        echo -e "sample-id\ttotal_count"
        awk '
        BEGIN { FS=OFS="\t" }

        \$1 == "#OTU ID" {
            delete cols
            for (i=2; i<=NF; i++) cols[i] = \$i
            next
        }

        /^#/ { next }

        {
            for (i=2; i<=NF; i++) sum[cols[i]] += \$i
        }

        END {
            for (s in sum) {
                if (s != "") print s, sum[s]
            }
        }
        ' *_feature-table.tsv | sort -k1,1
    } > total_counts.tsv 

    awk -v min_counts="${minCounts}" '
    BEGIN { FS=OFS="\t" }
    NR==1 { print "sample-id"; next }
    \$2 >= min_counts { print \$1 }
    ' total_counts.tsv > kept_samples.tsv


    for region_table in ${table}; do
        region_table_base=\$(basename "\$region_table" .qza)
        qiime feature-table filter-samples \\
            --i-table "\$region_table" \\
            --m-metadata-file kept_samples.tsv \\
            --o-filtered-table "\${region_table_base}.filtered.qza"
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

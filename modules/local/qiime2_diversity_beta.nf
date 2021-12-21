process QIIME2_DIVERSITY_BETA {
    tag "${core.baseName}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    tuple path(metadata), path(core), val(category)

    output:
    path("beta_diversity/*"), emit: beta
    path "versions.yml"     , emit: versions

    script:
    if ( category.length() > 0 ) {
        """
        export XDG_CONFIG_HOME="\${PWD}/HOME"

        IFS=',' read -r -a metacategory <<< \"$category\"
        for j in \"\${metacategory[@]}\"
        do
            qiime diversity beta-group-significance \
                --i-distance-matrix ${core} \
                --m-metadata-file ${metadata} \
                --m-metadata-column \"\$j\" \
                --o-visualization ${core.baseName}-\$j.qzv \
                --p-pairwise
            qiime tools export --input-path ${core.baseName}-\$j.qzv \
                --output-path beta_diversity/${core.baseName}-\$j
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
        END_VERSIONS
        """
    } else {
        """
        mkdir beta_diversity
        echo "" > "beta_diversity/WARNING No column in ${metadata.baseName} seemed suitable.txt"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
        END_VERSIONS
        """
    }
}

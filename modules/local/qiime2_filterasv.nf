process QIIME2_FILTERASV {
    tag "${category}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    path(metadata)
    path(table)
    val(category)

    output:
    path("*.qza")       , emit: qza
    path "versions.yml" , emit: versions

    script:
    if ( category.length() > 0 ) {
        """
        export XDG_CONFIG_HOME="\${PWD}/HOME"

        IFS=',' read -r -a metacategory <<< \"$category\"

        #remove samples that do not have any value
        for j in \"\${metacategory[@]}\"
        do
            qiime feature-table filter-samples \
                --i-table ${table} \
                --m-metadata-file ${metadata} \
                --p-where \"\$j<>\'\'\" \
                --o-filtered-table \$j.qza
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
        END_VERSIONS
        """
    } else {
        """
        mkdir beta_diversity
        echo "" > "WARNING No column in ${metadata.baseName} seemed suitable.qza"
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
        END_VERSIONS
        """
    }
}

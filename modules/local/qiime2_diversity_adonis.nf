process QIIME2_DIVERSITY_ADONIS {
    tag "${core.baseName}"
    label 'process_low'

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.8"

    input:
    tuple path(metadata), path(core), val(category)

    output:
    path("adonis/*")     , emit: html
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    def formula = params.qiime_adonis_formula ?: ''
    if ( category.length() > 0 || params.qiime_adonis_formula ) {
        """
        export XDG_CONFIG_HOME="\${PWD}/HOME"

        if [ "${formula}" == '' ]; then
            adonisformula=\$( echo $category | sed "s/,/+/g" )
        else
            adonisformula=${formula}
        fi

        qiime diversity adonis \\
            --p-n-jobs $task.cpus \\
            --i-distance-matrix ${core} \\
            --m-metadata-file ${metadata} \\
            --o-visualization ${core.baseName}_adonis.qzv \\
            $args \\
            --p-formula \$adonisformula
        qiime tools export --input-path ${core.baseName}_adonis.qzv \\
            --output-path adonis/${core.baseName}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
        END_VERSIONS
        """
    } else {
        """
        mkdir adonis
        echo "" > "adonis/WARNING No formula was given with --qiime_adonis_formula and no column in ${metadata.baseName} seemed suitable.txt"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            qiime2: \$( qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g" )
        END_VERSIONS
        """
    }
}

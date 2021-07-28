// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_FILTERASV {
    tag "${category}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(metadata)
    path(table)
    val(category)

    output:
    path("*.qza")       , emit: qza
    path "*.version.txt", emit: version

    script:
    def software     = getSoftwareName(task.process)
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

        echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
        """
    } else {
        """
        mkdir beta_diversity
        echo "" > "WARNING No column in ${metadata.baseName} seemed suitable.qza"
        echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
        """
    }
}

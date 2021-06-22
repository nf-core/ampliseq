// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_DIVERSITY_ALPHA {
    tag "${core.baseName}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    tuple path(metadata), path(core), val(category)

    output:
    path("alpha_diversity/*"), emit: alpha
    path "*.version.txt"     , emit: version

    script:
    def software     = getSoftwareName(task.process)
    if ( category.length() > 0 ) {
        """
        export XDG_CONFIG_HOME="\${PWD}/HOME"

        qiime diversity alpha-group-significance \
            --i-alpha-diversity ${core} \
            --m-metadata-file ${metadata} \
            --o-visualization ${core.baseName}-vis.qzv
        qiime tools export --input-path ${core.baseName}-vis.qzv \
            --output-path "alpha_diversity/${core.baseName}"

        echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
        """
    } else {
        """
        mkdir alpha_diversity
        echo "" > "alpha_diversity/WARNING No column in ${metadata.baseName} seemed suitable.txt"
        echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
        """        
    }
}
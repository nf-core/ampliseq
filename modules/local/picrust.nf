// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process PICRUST {
    tag "${seq},${abund}"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::picrust2=2.4.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/picrust2:2.4.1--py_0"
    } else {
        container "quay.io/biocontainers/picrust2:2.4.1--py_0"
    }

    input:
    path(seq)
    path(abund)
    val(source)
    val(message)

    output:
    path("all_output/*") , emit: outfolder
    path("*_descrip.tsv"), emit: pathways
    path "*.version.txt" , emit: version
    path "*.args.txt"    , emit: args
    path "${message}.txt"

    script:
    def software      = getSoftwareName(task.process)
    """
    #If input is QIIME2 file, than (1) the first line and (2) the first character (#) of the second line need to be removed
    if [ "$source" == 'QIIME2' ]
    then
        tail -n +2 "$abund" > "${abund}.tmp" && mv "${abund}.tmp" "$abund"
    fi

    picrust2_pipeline.py \\
        $options.args \\
        -s $seq \\
        -i $abund \\
        -o all_output \\
        -p $task.cpus \\
        --in_traits EC,KO \\
        --verbose

    #Add descriptions to identifiers
    add_descriptions.py -i all_output/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_pred_metagenome_unstrat_descrip.tsv
    add_descriptions.py -i all_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_pred_metagenome_unstrat_descrip.tsv
    add_descriptions.py -i all_output/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o METACYC_path_abun_unstrat_descrip.tsv

    echo "$message" > "${message}.txt"
    echo -e "picrust\t$options.args" > "picrust.args.txt"
    picrust2_pipeline.py -v | sed -e "s/picrust2_pipeline.py //g" > "${software}.version.txt"
    """
}

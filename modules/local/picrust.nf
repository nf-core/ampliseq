process PICRUST {
    tag "${seq},${abund}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::picrust2=2.4.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picrust2:2.4.2--pyhdfd78af_0' :
        'quay.io/biocontainers/picrust2:2.4.2--pyhdfd78af_0' }"

    input:
    path(seq)
    path(abund)
    val(source)
    val(message)

    output:
    path("all_output/*") , emit: outfolder
    path("*_descrip.tsv"), emit: pathways
    path "versions.yml"  , emit: versions
    path "*.args.txt"    , emit: args
    path "${message}.txt"

    script:
    def args = task.ext.args ?: ''
    """
    #If input is QIIME2 file, than (1) the first line and (2) the first character (#) of the second line need to be removed
    if [ "$source" == 'QIIME2' ]
    then
        tail -n +2 "$abund" > "${abund}.tmp" && mv "${abund}.tmp" "$abund"
    fi

    picrust2_pipeline.py \\
        $args \\
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
    echo -e "picrust\t$args" > "picrust.args.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        picrust2: \$( picrust2_pipeline.py -v | sed -e "s/picrust2_pipeline.py //g" )
    END_VERSIONS
    """
}

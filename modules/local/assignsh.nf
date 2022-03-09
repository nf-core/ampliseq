process ASSIGNSH {
    tag "${asvtable}"
    label 'process_low'
    
    conda (params.enable_conda ? "pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5':
        'quay.io/biocontainers/pandas:1.1.5' }"

    input:
    path asvtable
    val sh_info
    path blastfile
    val  outtable
    
    output:
    path outtable        , emit: tsv
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    add_sh_to_taxonomy.py ${sh_info.join(' ')} $asvtable $blastfile $outtable $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
    END_VERSIONS
    """
}
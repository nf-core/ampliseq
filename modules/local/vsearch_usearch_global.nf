process VSEARCH_USEARCH_GLOBAL {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::vsearch=2.21.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--h95f258a_0':
        'quay.io/biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    path asvfasta
    path dbfasta
    val  outfile

    output:
    path outfile         , emit: tsv
    path "versions.yml"  , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    vsearch \\
        --usearch_global $asvfasta \\
        --db $dbfasta \\
        --threads $task.cpus \\
        $args \\
        --blast6out $outfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}

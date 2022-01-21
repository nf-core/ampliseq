process ITSX_CUTASV {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::itsx=1.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/itsx:1.1.3--hdfd78af_1' :
        'quay.io/biocontainers/itsx:1.1.3--hdfd78af_1' }"

    input:
    path fasta

    output:
    path "ASV_ITS_seqs.full.fasta", emit: fasta
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    ITSx \\
        -i $fasta \\
        $args \\
        --cpu $task.cpus \\
        -o ASV_ITS_seqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ITSx: \$( ITSx -h 2>&1 > /dev/null | tail -n 2 | head -n 1 | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}

process EMU_ABUNDANCE {
    debug true
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::emu=3.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emu:3.4.4--hdfd78af_1':
        'quay.io/biocontainers/emu:3.4.4--hdfd78af_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*abundance.tsv"), emit: report
    tuple val(meta), path("*read-assignment-distributions.tsv"), emit: assignment_report, optional:true
    path "versions.yml"           , emit: versions
    tuple val(meta), path("*.sam"), emit: samfile, optional:true
    tuple val(meta), path("*.fa"), emit: unclassified_fa , optional:true


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    emu \\
        abundance \\
        $args \\
        --threads $task.cpus \\
        $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emu: \$(echo \$(emu --version 2>&1) | sed 's/^.*emu //; s/Using.*\$//' )
    END_VERSIONS
    """
}

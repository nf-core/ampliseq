process COMPARE_PROFILE {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(observed_profile, stageAs: 'observed/*') // Exported taxonomic profile table from within pipeline
    path(expected_profile, stageAs: 'expected/*')

    output:
    path("*.svg")                        , emit: svg
    path("*.png")                        , emit: png
    path("profile_summary.tsv")          , emit: performance
    path("profile_per-sample.tsv")       , emit: performance_per_sample
    path("profile.log")                  , emit: log
    path "versions.yml"                  , emit: versions, topic: versions

    script:
    def fbeta = task.ext.fbeta ?: 2
    def mode = task.ext.mode ?: "complete"
    """
    compare_profile.r \\
        "${meta.id}" \\
        "$observed_profile" \\
        "$expected_profile" \\
        "$fbeta" \\
        "$mode" \\
        >"profile.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
    END_VERSIONS
    """
}

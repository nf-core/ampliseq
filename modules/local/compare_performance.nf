process COMPARE_PERFORMANCE {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(alignment) // alignment file with "query,target,mismatch_final"
    path(observed_abundances)        // ASV table from within pipeline
    path(expected_abundances)        // ASV table from params.expected_abundances

    output:
    path("*.svg")                            , emit: svg
    path("*.png")                            , emit: png
    path("performance_summary.tsv")          , emit: performance
    path("performance_per-sample.tsv")       , emit: performance_per_sample
    path("*_abundances.tsv")                 , emit: abundance_table_per_sample
    path("performance.log")                  , emit: log
    path("Warnings.txt")                     , emit: warnings
    path "versions.yml"                      , emit: versions, topic: versions

    script:
    def fbeta = task.ext.fbeta ?: 2
    def mismatch_threshold = task.ext.mismatch_threshold ?: 0
    """
    compare_performance.r \\
        "${meta.id}" \\
        "$alignment" \\
        "$observed_abundances" \\
        "$expected_abundances" \\
        "$fbeta" \\
        "$mismatch_threshold" \\
        >"performance.log"

    # isolate warnings (when grep find no match, exit code is 1, "|| true" fixes that)
    grep "WARN -" "performance.log" | sed 's/^\\[1\\] //g' | sed 's/"//g' > Warnings.txt || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
    END_VERSIONS
    """
}

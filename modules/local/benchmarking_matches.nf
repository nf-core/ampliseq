process BENCHMARKING_MATCHES {
    tag "${prefix}"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(blast6out) // VSEARCH_USEARCHGLOBAL
    path(detected_abundances)        // ASV table from within pipeline
    path(expected_abundances)        // ASV table from params.benchmark_abundances
    val(similarity_threshold)        // Similarity threshold for maches of VSEARCH
    val(query_or_target)             // what reagion to evaluate (query,target,alignment) from params.benchmark_region

    output:
    path("*.svg")                            , emit: svg
    path("*.png")                            , emit: png
    path("*_nucleotide-differences.tsv")     , emit: matches
    path("*_per-sample.tsv")                 , emit: matches_per_sample
    path("*nucleotide-differences.log")      , emit: log
    path("*.md5sum_version")                 , emit: md5sum_version
    path "versions.yml"                      , emit: versions, topic: versions

    script:
    prefix = "${meta.id}"
    """
    benchmarking_matches.r \\
        "$blast6out" \\
        "$detected_abundances" \\
        "$expected_abundances" \\
        "$prefix" \\
        "$similarity_threshold" \\
        "$query_or_target" \\
        >"${prefix}_nucleotide-differences.log"
    echo "md5sum_version $prefix" > "${prefix}.md5sum_version"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
    END_VERSIONS
    """
}

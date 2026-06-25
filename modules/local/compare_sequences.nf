process COMPARE_SEQUENCES {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::r-base=4.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.2.1' :
        'biocontainers/r-base:4.2.1' }"

    input:
    tuple val(meta), path(blast6out) // VSEARCH_USEARCHGLOBAL
    path(detected_abundances)        // ASV table from within pipeline
    path(expected_abundances)        // ASV table from params.expected_abundances
    val(similarity_threshold)        // Similarity threshold for maches of VSEARCH
    val(query_or_target)             // what reagion to evaluate (query,target,alignment) from params.expected_sequences_region

    output:
    path("*.svg")                            , emit: svg
    path("*.png")                            , emit: png
    path("*nucleotide-differences.tsv")      , emit: matches
    path("*_per-sample.tsv")                 , emit: matches_per_sample
    path("*nucleotide-differences.log")      , emit: log
    path("md5sum_version.txt")               , emit: md5sum_version
    path("Warnings.txt")                     , emit: warnings
    path "versions.yml"                      , emit: versions, topic: versions

    script:
    """
    compare_sequences.r \\
        "$blast6out" \\
        "$detected_abundances" \\
        "$expected_abundances" \\
        "" \\
        "$similarity_threshold" \\
        "$query_or_target" \\
        >"nucleotide-differences.log"
    echo "md5sum_version ${meta.id}" > "md5sum_version.txt"

    # isolate warnings in plot (when grep find no match, exit code is 1, "|| true" fixes that)
    grep "WARN -" "nucleotide-differences.log" | sed 's/^\\[1\\] //g' | sed 's/"//g' > Warnings.txt || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | sed -n 1p | sed 's/R version //g' | sed 's/\\s.*\$//')
    END_VERSIONS
    """
}

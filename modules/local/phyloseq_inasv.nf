process PHYLOSEQ_INASV {
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path(biom_file)

    output:
    path( "*.tsv" )          , emit: tsv
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tail $biom_file -n +2 | sed '1s/#OTU ID/ASV_ID/' > reformat_$biom_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | sed -n 1p | grep -Eo 'version [[:alnum:].]+' | sed 's/version //')
    END_VERSIONS
    """
}

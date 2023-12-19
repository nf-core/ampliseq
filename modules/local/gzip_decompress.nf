process GZIP_DECOMPRESS {
    tag "$file"
    label 'process_single'

    conda "conda-forge::sed=4.7 conda-forge::gzip=1.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    path(file)

    output:
    path("$outfile"), emit: ungzip
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    outfile = task.ext.outfile ?: file.baseName.toString().replaceFirst(/\.gz$/, "")

    """
    gzip $args -c -d $file > $outfile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/gzip //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}

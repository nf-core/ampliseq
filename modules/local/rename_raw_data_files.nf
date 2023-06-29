process RENAME_RAW_DATA_FILES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}{_1,_2,}.fastq.gz", includeInputs: true), emit: fastq
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def args        = task.ext.args ?: 'ln -s'
    if (meta.single_end) {
        """
        if [ ! -f  ${meta.id}.fastq.gz ]; then
            $args $reads ${meta.id}.fastq.gz
        else
            touch ${meta.id}.fastq.gz
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
        END_VERSIONS
        """
    } else {
        """
        [ -f "${meta.id}_1.fastq.gz" ] || $args "${reads[0]}" "${meta.id}_1.fastq.gz"
        [ -f "${meta.id}_2.fastq.gz" ] || $args "${reads[1]}" "${meta.id}_2.fastq.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
        END_VERSIONS
        """
    }
}

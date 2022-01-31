process RENAME_RAW_DATA_FILES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}{_1,_2,}.fastq.gz", includeInputs: true), emit: fastq
    path "versions.yml"                                                      , emit: versions

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    if (meta.single_end) {
        """
        if [ ! -f  ${meta.id}.fastq.gz ]; then
            ln -s $reads ${meta.id}.fastq.gz
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
        [ -f "${meta.id}_1.fastq.gz" ] || ln -s "${reads[0]}" "${meta.id}_1.fastq.gz"
        [ -f "${meta.id}_2.fastq.gz" ] || ln -s "${reads[1]}" "${meta.id}_2.fastq.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            sed: \$(sed --version 2>&1 | sed -n 1p | sed 's/sed (GNU sed) //')
        END_VERSIONS
        """
    }
}

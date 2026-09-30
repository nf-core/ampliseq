process RENAME_RAW_DATA_FILES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ac/ac8d2a429a6d54f2642e14f2334a12493b42ae7769071d3da79030f4b5e3fd66/data' :
        'community.wave.seqera.io/library/bash_coreutils:25e9437236fbb54f' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}{_1,_2,}.fastq.gz", includeInputs: true), emit: fastq
    path "versions.yml"                                                      , emit: versions_rename_raw_data_files, topic: versions

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
            bash: \$(echo \"\$BASH_VERSION\")
        END_VERSIONS
        """
    } else {
        """
        [ -f "${meta.id}_1.fastq.gz" ] || $args "${reads[0]}" "${meta.id}_1.fastq.gz"
        [ -f "${meta.id}_2.fastq.gz" ] || $args "${reads[1]}" "${meta.id}_2.fastq.gz"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bash: \$(echo \"\$BASH_VERSION\")
        END_VERSIONS
        """
    }
}

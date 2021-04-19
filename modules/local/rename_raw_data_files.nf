// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process RENAME_RAW_DATA_FILES {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::cutadapt=3.2=py38h0213d0e_0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/cutadapt:3.2--py38h0213d0e_0'
    } else {
        container 'quay.io/biocontainers/cutadapt:3.2--py38h0213d0e_0'
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}{_1,_2,}.fastq.gz", includeInputs: true), emit: fastq

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    if (meta.single_end) {
        """
        [ ! -f  ${meta.id}.fastq.gz ] && ln -s $reads ${meta.id}.fastq.gz
        """
    } else {
        """
        [ -f "${meta.id}_1.fastq.gz" ] || ln -s "${reads[0]}" "${meta.id}_1.fastq.gz"
        [ -f "${meta.id}_2.fastq.gz" ] || ln -s "${reads[1]}" "${meta.id}_2.fastq.gz"
        """
    }
}
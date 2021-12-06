// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? 'bioconda::cutadapt=3.2' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/cutadapt:3.2--py38h0213d0e_0'
    } else {
        container 'quay.io/biocontainers/cutadapt:3.2--py38h0213d0e_0'
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trim.fast*')   , emit: reads
    tuple val(meta), path('*.log')          , emit: log
    path '*.version.txt'                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def fqsuff = ".fastq.gz"
    if (reads.toString().endsWith("fasta") || reads.toString().endsWith("fna")){
        fqsuff = ".fasta"
    }

    def trimmed  = meta.single_end ? "-o ${prefix}.trim${fqsuff}" : "-o ${prefix}_1.trim${fqsuff} -p ${prefix}_2.trim${fqsuff}"
    """
    cutadapt \\
        --cores $task.cpus \\
        $options.args \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log
    echo \$(cutadapt --version) > ${software}.version.txt
    """
}

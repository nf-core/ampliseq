// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process ITSX_CUTASV {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::itsx=1.1.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/itsx:1.1.3--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/itsx:1.1.3--hdfd78af_1"
    }

    input:
    path fasta
    
    output:
    path "ASV_ITS_seqs.full.fasta", emit: fasta
    path "*.version.txt"          , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    ITSx -i $fasta ${options}.args --cpu ${task.cpus} -o ASV_ITS_seqs

    ITSx -h 2>&1 > /dev/null | tail -n 2 | head -n 1 | cut -f 2 -d ' ' > ${software}.version.txt
    """
}

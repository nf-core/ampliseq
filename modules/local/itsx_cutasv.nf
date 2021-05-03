// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process ITSX_CUTASV {
    label 'process_low'
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
    path(fasta)
    
    output:
    path(  "ASV_ITS_seqs.full.fasta" ), emit: fasta
    //path( "*addSpecies.fna*"), emit: addspecies

    script:
    """
    ITSx -i $fasta -t all --preserve T --date F --positions F --graphical F --save_regions none --cpu ${task.cpus} -o ASV_ITS_seqs
    """
}

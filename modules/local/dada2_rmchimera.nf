// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_RMCHIMERA {
    tag "$meta.run"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconductor-dada2=1.18.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    }

    input:
    tuple val(meta), path(seqtab)

    output:
    tuple val(meta), path("*.ASVtable.rds"), emit: rds
    path "*.version.txt"                   , emit: version
    path "*.args.txt"                      , emit: args

    script:
    def software      = getSoftwareName(task.process)
    def no_samples    = meta.id.size()
    def first_sample  = meta.id.first()
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    seqtab = readRDS("${seqtab}")

    #remove chimera
    seqtab.nochim <- removeBimeraDenovo(seqtab, $options.args, multithread=$task.cpus, verbose=TRUE)
    if ( ${no_samples} == 1 ) { rownames(seqtab.nochim) <- "${first_sample}" }
    saveRDS(seqtab.nochim,"${meta.run}.ASVtable.rds")

    write.table('removeBimeraDenovo\t$options.args', file = "removeBimeraDenovo.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}

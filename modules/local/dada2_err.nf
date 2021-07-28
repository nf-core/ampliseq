// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_ERR {
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
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.err.rds"), emit: errormodel
    tuple val(meta), path("*.err.pdf"), emit: pdf
    tuple val(meta), path("*.err.log"), emit: log
    tuple val(meta), path("*.err.convergence.txt"), emit: convergence
    path "*.version.txt"              , emit: version
    path "*.args.txt"                 , emit: args

    script:
    def software      = getSoftwareName(task.process)
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        fnFs <- sort(list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE))
        fnRs <- sort(list.files(".", pattern = "_2.filt.fastq.gz", full.names = TRUE))

        sink(file = "${meta.run}.err.log")
        errF <- learnErrors(fnFs, $options.args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${meta.run}_1.err.rds")
        errR <- learnErrors(fnRs, $options.args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errR, "${meta.run}_2.err.rds")
        sink(file = NULL)

        pdf("${meta.run}_1.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        pdf("${meta.run}_2.err.pdf")
        plotErrors(errR, nominalQ = TRUE)
        dev.off()

        sink(file = "${meta.run}_1.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        sink(file = "${meta.run}_2.err.convergence.txt")
        dada2:::checkConvergence(errR)
        sink(file = NULL)

        write.table('learnErrors\t$options.args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        fnFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        sink(file = "${meta.run}.err.log")
        errF <- learnErrors(fnFs, $options.args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${meta.run}.err.rds")
        sink(file = NULL)

        pdf("${meta.run}.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        sink(file = "${meta.run}.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        write.table('learnErrors\t$options.args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
    }
}

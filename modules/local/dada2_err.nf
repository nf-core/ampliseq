process DADA2_ERR {
    tag "$meta.run"
    label 'process_medium'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.err.rds"), emit: errormodel
    tuple val(meta), path("*.err.pdf"), emit: pdf
    tuple val(meta), path("*.err.log"), emit: log
    tuple val(meta), path("*.err.convergence.txt"), emit: convergence
    path "versions.yml"               , emit: versions
    path "*.args.txt"                 , emit: args

    script:
    def args = task.ext.args ?: ''
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        fnFs <- sort(list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE))
        fnRs <- sort(list.files(".", pattern = "_2.filt.fastq.gz", full.names = TRUE))

        sink(file = "${meta.run}.err.log")
        errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${meta.run}_1.err.rds")
        errR <- learnErrors(fnRs, $args, multithread = $task.cpus, verbose = TRUE)
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

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        fnFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        sink(file = "${meta.run}.err.log")
        errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${meta.run}.err.rds")
        sink(file = NULL)

        pdf("${meta.run}.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        sink(file = "${meta.run}.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
}

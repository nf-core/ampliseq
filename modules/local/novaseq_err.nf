process NOVASEQ_ERR {
    tag "$meta.run"
    label 'process_medium'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"
    input:
    tuple val(meta), path(errormodel)

    output:
    tuple val(meta), path("*.md.err.rds"), emit: errormodel
    tuple val(meta), path("*.md.err.pdf"), emit: pdf
    tuple val(meta), path("*.md.err.log"), emit: log
    tuple val(meta), path("*.md.err.convergence.txt"), emit: convergence

    script:
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        errF = readRDS("${errormodel[0]}")
        errR = readRDS("${errormodel[1]}")

        #monotone decreasing
        sink(file = "${meta.run}.md.err.log")
        make.monotone.decreasing <- function(v) sapply(seq_along(v), function(i) max(v[i:length(v)]))

        errF.md <- t(apply(getErrors(errF), 1, make.monotone.decreasing))
        errF.md.full <- errF
        errF.md.full$err_out <- errF.md
        saveRDS(errF.md.full, "${meta.run}_1.md.err.rds")
        
        errR.md <- t(apply(getErrors(errR), 1, make.monotone.decreasing))
        errR.md.full <- errR
        errR.md.full$err_out <- errR.md
        saveRDS(errR.md.full, "${meta.run}_2.md.err.rds")
        sink(file = NULL)

        pdf("${meta.run}_1.md.err.pdf")
        plotErrors(errF.md.full, nominalQ = TRUE)
        dev.off()

        pdf("${meta.run}_2.md.err.pdf")
        plotErrors(errR.md.full, nominalQ = TRUE)
        dev.off()

        sink(file = "${meta.run}_1.md.err.convergence.txt")
        dada2:::checkConvergence(errF.md.full)
        sink(file = NULL)

        sink(file = "${meta.run}_2.md.err.convergence.txt")
        dada2:::checkConvergence(errR.md.full)
        sink(file = NULL)
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        errF = readRDS("${errormodel}")

        #monotone decreasing
        sink(file = "${meta.run}.md.err.log")
        make.monotone.decreasing <- function(v) sapply(seq_along(v), function(i) max(v[i:length(v)]))

        errF.md <- t(apply(getErrors(errF), 1, make.monotone.decreasing))
        errF.md.full <- errF
        errF.md.full$err_out <- errF.md
        saveRDS(errF.md.full, "${meta.run}.md.err.rds")
        sink(file = NULL)

        pdf("${meta.run}.md.err.pdf")
        plotErrors(errF.md.full, nominalQ = TRUE)
        dev.off()

        sink(file = "${meta.run}.md.err.convergence.txt")
        dada2:::checkConvergence(errF.md.full)
        sink(file = NULL)
        """
    }
}
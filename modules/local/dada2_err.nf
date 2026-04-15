process DADA2_ERR {
    tag "$meta.run"
    label 'process_medium'

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data' :
        'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.err.rds"), emit: errormodel
    tuple val(meta), path("*.err.pdf"), emit: pdf
    tuple val(meta), path("*.err.svg"), emit: svg
    tuple val(meta), path("*.err.log"), emit: log
    tuple val(meta), path("*.err.convergence.txt"), emit: convergence
    path "versions.yml"               , emit: versions
    path "*.args.txt"                 , emit: args


    script:
    def cmd_errfun = task.ext.cmd_errfun ?: ""
    def prefix = task.ext.prefix ?: "prefix"
    def args = task.ext.args ?: ''
    def seed = task.ext.seed ?: '100'
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))
        set.seed($seed) # Initialize random number generator for reproducibility

        ${cmd_errfun}

        fnFs <- list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE)
        fnRs <- sub("_1.filt.fastq.gz", "_2.filt.fastq.gz", fnFs)

        sink(file = "${prefix}.err.log")
        errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${prefix}_1.err.rds")
        errR <- learnErrors(fnRs, $args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errR, "${prefix}_2.err.rds")
        sink(file = NULL)

        pdf("${prefix}_1.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()
        svg("${prefix}_1.err.svg")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        pdf("${prefix}_2.err.pdf")
        plotErrors(errR, nominalQ = TRUE)
        dev.off()
        svg("${prefix}_2.err.svg")
        plotErrors(errR, nominalQ = TRUE)
        dev.off()

        sink(file = "${prefix}_1.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        sink(file = "${prefix}_2.err.convergence.txt")
        dada2:::checkConvergence(errR)
        sink(file = NULL)

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))
        set.seed($seed) # Initialize random number generator for reproducibility

        ${cmd_errfun}

        fnFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        sink(file = "${prefix}.err.log")
        errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${prefix}.err.rds")
        sink(file = NULL)

        pdf("${prefix}.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()
        svg("${prefix}.err.svg")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        sink(file = "${prefix}.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
}

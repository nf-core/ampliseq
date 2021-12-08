// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_DENOISING {
    tag "$meta.run"
    label 'process_medium'
    label 'process_long'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconductor-dada2=1.20.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.20.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.20.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path("filtered/*"), path(errormodel)

    output:
    tuple val(meta), path("*.dada.rds")   , emit: denoised
    tuple val(meta), path("*.seqtab.rds") , emit: seqtab
    tuple val(meta), path("*.mergers.rds"), emit: mergers
    tuple val(meta), path("*.log")        , emit: log
    path "*.version.txt"                  , emit: version
    path "*.args.txt"                     , emit: args

    script:
    def software      = getSoftwareName(task.process)
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        errF = readRDS("${errormodel[0]}")
        errR = readRDS("${errormodel[1]}")

        filtFs <- sort(list.files("./filtered/", pattern = "_1.filt.fastq.gz", full.names = TRUE))
        filtRs <- sort(list.files("./filtered/", pattern = "_2.filt.fastq.gz", full.names = TRUE))

        #denoising
        sink(file = "${meta.run}.dada.log")
        dadaFs <- dada(filtFs, err = errF, $options.args, multithread = $task.cpus)
        saveRDS(dadaFs, "${meta.run}_1.dada.rds")
        dadaRs <- dada(filtRs, err = errR, $options.args, multithread = $task.cpus)
        saveRDS(dadaRs, "${meta.run}_2.dada.rds")
        sink(file = NULL)

        #make table
        mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, $options.args2, verbose=TRUE)
        saveRDS(mergers, "${meta.run}.mergers.rds")
        seqtab <- makeSequenceTable(mergers)
        saveRDS(seqtab, "${meta.run}.seqtab.rds")

        write.table('dada\t$options.args', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        write.table('mergePairs\t$options.args2', file = "mergePairs.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        errF = readRDS("${errormodel}")

        filtFs <- sort(list.files("./filtered/", pattern = ".fastq.gz", full.names = TRUE))

        #denoising
        sink(file = "${meta.run}.dada.log")
        dadaFs <- dada(filtFs, err = errF, $options.args, multithread = $task.cpus)
        saveRDS(dadaFs, "${meta.run}.dada.rds")
        sink(file = NULL)

        #make table
        seqtab <- makeSequenceTable(dadaFs)
        saveRDS(seqtab, "${meta.run}.seqtab.rds")

        #dummy file to fulfill output rules
        saveRDS("dummy", "dummy_${meta.run}.mergers.rds")

        write.table('dada\t$options.args', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        """
    }
}

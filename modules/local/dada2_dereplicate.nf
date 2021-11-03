// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_DEREPLICATE {
    tag "$meta.run"
    label 'process_medium'

    conda (params.enable_conda ? "bioconductor-dada2=1.20.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.20.0--r41h399db7b_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.20.0--r41h399db7b_0"
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.derep.rds"), emit: dereplicated
    path "*.version.txt"                , emit: version

    script:
    def software      = getSoftwareName(task.process)
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        filtFs <- sort(list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE))
        filtRs <- sort(list.files(".", pattern = "_2.filt.fastq.gz", full.names = TRUE))

        derepFs <- derepFastq(filtFs, n = 10000, verbose = TRUE)
        saveRDS(derepFs, "${meta.run}_1.derep.rds")
        derepRs <- derepFastq(filtRs, n = 10000, verbose = TRUE)
        saveRDS(derepRs, "${meta.run}_2.derep.rds")

        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        filtFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        derepFs <- derepFastq(filtFs, n = 10000, verbose = TRUE)
        saveRDS(derepFs, "${meta.run}.derep.rds")

        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
    }
}

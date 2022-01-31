process DADA2_FILTNTRIM {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    tuple val(meta), path(reads), val(trunclenf), val(trunclenr)

    output:
    tuple val(meta), path("*.filter_stats.tsv"), emit: log
    tuple val(meta), path("*.filt.fastq.gz")   , emit: reads
    path "versions.yml"                        , emit: versions
    path "*.args.txt"                          , emit: args

    script:
    def args        = task.ext.args ?: ''
    def in_and_out  = meta.single_end ? "\"${reads}\", \"${meta.id}.filt.fastq.gz\"" : "\"${reads[0]}\", \"${meta.id}_1.filt.fastq.gz\", \"${reads[1]}\", \"${meta.id}_2.filt.fastq.gz\""
    def trunclenf   = trunclenf[1].toInteger()
    def trunclenr   = trunclenr[1].toInteger()
    def trunc_args  = meta.single_end ? "truncLen = $trunclenf" : "truncLen = c($trunclenf, $trunclenr)"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    out <- filterAndTrim($in_and_out,
        $trunc_args,
        $args,
        compress = TRUE,
        multithread = $task.cpus,
        verbose = TRUE)
    out <- cbind(out, ID = row.names(out))

    write.table( out, file = "${meta.id}.filter_stats.tsv", sep = "\\t", row.names = FALSE, quote = FALSE, na = '')
    write.table(paste('filterAndTrim\t$trunc_args','$args',sep=","), file = "filterAndTrim.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

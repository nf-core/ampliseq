process DADA2_FILTNTRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(reads), val(trunclenf), val(trunclenr)

    output:
    tuple val(meta), path("*.filt.fastq.gz"), path("*.filter_stats.tsv"), path("*.args.txt"), emit: reads_logs_args
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def in_and_out  = meta.single_end ? "\"${reads}\", \"${meta.id}.filt.fastq.gz\"" : "\"${reads[0]}\", \"${meta.id}_1.filt.fastq.gz\", \"${reads[1]}\", \"${meta.id}_2.filt.fastq.gz\""
    def outfiles    = meta.single_end ? "\"${meta.id}.filt.fastq.gz\"" : "\"${meta.id}_1.filt.fastq.gz\", \"${meta.id}_2.filt.fastq.gz\""
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

    # If no reads passed the filter, write an empty GZ file
    if(out[2] == '0'){
        for(fp in c($outfiles)){
            print(paste("Writing out an empty file:", fp))
            handle <- gzfile(fp, "w")
            write("", handle)
            close(handle)
        }
    }

    write.table( out, file = "${meta.id}.filter_stats.tsv", sep = "\\t", row.names = FALSE, quote = FALSE, na = '')
    write.table(paste('filterAndTrim\t$trunc_args','$args',sep=","), file = "filterAndTrim.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

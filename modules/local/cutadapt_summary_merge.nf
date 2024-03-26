process CUTADAPT_SUMMARY_MERGE {
    tag "${files}"
    label 'process_low'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    val(action)
    path(files)

    output:
    path("cutadapt_summary.tsv")      , emit: tsv
    path "versions.yml", optional:true, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if (action == "merge") {
        """
        #!/usr/bin/env Rscript
        standard <- read.table(\"${files[0]}\", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
        doubleprimer <- read.table(\"${files[1]}\", header = TRUE, sep = "\\t", stringsAsFactors = FALSE)
        colnames(doubleprimer) <- c("sample", "cutadapt_doubleprimer_total_processed", "cutadapt_doubleprimer_reverse_complemented", "cutadapt_doubleprimer_passing_filters", "cutadapt_doubleprimer_passing_filters_percent")

        #merge
        df <- merge(standard, doubleprimer, by = "sample")

        #filter columns
        remove_columns <- c("cutadapt_doubleprimer_total_processed")
        for(column in remove_columns) df[column]<-NULL

        #write
        write.table(df, file = \"cutadapt_summary.tsv\", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\\t")

        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")) ), "versions.yml")
        """
    } else {
        """
        cp $files cutadapt_summary.tsv
        """
    }
}

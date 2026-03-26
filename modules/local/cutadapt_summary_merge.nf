process CUTADAPT_SUMMARY_MERGE {
    tag "${files}"
    label 'process_single'

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data' :
        'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    val(action)
    path(files)

    output:
    path("cutadapt_summary.tsv")      , emit: tsv
    path "versions.yml", optional:true, emit: versions


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

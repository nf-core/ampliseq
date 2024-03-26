process DADA2_STATS {
    tag "$meta.run"
    label 'process_low'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path("filter_and_trim_files/*"), path(denoised), path(mergers), path(seqtab_nochim)

    output:
    tuple val(meta), path("*.stats.tsv"), emit: stats
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "prefix"
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        #combine filter_and_trim files
        for (data in list.files("./filter_and_trim_files", full.names=TRUE)){
            if (!exists("filter_and_trim")){ filter_and_trim <- read.csv(data, header=TRUE, sep="\\t") }
            if (exists("filter_and_trim")){
                tempory <-read.csv(data, header=TRUE, sep="\\t")
                filter_and_trim <-unique(rbind(filter_and_trim, tempory))
                rm(tempory)
            }
        }
        rownames(filter_and_trim) <- filter_and_trim\$ID
        filter_and_trim["ID"] <- NULL
        #write.table( filter_and_trim, file = "${prefix}.filter_and_trim.tsv", sep = "\\t", row.names = TRUE, quote = FALSE, na = '')

        #read data
        dadaFs = readRDS("${denoised[0]}")
        dadaRs = readRDS("${denoised[1]}")
        mergers = readRDS("$mergers")
        seqtab.nochim = readRDS("$seqtab_nochim")

        #track reads through pipeline
        getN <- function(x) sum(getUniques(x))
        if ( nrow(filter_and_trim) == 1 ) {
            track <- cbind(filter_and_trim, getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
        } else {
            track <- cbind(filter_and_trim, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
        }
        colnames(track) <- c("DADA2_input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
        rownames(track) <- sub(pattern = "_1.fastq.gz\$", replacement = "", rownames(track)) #this is when cutadapt is skipped!
        track <- cbind(sample = sub(pattern = "(.*?)\\\\..*\$", replacement = "\\\\1", rownames(track)), track)
        write.table( track, file = "${prefix}.stats.tsv", sep = "\\t", row.names = FALSE, quote = FALSE, na = '')

        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        #combine filter_and_trim files
        for (data in list.files("./filter_and_trim_files", full.names=TRUE)){
            if (!exists("filter_and_trim")){ filter_and_trim <- read.csv(data, header=TRUE, sep="\\t") }
            if (exists("filter_and_trim")){
                tempory <-read.csv(data, header=TRUE, sep="\\t")
                filter_and_trim <-unique(rbind(filter_and_trim, tempory))
                rm(tempory)
            }
        }
        rownames(filter_and_trim) <- filter_and_trim\$ID
        filter_and_trim["ID"] <- NULL
        #write.table( filter_and_trim, file = "${prefix}.filter_and_trim.tsv", sep = "\\t", row.names = TRUE, quote = FALSE, na = '')

        #read data
        dadaFs = readRDS("${denoised[0]}")
        seqtab.nochim = readRDS("$seqtab_nochim")

        #track reads through pipeline
        getN <- function(x) sum(getUniques(x))
        if ( nrow(filter_and_trim) == 1 ) {
            track <- cbind(filter_and_trim, getN(dadaFs), rowSums(seqtab.nochim))
        } else {
            track <- cbind(filter_and_trim, sapply(dadaFs, getN), rowSums(seqtab.nochim))
        }
        colnames(track) <- c("DADA2_input", "filtered", "denoised", "nonchim")
        track <- cbind(sample = sub(pattern = "(.*?)\\\\..*\$", replacement = "\\\\1", rownames(track)), track)
        write.table( track, file = "${prefix}.stats.tsv", sep = "\\t", row.names = FALSE, quote = FALSE, na = '')

        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
}

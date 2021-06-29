// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_STATS {
    tag "$meta.run"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconductor-dada2=1.18.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    }

    input:
    tuple val(meta), path("filter_and_trim_files/*"), path(denoised), path(mergers), path(seqtab_nochim)
    
    output:
    tuple val(meta), path("*.stats.tsv"), emit: stats
    path "*.version.txt"                , emit: version

    script:
    def software      = getSoftwareName(task.process)
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        #combine filter_and_trim files
        for (data in list.files("./filter_and_trim_files", full.names=TRUE)){
            if (!exists("filter_and_trim")){ filter_and_trim <- read.csv(data, header=TRUE, sep="\t") }
            if (exists("filter_and_trim")){
                tempory <-read.csv(data, header=TRUE, sep="\t")
                filter_and_trim <-unique(rbind(filter_and_trim, tempory))
                rm(tempory)
            }
        }
        rownames(filter_and_trim) <- filter_and_trim\$ID
        filter_and_trim["ID"] <- NULL
        #write.table( filter_and_trim, file = "${meta.run}.filter_and_trim.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

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
        track <- cbind(sample = sub(pattern = "(.*?)\\\\..*\$", replacement = "\\\\1", rownames(track)), track)
        track\$sample <- sub(pattern = "_1\$", replacement = "", track\$sample)
        write.table( track, file = "${meta.run}.stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        #combine filter_and_trim files
        for (data in list.files("./filter_and_trim_files", full.names=TRUE)){
            if (!exists("filter_and_trim")){ filter_and_trim <- read.csv(data, header=TRUE, sep="\t") }
            if (exists("filter_and_trim")){
                tempory <-read.csv(data, header=TRUE, sep="\t")
                filter_and_trim <-unique(rbind(filter_and_trim, tempory))
                rm(tempory)
            }
        }
        rownames(filter_and_trim) <- filter_and_trim\$ID
        filter_and_trim["ID"] <- NULL
        #write.table( filter_and_trim, file = "${meta.run}.filter_and_trim.tsv", sep = "\t", row.names = TRUE, quote = FALSE)

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
        colnames(track) <- c("input", "filtered", "denoised", "nonchim")
        track <- cbind(sample = sub(pattern = "(.*?)\\\\..*\$", replacement = "\\\\1", rownames(track)), track)
        write.table( track, file = "${meta.run}.stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """        
    }
}
// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_QUALITY {
    tag "$meta"
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
    tuple val(meta), path(reads)
    
    output:
    path "${meta}_qual_stats.pdf"            , emit: pdf
    tuple val(meta), path("*_qual_stats.tsv"), emit: tsv
    path "*.args.txt"                        , emit: args

    script:
    """
    dada_quality.r "${meta}_qual_stats" $options.args
    echo 'plotQualityProfile\t$options.args' > "plotQualityProfile.args.txt"
    """
}

process DADA2_FILTNTRIM {
    tag "$meta.id"
    label 'process_medium'
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
    tuple val(meta), path(reads), val(trunclenf), val(trunclenr)
    
    output:
    tuple val(meta), path("*.filter_stats.tsv"), emit: log
    tuple val(meta), path("*.filt.fastq.gz")   , emit: reads
    path "*.version.txt"                       , emit: version
    path "*.args.txt"                          , emit: args

    script:
    def software    = getSoftwareName(task.process)
    def in_and_out  = meta.single_end ? "\"${reads}\", \"${meta.id}.filt.fastq.gz\"" : "\"${reads[0]}\", \"${meta.id}_1.filt.fastq.gz\", \"${reads[1]}\", \"${meta.id}_2.filt.fastq.gz\""
    def trunclenf   = trunclenf[1].toInteger()
    def trunclenr   = trunclenr[1].toInteger()
    def trunc_args  = meta.single_end ? "truncLen = $trunclenf" : "truncLen = c($trunclenf, $trunclenr)"
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    out <- filterAndTrim($in_and_out, 
        $trunc_args,
        $options.args,
        compress = TRUE, 
        multithread = $task.cpus, 
        verbose = TRUE)
    out <- cbind(out, ID = row.names(out))

    write.table( out, file = "${meta.id}.filter_stats.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(paste('filterAndTrim\t$trunc_args','$options.args',sep=","), file = "filterAndTrim.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}

process DADA2_ERR {
    tag "$meta.run"
    label 'process_medium'
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
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("*.err.rds"), emit: errormodel
    tuple val(meta), path("*.err.pdf"), emit: pdf
    tuple val(meta), path("*.err.log"), emit: log
    tuple val(meta), path("*.err.convergence.txt"), emit: convergence
    path "*.version.txt"              , emit: version
    path "*.args.txt"                 , emit: args

    script:
    def software      = getSoftwareName(task.process)
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        fnFs <- sort(list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE))
        fnRs <- sort(list.files(".", pattern = "_2.filt.fastq.gz", full.names = TRUE))

        sink(file = "${meta.run}.err.log")
        errF <- learnErrors(fnFs, $options.args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${meta.run}_1.err.rds")
        errR <- learnErrors(fnRs, $options.args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errR, "${meta.run}_2.err.rds")
        sink(file = NULL)

        pdf("${meta.run}_1.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        pdf("${meta.run}_2.err.pdf")
        plotErrors(errR, nominalQ = TRUE)
        dev.off()
        
        sink(file = "${meta.run}_1.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        sink(file = "${meta.run}_2.err.convergence.txt")
        dada2:::checkConvergence(errR)
        sink(file = NULL)

        write.table('learnErrors\t$options.args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        fnFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        sink(file = "${meta.run}.err.log")
        errF <- learnErrors(fnFs, $options.args, multithread = $task.cpus, verbose = TRUE)
        saveRDS(errF, "${meta.run}.err.rds")
        sink(file = NULL)

        pdf("${meta.run}.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        sink(file = "${meta.run}.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        write.table('learnErrors\t$options.args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """        
    }
}

process DADA2_DEREPLICATE {
    tag "$meta.run"
    label 'process_medium'

    conda (params.enable_conda ? "bioconductor-dada2=1.18.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h5f743cb_0"
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

        derepFs <- derepFastq(filtFs, verbose = TRUE)
        saveRDS(derepFs, "${meta.run}_1.derep.rds")
        derepRs <- derepFastq(filtRs, verbose = TRUE)
        saveRDS(derepRs, "${meta.run}_2.derep.rds")

        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        filtFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        derepFs <- derepFastq(filtFs, verbose = TRUE)
        saveRDS(derepFs, "${meta.run}.derep.rds")

        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """        
    }
}

process DADA2_DENOISING {
    tag "$meta.run"
    label 'process_medium'
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
    tuple val(meta), path(dereplicated), path(errormodel)
    
    output:
    tuple val(meta), path("*.dada.rds")   , emit: denoised
    tuple val(meta), path("*.seqtab.rds") , emit: seqtab
    tuple val(meta), path("*.mergers.rds"), emit: mergers
    tuple val(meta), path("*.dada.log")   , emit: log
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

        derepFs = readRDS("${dereplicated[0]}")
        derepRs = readRDS("${dereplicated[1]}")

        #denoising
        sink(file = "${meta.run}.dada.log")
        dadaFs <- dada(derepFs, err = errF, $options.args, multithread = $task.cpus)
        saveRDS(dadaFs, "${meta.run}_1.dada.rds")
        dadaRs <- dada(derepRs, err = errR, $options.args, multithread = $task.cpus)
        saveRDS(dadaRs, "${meta.run}_2.dada.rds")
        sink(file = NULL)

        #make table
        mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, $options.args2, verbose=TRUE)
        saveRDS(mergers, "${meta.run}.mergers.rds")
        seqtab <- makeSequenceTable(mergers)
        saveRDS(seqtab, "${meta.run}.seqtab.rds")

        write.table('dada\t$options.args', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table('mergePairs\t$options.args2', file = "mergePairs.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        errF = readRDS("${errormodel}")

        derepFs = readRDS("${dereplicated}")

        #denoising
        sink(file = "${meta.run}.dada.log")
        dadaFs <- dada(derepFs, err = errF, $options.args, multithread = $task.cpus)
        saveRDS(dadaFs, "${meta.run}.dada.rds")
        sink(file = NULL)

        #make table
        seqtab <- makeSequenceTable(dadaFs)
        saveRDS(seqtab, "${meta.run}.seqtab.rds")

        #dummy file to fulfill output rules
        saveRDS("dummy", "dummy_${meta.run}.mergers.rds")

        write.table('dada\t$options.args', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
        """        
    }
}

process DADA2_RMCHIMERA {
    tag "$meta.run"
    label 'process_medium'
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
    tuple val(meta), path(seqtab)
    
    output:
    tuple val(meta), path("*.ASVtable.rds"), emit: rds
    path "*.version.txt"                   , emit: version
    path "*.args.txt"                      , emit: args

    script:
    def software      = getSoftwareName(task.process)
    def no_samples    = meta.id.size()
    def first_sample  = meta.id.first()
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    seqtab = readRDS("${seqtab}")

    #remove chimera
    seqtab.nochim <- removeBimeraDenovo(seqtab, $options.args, multithread=$task.cpus, verbose=TRUE)
    if ( ${no_samples} == 1 ) { rownames(seqtab.nochim) <- "${first_sample}" }
    saveRDS(seqtab.nochim,"${meta.run}.ASVtable.rds")

    write.table('removeBimeraDenovo\t$options.args', file = "removeBimeraDenovo.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}

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
        colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
        track <- cbind(sample = sub(pattern = "(.*?)\\\\..*\$", replacement = "\\\\1", rownames(track)), track)
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

process DADA2_MERGE {
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
    path(files)
    path(rds)
    
    output:
    path( "DADA2_stats.tsv" ), emit: dada2stats
    path( "DADA2_table.tsv" ), emit: dada2asv
    path( "ASV_table.tsv" )  , emit: asv
    path( "ASV_seqs.fasta" ) , emit: fasta
    path( "DADA2_table.rds" ), emit: rds
    path "*.version.txt"     , emit: version

    script:
    def software      = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))

    #combine stats files
    for (data in sort(list.files(".", pattern = ".stats.tsv", full.names = TRUE))) {
        if (!exists("stats")){ stats <- read.csv(data, header=TRUE, sep="\t") }
        if (exists("stats")){
            temp <-read.csv(data, header=TRUE, sep="\t")
            stats <-unique(rbind(stats, temp))
            rm(temp)
        }
    }
    write.table( stats, file = "DADA2_stats.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    #combine dada-class objects
    files <- sort(list.files(".", pattern = ".ASVtable.rds", full.names = TRUE))
    if ( length(files) == 1 ) {
        ASVtab = readRDS(files[1])
    } else {
        ASVtab <- mergeSequenceTables(tables=files, repeats = "error", orderBy = "abundance", tryRC = FALSE)
    }
    saveRDS(ASVtab, "DADA2_table.rds")

    df <- t(ASVtab)
    colnames(df) <- gsub('_1.filt.fastq.gz', '', colnames(df))
    colnames(df) <- gsub('.filt.fastq.gz', '', colnames(df))
    df <- data.frame(sequence = rownames(df), df)
    row.names(df) <- paste0("ASV_", seq(nrow(df)))
    df <- data.frame(ASV_ID=row.names(df), df)

    # file to publish
    write.table( df, file = "DADA2_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

    # Write fasta file with ASV sequences to file
    fasta.tab <- df[,c("ASV_ID","sequence")]
    fasta.tab\$ASV_ID <- gsub("ASV_",">ASV_",fasta.tab\$ASV_ID)
    fasta.tab.join <- c(rbind( fasta.tab\$ASV_ID, fasta.tab\$sequence ))
    write( fasta.tab.join, file = 'ASV_seqs.fasta' )

    # Write ASV file with ASV abundances to file
    df\$sequence <- NULL
    df\$ASV_ID <- NULL
    write.table( df, file = "ASV_table.tsv", sep="\t", row.names = TRUE, col.names = NA, quote = FALSE)

    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}

process DADA2_TAXONOMY {
    label 'process_high'
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
    path(fasta)
    path(database)
    
    output:
    path( "ASV_tax.tsv" ), emit: tsv
    path( "ASV_tax.rds" ), emit: rds
    path "*.version.txt" , emit: version
    path "*.args.txt"    , emit: args

    script:
    def software      = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    seq <- getSequences(\"$fasta\", collapse = TRUE, silence = FALSE)
    taxa <- assignTaxonomy(seq, \"$database\", $options.args, multithread = $task.cpus, verbose=TRUE)
    saveRDS(taxa, "ASV_tax.rds")

    taxa <- cbind(sequence = rownames(taxa), taxa)
    write.table(taxa, file = "ASV_tax.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    write.table('assignTaxonomy\t$options.args', file = "assignTaxonomy.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}

process DADA2_ADDSPECIES {
    label 'process_high'
    label 'single_cpu'
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
    path(taxtable)
    path(database)
    
    output:
    path( "ASV_tax_species.tsv" ), emit: tsv
    path "*.version.txt"         , emit: version
    path "*.args.txt"            , emit: args

    script:
    def software      = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    taxtable <- readRDS(\"$taxtable\")
    taxa <- addSpecies(taxtable, \"$database\", $options.args, verbose=TRUE)
    taxa <- cbind(sequence = rownames(taxa), taxa)
    write.table(taxa, file = "ASV_tax_species.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    write.table('addSpecies\t$options.args', file = "addSpecies.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}
process DADA2_DENOISING {
    tag "$meta.run"
    label 'process_medium'
    label 'process_long'

    conda "bioconda::bioconductor-dada2=1.30.0 conda-forge::r-base=4.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path("filtered/*"), path(errormodel)

    output:
    tuple val(meta), path("*.dada.rds")   , emit: denoised
    tuple val(meta), path("*.seqtab.rds") , emit: seqtab
    tuple val(meta), path("*.mergers.rds"), emit: mergers
    tuple val(meta), path("*.log")        , emit: log
    path "versions.yml"                   , emit: versions
    path "*.args.txt"                     , emit: args


    script:
    def prefix = task.ext.prefix ?: "prefix"
    def quality_type = task.ext.quality_type ?: "Auto"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def quantile = task.ext.quantile ?: 0.001
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        errF <- readRDS("${errormodel[0]}")
        errR <- readRDS("${errormodel[1]}")

        filtFs <- sort(list.files("./filtered/", pattern = "_1.filt.fastq.gz", full.names = TRUE), method = "radix")
        filtRs <- sort(list.files("./filtered/", pattern = "_2.filt.fastq.gz", full.names = TRUE), method = "radix")

        #denoising
        sink(file = "${prefix}.dada.log")
        if ("${quality_type}" == "Auto") {
            # Avoid using memory-inefficient derepFastq() if not necessary
            dadaFs <- dada(filtFs, err = errF, $args, multithread = $task.cpus)
            dadaRs <- dada(filtRs, err = errR, $args, multithread = $task.cpus)
        } else {
            derepFs <- derepFastq(filtFs, qualityType="${quality_type}")
            dadaFs <- dada(derepFs, err = errF, $args, multithread = $task.cpus)
            derepRs <- derepFastq(filtRs, qualityType="${quality_type}")
            dadaRs <- dada(derepRs, err = errR, $args, multithread = $task.cpus)
        }
        saveRDS(dadaFs, "${prefix}_1.dada.rds")
        saveRDS(dadaRs, "${prefix}_2.dada.rds")
        sink(file = NULL)

        # merge
        if ("${params.mergepairs_strategy}" == "consensus") {
            mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, $args2, justConcatenate = FALSE, verbose=TRUE)
            concats <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, $args2, justConcatenate = TRUE, verbose=TRUE)

            # in case there is only one sample in the entire run
            if (is.data.frame(mergers)) {
                mergers <- list(sample = mergers)
                concats <- list(sample = concats)
            }

            # define the overlap threshold to decide if concatenation or not
            min_overlap_obs <- lapply(mergers, function(X) {
                mergers_accepted <- X[["accept"]]
                if (sum(mergers_accepted) > 0) {
                    min_overlap_obs <- X[["nmatch"]][mergers_accepted] + X[["nmismatch"]][mergers_accepted]
                    rep(min_overlap_obs, X[["abundance"]][mergers_accepted])
                } else {
                    NA
                }
            })

            min_overlap_obs <- Reduce(c, min_overlap_obs)
            min_overlap_obs <- min_overlap_obs[!is.na(min_overlap_obs)]
            min_overlap_obs <- quantile(min_overlap_obs, $quantile)

            for (x in names(mergers)) {
                to_concat <- !mergers[[x]][["accept"]] & (mergers[[x]][["nmismatch"]] + mergers[[x]][["nmatch"]]) < min_overlap_obs

                if (sum(to_concat) > 0) {
                    mergers[[x]][to_concat, ] <- concats[[x]][to_concat, ]
                    # filter out unaccepted non concatenated sequences
                    mergers[[x]] <- mergers[[x]][mergers[[x]][["accept"]], ]
                }

            }

            # if one sample, need to convert back to df for next steps

            if(length(mergers) == 1) {
                mergers <- mergers[[1]]
            }

        } else {
            mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, $args2, verbose=TRUE)
        }

        saveRDS(mergers, "${prefix}.mergers.rds")

        # make table
        seqtab <- makeSequenceTable(mergers)
        saveRDS(seqtab, "${prefix}.seqtab.rds")

        write.table('dada\t$args', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        write.table('mergePairs\t$args2', file = "mergePairs.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))

        errF = readRDS("${errormodel}")

        filtFs <- sort(list.files("./filtered/", pattern = ".fastq.gz", full.names = TRUE))

        #denoising
        sink(file = "${prefix}.dada.log")
        if ("${quality_type}" == "Auto") {
            # Avoid using memory-inefficient derepFastq() if not necessary
            dadaFs <- dada(filtFs, err = errF, $args, multithread = $task.cpus)
        } else {
            derepFs <- derepFastq(filtFs, qualityType="${quality_type}")
            dadaFs <- dada(derepFs, err = errF, $args, multithread = $task.cpus)
        }
        saveRDS(dadaFs, "${prefix}.dada.rds")
        sink(file = NULL)

        #make table
        seqtab <- makeSequenceTable(dadaFs)
        saveRDS(seqtab, "${prefix}.seqtab.rds")

        #dummy file to fulfill output rules
        saveRDS("dummy", "dummy_${prefix}.mergers.rds")

        write.table('dada\t$args', file = "dada.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
}

process DECONTAM {
    label 'process_single'

    conda "bioconda::bioconductor-decontam=1.30.0 conda-forge::r-base=4.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/375f1cd63fd980cd23dd3ed56e77fbb40d0bca26bfb8137eef60d473123f2426/data' :
        'community.wave.seqera.io/library/bioconductor-decontam_r-base:b0f82fb5ac5f1dc2' }"

    input:
    path(abundances_tsv)                   // first column with feature ids, other columns with numbers; first row with header
    path(metadata_tsv)                     // with columns 'control' & 'quant_reading' (e.g. PicoGreen fluroescent intensity) from samplesheet
    val(iscontaminant_method)              // one of "auto", "frequency", "prevalence", "combined", "minimum", "either", "both"
    val(iscontaminant_threshold)           // default: 0.1
    val(isnotcontaminant_threshold)        // default: 0.5

    output:
    path("decontaminated.tsv")        , emit: decontaminated_abundances
    path("decontaminated_counts.tsv") , emit: decontaminated_counts
    path("decontaminated_log.tsv")    , emit: decontaminated_log
    path("decontaminated_details.tsv"), emit: decontaminated_details
    path("notcontaminant.tsv")        , emit: notcontaminant_abundances, optional: true
    path("notcontaminant_counts.tsv") , emit: notcontaminant_counts, optional: true
    path("notcontaminant_log.tsv")    , emit: notcontaminant_log, optional: true
    path("notcontaminant_details.tsv"), emit: notcontaminant_details, optional: true
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args   ?: 'detailed=TRUE, normalize=TRUE'
    def method  = iscontaminant_method ?: "auto"
    def seed    = task.ext.seed ?: '100'
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(decontam))
    set.seed($seed) # Initialize random number generator for reproducibility

    # get controls and concentrations
    metadata <- read.csv("$metadata_tsv", header=TRUE, sep="\\t")
    metadata <- metadata[order(metadata[,1]),] # make sure samples in metadata and abundances are in same order
    if("quant_reading" %in% colnames(metadata)) {
        concentrations <- metadata\$quant_reading
    } else { concentrations <- list() }
    if("control" %in% colnames(metadata)) {
        negative_controls <- metadata\$control == 'control'
    } else { negative_controls <- list() }

    # get abundance table
    abundances_tsv <- read.csv("$abundances_tsv", header=TRUE, sep="\\t")
    abundances <- abundances_tsv[,-1]
    rownames(abundances) <- abundances_tsv[,1]
	abundances <- t(abundances)
	abundances <- abundances[order(rownames(abundances)),] # make sure samples in metadata and abundances are in same order

    # make sure that the metadata doesnt contain more samples than the abundance table (samples might have been lost before)
    metadata <- metadata[metadata[,1] %in% rownames(abundances),]

    # start dataframe for logging feature numbers
    df <- data.frame(input = nrow(abundances_tsv))

    # find contaminats: null hypothesis here is that sequences are **not** contaminants.
    # Requires sufficient positive proof an ASV is a contaminant before calling it so.
    # Uses prevalence, frequency, or a combination of both
    if(length(negative_controls) > 0 && length(concentrations) > 0) {
        contam_table <- isContaminant(as.matrix(abundances), conc=concentrations, neg=negative_controls, threshold=$iscontaminant_threshold, method='$method', $args)
    } else if(length(concentrations) > 0) {
        contam_table <- isContaminant(as.matrix(abundances), conc=concentrations, threshold=$iscontaminant_threshold, method='$method', $args)
    } else if(length(negative_controls) > 0) {
        contam_table <- isContaminant(as.matrix(abundances), neg=negative_controls, threshold=$iscontaminant_threshold, method='$method', $args)
    } else { stop("Neither negative controls nor concentration values were provided.") }
    contam_table <- cbind(ID = rownames(contam_table), contam_table)
    rownames(contam_table)[1] <- colnames(abundances_tsv)[1]
    write.table(contam_table, file = "decontaminated_details.tsv", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    # decontaminate abundance table
    iscontaminant_id <- contam_table[!(contam_table\$contaminant),][,1]
    abundances_filtered <- abundances_tsv[abundances_tsv[,1] %in% iscontaminant_id,]
    write.table(abundances_filtered, file = "decontaminated.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    # log retained features
    write.table(cbind(df, data.frame(isContaminant = nrow(abundances_filtered))), file = "decontaminated_log.tsv", row.names=FALSE, sep="\t")

    # log retained counts
    counts <- as.data.frame( t( rbind( colSums(abundances_tsv[-1]), colSums(abundances_filtered[-1]) ) ) )
    counts\$ID <- rownames(counts)
    colnames(counts) <- c("isContaminant_input","isContaminant_output", "sample")
    write.table(counts, file = "decontaminated_counts.tsv", row.names=FALSE, sep="\t")

    # Find non-contaminants: null hypothesis here is that sequences **are** contaminants.
    # Requires sufficient positive proof an ASV is not a contaminant before calling it so.
    # Only uses prevalence, i.e. method='prevalence' is fixed
    if(length(negative_controls) > 0) {
        # find non-contaminants
        notcontam_table <- isNotContaminant(as.matrix(abundances), neg=negative_controls, threshold=$isnotcontaminant_threshold, method='prevalence', $args)
        notcontam_table <- cbind(ID = rownames(notcontam_table), notcontam_table)
        rownames(notcontam_table)[1] <- colnames(abundances_tsv)[1]
        write.table( notcontam_table, file = "notcontaminant_details.tsv", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

        # non-contaminants abundance table
        isnotcontaminant_id <- notcontam_table[notcontam_table\$not.contaminant,][,1]
        abundances_filtered <- abundances_tsv[abundances_tsv[,1] %in% isnotcontaminant_id,]
        if(nrow(abundances_filtered) > 0) {
            write.table( abundances_filtered, file = "notcontaminant.tsv", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')
        }

        # log retained features
        write.table( cbind(df, data.frame(isNotContaminant = nrow(abundances_filtered))), file = "notcontaminant_log.tsv", row.names=FALSE, sep="\t")

        # log retained counts
        counts <- as.data.frame( t( rbind( colSums(abundances_tsv[-1]), colSums(abundances_filtered[-1]) ) ) )
        counts\$ID <- rownames(counts)
        colnames(counts) <- c("isNotContaminant_input","isNotContaminant_output", "sample")
        write.table(counts, file = "notcontaminant_counts.tsv", row.names=FALSE, sep="\t")
    }

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    decontam: ", packageVersion("decontam")) ), "versions.yml")
    """
}

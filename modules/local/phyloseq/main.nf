process PHYLOSEQ {
    tag "$prefix"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.50.0--r44hdfd78af_0' :
        'biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0' }"

    input:
    tuple val(prefix), path(tax_tsv), path(otu_tsv), path(sam_tsv), path(tree)

    output:
    tuple val(prefix), path("*phyloseq.rds"), emit: rds
    path "versions.yml"                     , emit: versions_phyloseq, topic: versions

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(phyloseq))

    otu_df  <- read.table("$otu_tsv", sep="\\t", header=TRUE, row.names=1)
    tax_df  <- read.table("$tax_tsv", sep="\\t", header=TRUE, row.names=1, comment.char = "")
    # per-rank confidences and the source database are not ranks, and phyloseq treats every column as one
    tax_df  <- tax_df[, !grepl("_confidence\$", colnames(tax_df)) & colnames(tax_df) != "database", drop=FALSE]
    otu_mat <- as.matrix(otu_df)
    tax_mat <- as.matrix(tax_df)
    storage.mode(otu_mat) <- "integer"

    OTU     <- otu_table(otu_mat, taxa_are_rows=TRUE)
    TAX     <- tax_table(tax_mat)
    phy_obj <- phyloseq(OTU, TAX)

    if (file.exists("$sam_tsv")) {
        sam_df  <- read.table("$sam_tsv", sep="\\t", header=TRUE, row.names=1)
        SAM     <- sample_data(sam_df)
        phy_obj <- merge_phyloseq(phy_obj, SAM)
    }

    if (file.exists("$tree")) {
        TREE    <- read_tree("$tree")
        phy_obj <- merge_phyloseq(phy_obj, TREE)
    }

    saveRDS(phy_obj, file = paste0("$prefix", "_phyloseq.rds"))

    # Version information
    writeLines(c("\\"${task.process}\\":",
        paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
        paste0("    phyloseq: ", packageVersion("phyloseq"))),
        "versions.yml"
    )
    """
}

process TREESUMMARIZEDEXPERIMENT {
    tag "$prefix"
    label 'process_low'

    conda "bioconda::bioconductor-treesummarizedexperiment=2.10.0 conda-forge::r-base=4.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-treesummarizedexperiment:2.10.0--r43hdfd78af_0' :
        'biocontainers/bioconductor-treesummarizedexperiment:2.10.0--r43hdfd78af_0' }"

    input:
    tuple val(prefix), path(tax_tsv), path(otu_tsv), path(sam_tsv), path(tree)

    output:
    tuple val(prefix), path("*TreeSummarizedExperiment.rds"), emit: rds
    path "versions.yml"                                     , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    suppressPackageStartupMessages(library(TreeSummarizedExperiment))

    # Read otu table. It must be in a SimpleList as a matrix where rows
    # represent taxa and columns samples.
    otu_mat  <- read.table("$otu_tsv", sep="\\t", header=TRUE, row.names=1)
    otu_mat <- as.matrix(otu_mat)
    assays <- SimpleList(counts = otu_mat)
    # Read taxonomy table. Correct format for it is DataFrame.
    taxonomy_table  <- read.table("$tax_tsv", sep="\\t", header=TRUE, row.names=1)
    taxonomy_table <- DataFrame(taxonomy_table)

    # Match rownames between taxonomy table and abundance matrix.
    taxonomy_table <- taxonomy_table[match(rownames(otu_mat), rownames(taxonomy_table)), ]

    # Create TreeSE object.
    tse <- TreeSummarizedExperiment(
        assays = assays,
        rowData = taxonomy_table
    )

    # If taxonomy table contains sequences, move them to referenceSeq slot
    if (!is.null(rowData(tse)[["sequence"]])) {
        referenceSeq(tse) <- DNAStringSet( rowData(tse)[["sequence"]] )
        rowData(tse)[["sequence"]] <- NULL
    }

    # If provided, we add sample metadata as DataFrame object. rownames of
    # sample metadata must match with colnames of abundance matrix.
    if (file.exists("$sam_tsv")) {
        sample_meta  <- read.table("$sam_tsv", sep="\\t", header=TRUE, row.names=1)
        sample_meta <- sample_meta[match(colnames(tse), rownames(sample_meta)), ]
        sample_meta  <- DataFrame(sample_meta)
        colData(tse) <- sample_meta
    }

    # If provided, we add phylogeny. The rownames in abundance matrix must match
    # with node labels in phylogeny.
    if (file.exists("$tree")) {
        phylogeny <- ape::read.tree("$tree")
        rowTree(tse) <- phylogeny
    }

    saveRDS(tse, file = paste0("$prefix", "_TreeSummarizedExperiment.rds"))

    # Version information
    writeLines(c("\\"${task.process}\\":",
        paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
        paste0("    TreeSummarizedExperiment: ", packageVersion("TreeSummarizedExperiment"))),
        "versions.yml"
    )
    """
}

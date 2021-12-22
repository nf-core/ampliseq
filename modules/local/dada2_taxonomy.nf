process DADA2_TAXONOMY {
    tag "${fasta},${database}"
    label 'process_high'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    path(fasta)
    path(database)
    val(outfile)

    output:
    path(outfile), emit: tsv
    path( "ASV_tax.rds" ), emit: rds
    path "versions.yml"  , emit: versions
    path "*.args.txt"    , emit: args

    script:
    def args = task.ext.args ?: ''
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    seq <- getSequences(\"$fasta\", collapse = TRUE, silence = FALSE)
    ranks = c('${params.dada_ref_databases[params.dada_ref_taxonomy]["ranks"].join("', '")}')
    taxa <- assignTaxonomy(
        seq,
        \"$database\",
        taxLevels = ranks,
        $args,
        multithread = $task.cpus,
        verbose=TRUE,
        outputBootstraps = TRUE
    )

    # Make a data frame, add ASV_ID from seq, set confidence to the bootstrap for the most specific taxon and reorder columns before writing to file
    tx <- data.frame(ASV_ID = names(seq), taxa, sequence = row.names(taxa\$tax), row.names = names(seq))

    tx\$confidence <- 0
    for ( r in rev(ranks) ) {
        tx\$confidence <- ifelse(tx\$confidence == 0 & !is.na(tx[[sprintf("tax.%s", r)]]), tx[[sprintf("boot.%s", r)]]/100, tx\$confidence)
    }

    taxa_export <- data.frame(ASV_ID = tx\$ASV_ID)
    for ( r in ranks ) {
        taxa_export[[r]]    <- tx[[sprintf("tax.%s", r)]]
    }
    taxa_export\$confidence <- tx\$confidence
    taxa_export\$sequence   <- tx\$sequence 
    row.names(taxa_export)  <- tx\$sequence

    write.table(taxa_export, file = \"$outfile\", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    # Save a version with rownames for addSpecies
    saveRDS(taxa_export, "ASV_tax.rds")

    write.table('assignTaxonomy\t$args', file = "assignTaxonomy.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

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
    taxa <- assignTaxonomy(seq, \"$database\", taxLevels = c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), $args, multithread = $task.cpus, verbose=TRUE, outputBootstraps = TRUE)

    # Make a data frame, add ASV_ID from seq, set confidence to the bootstrap for the most specific taxon and reorder columns before writing to file
    tx <- data.frame(ASV_ID = names(seq), taxa, sequence = row.names(taxa\$tax), row.names = names(seq))
    tx\$confidence <- with(tx,
        ifelse(!is.na(tax.Genus), boot.Genus,
            ifelse(!is.na(tax.Family), boot.Family,
                ifelse(!is.na(tax.Order), boot.Order,
                    ifelse(!is.na(tax.Class), boot.Class,
                        ifelse(!is.na(tax.Phylum), boot.Phylum,
                            ifelse(!is.na(tax.Kingdom), boot.Kingdom,
                                ifelse(!is.na(tax.Domain), boot.Domain, 0)
                            )
                        )
                    )
                )
            )
        )
    )/100
    taxa_export <- data.frame(
        ASV_ID = tx\$ASV_ID,
        Domain = tx\$tax.Domain,
        Kingdom = tx\$tax.Kingdom,
        Phylum = tx\$tax.Phylum,
        Class = tx\$tax.Class,
        Order = tx\$tax.Order,
        Family = tx\$tax.Family,
        Genus = tx\$tax.Genus,
        confidence = tx\$confidence,
        sequence = tx\$sequence,
        row.names = names(seq)
    )

    write.table(taxa_export, file = \"$outfile\", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    # Save a version with rownames for addSpecies
    taxa_export <- cbind( ASV_ID = tx\$ASV_ID, taxa\$tax, confidence = tx\$confidence)
    saveRDS(taxa_export, "ASV_tax.rds")

    write.table('assignTaxonomy\t$args', file = "assignTaxonomy.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

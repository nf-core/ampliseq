process DADA2_TAXONOMY {
    tag "${fasta},${database}"
    label 'process_high'

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data' :
        'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    path(fasta)
    path(database)
    val(outfile)
    val(taxlevels_input)

    output:
    path("*${outfile}.tsv")    , emit: tsv
    path( "*${outfile}.rds" ), emit: rds
    path "versions.yml"  , emit: versions
    path "*.args.txt"    , emit: args

    script:
    def args = task.ext.args ?: ''
    def taxlevels = taxlevels_input ?
        'c("' + taxlevels_input.split(",").join('","') + '")' :
        'c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")'
    def seed = task.ext.seed ?: '100'
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed($seed) # Initialize random number generator for reproducibility

    taxlevels <- $taxlevels

    seq <- getSequences(\"$fasta\", collapse = TRUE, silence = FALSE)
    taxa <- assignTaxonomy(seq, \"$database\", taxLevels = taxlevels, $args, multithread = $task.cpus, verbose=TRUE, outputBootstraps = TRUE)

    # (1) Make a data frame, add ASV_ID from seq
    tx <- data.frame(ASV_ID = names(seq), taxa, sequence = row.names(taxa\$tax), row.names = names(seq))

    # (2) Set confidence to the bootstrap for the most specific taxon
    # extract columns with taxonomic values
    tax <- tx[ , grepl( "tax." , names( tx ) ) ]
    # find first occurrence of NA
    res <- max.col(is.na(tax), ties = "first")
    # correct if no NA is present in column to NA
    if(any(res == 1)) is.na(res) <- (res == 1) & !is.na(tax[[1]])
    # find index of last entry before NA, which is the bootstrap value
    res <- res-1
    # if NA choose last entry
    res[is.na(res)] <- ncol(tax)
    # extract bootstrap values
    boot <- tx[ , grepl( "boot." , names( tx ) ) ]
    boot\$last_tax <- res
    valid_boot <- apply(boot,1,function(x) x[x[length(x)]][1]/100 )
    # replace missing bootstrap values (NA) with 0
    valid_boot[is.na(valid_boot)] <- 0
    # add bootstrap values to column confidence
    tx\$confidence <- valid_boot

    # (3) Reorder columns before writing to file
    expected_order <- c("ASV_ID",paste0("tax.",taxlevels),"confidence","sequence")
    expected_order <- intersect(expected_order,colnames(tx))
    taxa_export <- subset(tx, select = expected_order)
    colnames(taxa_export) <- sub("tax.", "", colnames(taxa_export))
    rownames(taxa_export) <- names(seq)

    write.table(taxa_export, file = \"${fasta.baseName}${outfile}.tsv\", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    # Save a version with rownames for addSpecies
    taxa_export <- cbind( ASV_ID = tx\$ASV_ID, taxa\$tax, confidence = tx\$confidence)
    saveRDS(taxa_export, "${fasta.baseName}${outfile}.rds")

    write.table('assignTaxonomy\t$args\ntaxlevels\t$taxlevels\nseed\t$seed', file = "assignTaxonomy.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

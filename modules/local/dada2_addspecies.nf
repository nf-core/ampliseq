process DADA2_ADDSPECIES {
    tag "${taxtable},${database}"
    label 'process_cpu_single'
    label 'process_medium_memory'
    label 'process_long'

    conda "bioconda::bioconductor-dada2=1.38.0 conda-forge::r-base=4.5.2 conda-forge::r-digest=0.6.39 conda-forge::tbb=2022.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/81/81153df5d53322e6d91b2c4c9bc4da50774fb1d101ead002fe75bb75fc6f036c/data' :
        'community.wave.seqera.io/library/bioconductor-dada2_r-base_r-digest_tbb:38acac09bac46f36' }"

    input:
    path(taxtable)
    path(database)
    val(taxlevels_input)

    output:
    path("${taxtable.baseName}.species.tsv")   , emit: tsv
    path "versions.yml" , emit: versions
    path "*.args.txt"   , emit: args


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

    #add "Species" if not already in taxlevels
    taxlevels <- $taxlevels
    if ( !"Species" %in% taxlevels ) { taxlevels <- c(taxlevels,"Species") }

    taxtable <- readRDS(\"$taxtable\")

    #remove Species annotation from assignTaxonomy
    taxa_nospecies <- taxtable[,!colnames(taxtable) %in% 'Species']

    tx <- addSpecies(taxa_nospecies, \"$database\", $args, verbose=TRUE)

    # Create a table with specified column order
    tmp <- data.frame(row.names(tx)) # To separate ASV_ID from sequence
    expected_order <- c("ASV_ID",taxlevels,"confidence")
    taxa <- as.data.frame( subset(tx, select = expected_order) )
    taxa\$sequence <- tmp[,1]
    row.names(taxa) <- row.names(tmp)

    #rename Species annotation to Species_exact
    colnames(taxa)[which(names(taxa) == "Species")] <- "Species_exact"

    #add Species annotation from assignTaxonomy again, after "Genus" column
    if ( "Species" %in% colnames(taxtable) ) {
        taxtable <- data.frame(taxtable)
        taxa_export <- data.frame(append(taxa, list(Species=taxtable\$Species), after=match("Genus", names(taxa))))
    } else {
        taxa_export <- taxa
    }

    write.table(taxa_export, file = \"${taxtable.baseName}.species.tsv\", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    write.table('addSpecies\t$args\ntaxlevels\t$taxlevels\nseed\t$seed', file = "addSpecies.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

process DADA2_ADDSPECIES {
    tag "${taxtable},${database}"
    label 'process_high'
    label 'single_cpu'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    path(taxtable)
    path(database)
    val(outfile)
    val(taxlevels_input)

    output:
    path(outfile)       , emit: tsv
    path "versions.yml" , emit: versions
    path "*.args.txt"   , emit: args

    when:
    task.ext.when == null || task.ext.when

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

    write.table(taxa_export, file = \"$outfile\", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    write.table('addSpecies\t$args\ntaxlevels\t$taxlevels\nseed\t$seed', file = "addSpecies.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

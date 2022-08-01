process DADA2_ADDSPECIES {
    tag "${taxtable},${database}"
    label 'process_high'
    label 'single_cpu'

    conda (params.enable_conda ? "bioconductor-dada2=1.22.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.22.0--r41h399db7b_0' :
        'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

    input:
    path(taxtable)
    path(database)
    val(outfile)

    output:
    path(outfile)       , emit: tsv
    path "versions.yml" , emit: versions
    path "*.args.txt"   , emit: args

    script:
    def args = task.ext.args ?: ''
    def taxlevels = task.ext.taxlevels ? 
        'c("' + task.ext.taxlevels.split(",").join('","') + '")' : 
        'c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")'
    def seed = task.ext.seed ?: '100'
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed($seed) # Initialize random number generator for reproducibility

    #add "Species" if not already in taxlevels
    taxlevels <- $taxlevels
    if ( !"Species" %in% taxlevels ) { taxlevels <- c(taxlevels,"Species") }

    taxtable <- readRDS(\"$taxtable\")

    tx <- addSpecies(taxtable, \"$database\", $args, verbose=TRUE)

    # Create a table with specified column order
    tmp <- data.frame(row.names(tx)) # To separate ASV_ID from sequence
    expected_order <- c("ASV_ID",taxlevels,"confidence")
    taxa <- as.data.frame( subset(tx, select = expected_order) )
    taxa\$sequence <- tmp[,1]
    row.names(taxa) <- row.names(tmp)

    write.table(taxa, file = \"$outfile\", sep = "\\t", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '')

    write.table('addSpecies\t$args\ntaxlevels\t$taxlevels\nseed\t$seed', file = "addSpecies.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

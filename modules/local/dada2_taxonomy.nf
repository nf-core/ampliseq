// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process DADA2_TAXONOMY {
    tag "${fasta},${database}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconductor-dada2=1.18.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    } else {
        container "quay.io/biocontainers/bioconductor-dada2:1.18.0--r40h5f743cb_0"
    }

    input:
    path(fasta)
    path(database)
    val(outfile)
    
    output:
    path(outfile), emit: tsv
    path( "ASV_tax.rds" ), emit: rds
    path "*.version.txt" , emit: version
    path "*.args.txt"    , emit: args

    script:
    def software      = getSoftwareName(task.process)
    """
    #!/usr/bin/env Rscript
    suppressPackageStartupMessages(library(dada2))
    set.seed(100) # Initialize random number generator for reproducibility

    seq <- getSequences(\"$fasta\", collapse = TRUE, silence = FALSE)
    taxa <- assignTaxonomy(seq, \"$database\", $options.args, multithread = $task.cpus, verbose=TRUE, outputBootstraps = TRUE)

    # Make a data frame, add ASV_ID from seq, set confidence to the bootstrap for the most specific taxon and reorder columns before writing to file
    tx <- data.frame(ASV_ID = names(seq), taxa, sequence = row.names(taxa\$tax), row.names = names(seq))
    tx\$confidence <- with(tx, 
        ifelse(!is.na(tax.Genus), boot.Genus, 
            ifelse(!is.na(tax.Family), boot.Family,
                ifelse(!is.na(tax.Order), boot.Order,
                    ifelse(!is.na(tax.Class), boot.Class,
                        ifelse(!is.na(tax.Phylum), boot.Phylum,
                            ifelse(!is.na(tax.Kingdom), boot.Kingdom, 0)
                        )
                    )
                )
            )
        )
    )/100
    taxa_export <- data.frame(
        ASV_ID = tx\$ASV_ID,
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

    write.table(taxa_export, file = \"$outfile\", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

    # Save a version with rownames for addSpecies
    taxa_export <- cbind( ASV_ID = tx\$ASV_ID, taxa\$tax, confidence = tx\$confidence )
    saveRDS(taxa_export, "ASV_tax.rds")

    write.table('assignTaxonomy\t$options.args', file = "assignTaxonomy.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(packageVersion("dada2"), file = "${software}.version.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
    """
}

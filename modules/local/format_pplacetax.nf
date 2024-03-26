process FORMAT_PPLACETAX {
    tag "${tax.baseName}"
    label 'process_high'

    conda "bioconda::bioconductor-dada2=1.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.30.0--r43hf17093f_0' :
        'biocontainers/bioconductor-dada2:1.30.0--r43hf17093f_0' }"

    input:
    tuple val(meta), path(tax)

    output:
    path("*.per_query_unique.tsv"), emit: unique
    path("*.taxonomy.tsv")        , emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    tax = read.table("$tax", header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = '', quote = '')

    df <- data.frame(
        name=character(),
        LWR=character(),
        fract=character(),
        aLWR=character(),
        afract=character(),
        taxopath=character(),
        stringsAsFactors=FALSE)

    for (asvid in unique(tax\$name)) {
        # subset per ASV ID
        temp = subset(tax, tax\$name == asvid)
        maxLWR = max(temp\$LWR)
        temp = subset(temp, temp\$LWR == maxLWR)

        if ( nrow(temp) == 1 ) {
            # STEP 1: if there is just 1 maximum, choose this row
            print( paste ( asvid,"was added in STEP 1" ) )
            temp = temp[which.max(temp\$LWR),]
            df <- rbind(df, temp)
        } else {
            # STEP2: be conservative, choose the taxonomy with less entries
            # count number of entries in "taxopath", simplified: number of semicolons
            temp\$taxlevels = nchar( temp\$taxopath ) - nchar(gsub(';', '', temp\$taxopath )) + 1
            minTaxlevels = min(temp\$taxlevels)
            temp = subset(temp, temp\$taxlevels == minTaxlevels)
            temp\$taxlevels = NULL
            if ( nrow(temp) == 1 ) {
                df <- rbind(df, temp)
                print( paste ( asvid,"was added in STEP 2" ) )
            } else {
                # STEP 3: if number of entries is same, remove last entry
                # then, check if reduced entries are identical > if yes, choose any row, if no, repeat
                # at that step the taxonomies have same length
                print( paste ( asvid,"enters STEP 3" ) )
                list_taxonpath <- str_split( temp\$taxopath, ";")
                df_taxonpath <- as.data.frame(do.call(rbind, list_taxonpath))
                for (i in ncol(df_taxonpath):0) {
                    # choose first column and change taxon to reduced overlap
                    if ( length(unique(df_taxonpath[,i])) == 1 ) {
                        if (i>1) {
                            temp\$taxopath <- apply(df_taxonpath[1,1:i], 1, paste, collapse=";")[[1]]
                        } else if (i==1) {
                            temp\$taxopath <- df_taxonpath[1,1]
                        }
                        df <- rbind(df, temp[1,])
                        print( paste ( asvid,"was added with",temp\$taxopath[[1]],"in STEP 3a" ) )
                        break
                    } else if (i==0) {
                        # if all fails, i.e. there is no consensus on any taxonomic level, use NA
                        temp\$taxopath <- "NA"
                        df <- rbind(df, temp[1,])
                        print( paste ( asvid,"was added with",temp\$taxopath[[1]],"in STEP 3b" ) )
                    }
                }
            }
        }
    }

    # output the cleaned, unified, unique taxonomic classification per ASV
    write.table(df, file = "${meta.id}.taxonomy.per_query_unique.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep = "\\t")

    # output ASV taxonomy table for QIIME2
    df_clean <- subset(df, select = c("name","taxopath"))
    colnames(df_clean) <- c("ASV_ID","taxonomy")
    write.table(df_clean, file = "${meta.id}.taxonomy.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep = "\\t")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
    """
}

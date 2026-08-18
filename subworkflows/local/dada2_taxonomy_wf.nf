/*
 * Taxonomic classification with DADA2
 */

include { CUTADAPT as CUTADAPT_TAXONOMY  } from '../../modules/nf-core/cutadapt/main'
include { VSEARCH_USEARCHGLOBAL          } from '../../modules/nf-core/vsearch/usearchglobal/main'

include { DADA2_TAXONOMY                               } from '../../modules/local/dada2_taxonomy'
include { DADA2_ADDSPECIES                             } from '../../modules/local/dada2_addspecies'
include { FORMAT_TAXRESULTS as FORMAT_TAXRESULTS_STD   } from '../../modules/local/format_taxresults'
include { FORMAT_TAXRESULTS as FORMAT_TAXRESULTS_ADDSP } from '../../modules/local/format_taxresults'
include { ASSIGNSH                                     } from '../../modules/local/assignsh'

include { makeComplement                 } from '../../subworkflows/local/utils_nfcore_ampliseq_pipeline'

// sanitize a database key into the filesystem/tag-safe name used throughout this repo
def sanitize(key) {
    key.replace('=','_').replace('.','_')
}

workflow DADA2_TAXONOMY_WF {
    take:
    ch_assigntax                // tuple val(db_key), path(assigntax)   -- one entry per listed --dada_ref_taxonomy database
    ch_addspecies                // tuple val(db_key), path(addspecies) -- one entry per listed database
    val_dada_ref_taxonomy_list   // List<String>, raw (un-sanitized) database keys, in the order they were listed
    ch_fasta                     // ASV fasta, the same for every listed database
    ch_full_fasta                // full (uncut) ASV fasta, the same for every listed database
    val_dada_assign_chunksize

    main:
    // reverse lookup: sanitized name -> original database key (needed to recover the tag lost by collectFile)
    def dbKeyForSanitized = val_dada_ref_taxonomy_list.collectEntries { key -> [ (sanitize(key)): key ] }

    // Set cutoff to use for SH assignment and path(s) to each listed database's SH taxonomy file(s)
    if ( params.addsh ) {
        vsearch_cutoff = 0.985
        ch_shinfo = ch_assigntax
            .map { db_key, _db -> db_key }
            .flatMap { db_key -> params.dada_ref_databases[db_key]["shfile"].collect { url -> [ db_key, file(url) ] } }
            .groupTuple(by: 0)
    }

    //cut taxonomy to expected amplicon
    if (params.cut_dada_ref_taxonomy) {
        ch_assigntax =
            ch_assigntax
                .map {
                    db_key, db ->
                        def meta = [:]
                        meta.single_end = true
                        meta.id = "assignTaxonomy_${sanitize(db_key)}"
                        meta.db_key = db_key
                        meta.primer_fwd = params.primer_fwd
                        meta.primer_rev_revcomp = makeComplement ( "${params.primer_rev}".reverse() )
                        [ meta, db ] }
        ch_assigntax =
            CUTADAPT_TAXONOMY ( ch_assigntax ).reads
                .map { meta, db -> [ meta.db_key, db ] }
    }

    //set file name prefix
    if (params.cut_its == "none") {
        ASV_tax_name = "ASV_tax"
    } else {
        ASV_tax_name = "ASV_ITS_tax"
    }

    //split sequences into chunks (shared across every listed database)
    ch_fasta_chunks = ch_fasta.splitFasta( by: val_dada_assign_chunksize, file: true )

    //per-database outfile suffix + taxlevels, derived once per listed database
    ch_assigntax_taxlevels = ch_assigntax
        .map {
            db_key, db ->
                def taxlevels = params.dada_assign_taxlevels ? "${params.dada_assign_taxlevels}" :
                    ( params.dada_ref_databases[db_key]["taxlevels"] ?: "" )
                [ db_key, db, ".${ASV_tax_name}.${sanitize(db_key)}", taxlevels ]
        }

    //DADA2 assignTaxonomy -- one task per (database x fasta chunk)
    DADA2_TAXONOMY (
        ch_assigntax_taxlevels
            .combine( ch_fasta_chunks )
            .map { db_key, db, outfile, taxlevels, fasta -> [ db_key, fasta, db, outfile, taxlevels ] }
    )

    // collect all DADA2_TAXONOMY.out.tsv chunks into one file per database
    ch_dada2_taxonomy_tsv = DADA2_TAXONOMY.out.tsv
        .collectFile(newLine: false, cache: true, keepHeader: true, skip: 1, sort: true) { db_key, tsv ->
            [ "${ASV_tax_name}.${sanitize(db_key)}.tsv", tsv ]
        }
        .map { f -> [ dbKeyForSanitized[ (f.name - "${ASV_tax_name}." - ".tsv") ], f ] }
    ch_dada2_taxonomy_tsv.subscribe{ _db_key, f -> file(f).copyTo("${params.outdir}/dada2") }

    if (params.cut_its != "none") {
        FORMAT_TAXRESULTS_STD (
            ch_dada2_taxonomy_tsv
                .combine( ch_full_fasta.collect() )
                .map { db_key, tsv, full_fasta -> [ db_key, tsv, full_fasta, "ASV_tax.${sanitize(db_key)}.tsv" ] }
        )
    }

    //DADA2 addSpecies
    if (!params.skip_dada_addspecies) {
        DADA2_ADDSPECIES (
            DADA2_TAXONOMY.out.rds
                .combine( ch_addspecies, by: 0 )
                .combine( ch_assigntax_taxlevels.map { db_key, _db, _outfile, taxlevels -> [ db_key, taxlevels ] }, by: 0 )
        )

        // collect all DADA2_ADDSPECIES.out.tsv chunks into one file per database
        ch_dada2_addspecies_tsv = DADA2_ADDSPECIES.out.tsv
            .collectFile(newLine: false, cache: true, keepHeader: true, skip: 1, sort: true) { db_key, tsv ->
                [ "${ASV_tax_name}_species.${sanitize(db_key)}.tsv", tsv ]
            }
            .map { f -> [ dbKeyForSanitized[ (f.name - "${ASV_tax_name}_species." - ".tsv") ], f ] }
        ch_dada2_addspecies_tsv.subscribe{ _db_key, f -> file(f).copyTo("${params.outdir}/dada2") }

        if (params.cut_its == "none") {
            ch_dada2_tax1 = ch_dada2_addspecies_tsv
        } else {
            FORMAT_TAXRESULTS_ADDSP (
                ch_dada2_addspecies_tsv
                    .combine( ch_full_fasta.collect() )
                    .map { db_key, tsv, full_fasta -> [ db_key, tsv, full_fasta, "ASV_tax_species.${sanitize(db_key)}.tsv" ] }
            )
            ch_dada2_tax1 = FORMAT_TAXRESULTS_ADDSP.out.tsv
        }
    //no DADA2 addSpecies, use results from assignTaxonomy:
    } else {
        if (params.cut_its == "none") {
            ch_dada2_tax1 = ch_dada2_taxonomy_tsv
        } else {
            ch_dada2_tax1 = FORMAT_TAXRESULTS_STD.out.tsv
        }
    }

    //if addsh set: add SH assignments
    if ( params.addsh ) {
        //set file name prefix for SH assignments
        if (!params.skip_dada_addspecies) {
            ASV_SH_name = "ASV_tax_species_SH"
        } else {
            ASV_SH_name = "ASV_tax_SH"
        }

        // one VSEARCH_USEARCHGLOBAL run per listed database, querying the same ASV fasta against that database's own reference
        ch_vsearch_input = ch_assigntax
            .combine( ch_fasta.collect() )
            .map { db_key, db, fasta -> [ [ id: "${ASV_tax_name}.vsearch.${sanitize(db_key)}", db_key: db_key ], fasta, db ] }

        VSEARCH_USEARCHGLOBAL(
            ch_vsearch_input.map { meta, fasta, _db -> [ meta, fasta ] },
            ch_vsearch_input.map { _meta, _fasta, db -> db },
            vsearch_cutoff, 'blast6out', ""
        )
        ch_blastfile = VSEARCH_USEARCHGLOBAL.out.txt.map { meta, blastfile -> [ meta.db_key, blastfile ] }

        ASSIGNSH(
            ch_dada2_tax1
                .join( ch_shinfo )
                .join( ch_blastfile )
                .map { db_key, tax_tsv, sh_files, blastfile -> [ db_key, tax_tsv, sh_files, blastfile, "${ASV_SH_name}.${sanitize(db_key)}.tsv" ] }
        )
        ch_dada2_tax = ASSIGNSH.out.tsv
    } else {
        ch_dada2_tax = ch_dada2_tax1
    }

    emit:
    cut_tax  = params.cut_dada_ref_taxonomy ? CUTADAPT_TAXONOMY.out.log : [[],[]]
    tax      = ch_dada2_tax
}

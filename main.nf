#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/ampliseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/ampliseq
    Website: https://nf-co.re/ampliseq
    Slack  : https://nfcore.slack.com/channels/ampliseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AMPLISEQ                } from './workflows/ampliseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_ampliseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_ampliseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_AMPLISEQ {

    take:
    ch_metadata            // channel: [ path(metadata) ] or empty
    ch_report_template     // channel: [ path(report_template) ]
    ch_report_css          // channel: [ path(report_css) ]
    ch_report_logo         // channel: [ path(report_logo) ]
    ch_report_abstract     // channel: [ path(report_abstract) ] or []
    ch_pplace_sheet        // channel: parsed pplace_sheet rows, or empty
    ch_expected_sequences  // channel: [ path(expected_sequences) ] or empty
    ch_expected_abundances // channel: [ path(expected_abundances) ] or empty
    ch_expected_profile    // channel: [ path(expected_profile) ] or empty
    ch_metadata_category   // channel: tokenized metadata_category, or empty

    main:

    //
    // WORKFLOW: Run pipeline
    //
    AMPLISEQ (
        // samplesheet, // currently the samplesheet is parsed in workflows/ampliseq.nf
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
        ch_metadata,
        ch_report_template,
        ch_report_css,
        ch_report_logo,
        ch_report_abstract,
        ch_pplace_sheet,
        ch_expected_sequences,
        ch_expected_abundances,
        ch_expected_profile,
        ch_metadata_category,
    )
    emit:
    multiqc_report = AMPLISEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params {
    input: Path?
    input_fasta: Path?
    input_folder: Path?
    primer_fwd: String?
    primer_rev: String?
    metadata: Path?
    multiregion: Path?
    outdir: String
    ref_taxonomy_storage: String?
    save_intermediates: Boolean
    email: String?
    sequencing_type: String = 'illumina_pe'
    binned_quality: String?
    ignore_binned_quality: Boolean
    illumina_pe_readthrough: Boolean
    quality_type: String = 'Auto'
    multiple_sequencing_runs: Boolean
    input_folder_extensions: String = '/*_R{1,2}_001.fastq.gz'
    min_read_counts: Integer = 1
    ignore_empty_input_files: Boolean
    retain_untrimmed: Boolean
    cutadapt_min_overlap: Integer = 3
    cutadapt_max_error_rate: Float = 0.1
    double_primer: Boolean
    ignore_failed_trimming: Boolean
    truncq: Integer = 2
    trunclenf: Integer?
    trunclenr: Integer?
    trunc_qmin: Integer = 25
    trunc_rmin: Float = 0.75
    max_ee: Integer = 2
    min_len: Integer = 50
    max_len: Integer?
    ignore_failed_filtering: Boolean
    asv_calling: String = 'auto'
    sample_inference: String = 'pooled'
    savont_options: String = '--fl-16s'
    mergepairs_strategy: String = 'merge'
    mergepairs_consensus_match: Integer = 1
    mergepairs_consensus_mismatch: Integer = -2
    mergepairs_consensus_gap: Integer = -4
    mergepairs_consensus_minoverlap: Integer = 12
    mergepairs_consensus_maxmismatch: Integer = 0
    mergepairs_consensus_percentile_cutoff: Float = 0.001
    vsearch_cluster: Boolean
    vsearch_cluster_id: Float = 0.97
    raise_filter_stacksize: Boolean = true
    decontam: String = 'none'
    decontam_decontaminate_method: String = 'auto'
    decontam_decontaminate_threshold: Float = 0.1
    decontam_notcontaminant_threshold: Float = 0.5
    filter_ssu: String?
    min_len_asv: Integer?
    max_len_asv: Integer?
    filter_codons: Boolean
    orf_start: Integer = 1
    orf_end: Integer?
    stop_codons: String = 'TAA,TAG'
    dada_ref_taxonomy: String = 'sbdi-gtdb=R11-RS232-1'
    dada_ref_tax_custom: Path?
    dada_ref_tax_custom_sp: Path?
    dada_min_boot: Integer = 50
    dada_assign_taxlevels: String?
    cut_dada_ref_taxonomy: Boolean
    consolidate_taxonomies: String?
    dada_addspecies_allowmultiple: Boolean
    dada_taxonomy_rc: Boolean
    dada_assign_chunksize: Integer = 10000
    pplace_sheet: Path?
    pplace_tree: Path?
    pplace_aln: Path?
    pplace_model: String?
    pplace_alnmethod: String = 'clustalo'
    pplace_taxonomy: Path?
    pplace_name: String?
    qiime_ref_taxonomy: String?
    qiime_ref_tax_custom: String?
    qiime_classifier: Path?
    kraken2_ref_taxonomy: String?
    kraken2_ref_tax_custom: Path?
    kraken2_assign_taxlevels: String?
    kraken2_confidence: Float = 0.0
    sintax_ref_taxonomy: String?
    sintax_ref_tax_custom: Path?
    sintax_assign_taxlevels: String?
    vsearch_lca_ref_tax_custom: Path?
    vsearch_lca_ref_taxonomy: String?
    vsearch_lca_assign_taxlevels: String?
    vsearch_lca_id: Float = 0.9
    vsearch_lca_maxaccepts: Integer = 0
    vsearch_lca_maxrejects: Integer = 0
    vsearch_lca_lca_cutoff: Float = 0.9
    vsearch_lca_query_cov: Float = 1.0
    addsh: Boolean
    cut_its: String = 'none'
    its_partial: Integer = 0
    its_extractor: String = 'itsx'
    sidle_ref_taxonomy: String?
    sidle_ref_tax_custom: Path?
    sidle_ref_seq_custom: Path?
    sidle_ref_aln_custom: Path?
    sidle_ref_tree_custom: Path?
    sidle_ref_degenerates: Integer = 5
    sidle_ref_cleaning: String?
    exclude_taxa: String = 'mitochondria,chloroplast'
    min_frequency: Integer = 1
    min_samples: Integer = 1
    metadata_category: String?
    metadata_category_barplot: String?
    qiime_adonis_formula: String?
    picrust: Boolean
    sbdiexport: Boolean
    diversity_rarefaction_depth: Integer = 500
    tax_agglom_min: Integer = 2
    tax_agglom_max: Integer = 6
    ancom_sample_min_count: Integer = 1
    ancom: Boolean
    ancombc: Boolean
    ancombc_formula: String?
    ancombc_formula_reflvl: String?
    ancombc_effect_size: Float = 1.0
    ancombc_significance: Float = 0.05
    ancombc2: Boolean
    ancombc2_formula: String?
    ancombc2_formula_reflvl: String?
    expected_sequences: Path?
    expected_sequences_region: String = 'observed'
    expected_sequences_mismatches: Integer = 0
    expected_sequences_merge: String = 'all'
    expected_abundances: Path?
    expected_profile: Path?
    report_template: String = '${projectDir}/assets/report_template.Rmd'
    report_css: String = '${projectDir}/assets/nf-core_style.css'
    report_logo: String = '${projectDir}/assets/nf-core-ampliseq_logo_light_long.png'
    report_title: String = 'Summary of analysis results'
    report_abstract: Path?
    run_pplace: Boolean
    skip_fastqc: Boolean
    skip_porechop_abi: Boolean
    skip_chopper: Boolean
    skip_cutadapt: Boolean
    skip_dada_quality: Boolean
    skip_barrnap: Boolean
    skip_qiime: Boolean
    skip_qiime_downstream: Boolean
    skip_taxonomy: Boolean
    skip_dada_taxonomy: Boolean
    skip_dada_addspecies: Boolean
    skip_barplot: Boolean
    skip_abundance_tables: Boolean
    skip_alpha_rarefaction: Boolean
    skip_diversity_indices: Boolean
    skip_phyloseq: Boolean
    skip_tse: Boolean
    skip_multiqc: Boolean
    skip_report: Boolean
    skip_parquet_summary: Boolean
    seed: Integer = 100
    version: Boolean
    publish_dir_mode: String
    email_on_fail: String?
    plaintext_email: Boolean
    max_multiqc_email_size: String = '25.MB'
    monochrome_logs: Boolean
    multiqc_config: Path?
    multiqc_logo: Path?
    multiqc_methods_description: Path?
    validate_params: Boolean = true
    pipelines_testdata_base_path: String
    trace_report_suffix: String
    help: Boolean
    help_full: Boolean
    show_hidden: Boolean
    custom_config_version: String
    custom_config_base: String
    config_profile_name: String
    config_profile_description: String
    config_profile_contact: String
    config_profile_url: String
    multiqc_title: String?
}

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        [],
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_AMPLISEQ (
        PIPELINE_INITIALISATION.out.metadata,
        PIPELINE_INITIALISATION.out.report_template,
        PIPELINE_INITIALISATION.out.report_css,
        PIPELINE_INITIALISATION.out.report_logo,
        PIPELINE_INITIALISATION.out.report_abstract,
        PIPELINE_INITIALISATION.out.pplace_sheet,
        PIPELINE_INITIALISATION.out.expected_sequences,
        PIPELINE_INITIALISATION.out.expected_abundances,
        PIPELINE_INITIALISATION.out.expected_profile,
        PIPELINE_INITIALISATION.out.metadata_category,
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_AMPLISEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

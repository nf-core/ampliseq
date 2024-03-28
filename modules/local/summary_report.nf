process SUMMARY_REPORT  {
    label 'process_low'

    conda "conda-forge::r-base=4.2.3 conda-forge::r-rmarkdown=2.22 conda-forge::r-tidyverse=2.0.0 conda-forge::r-knitr=1.43 conda-forge::r-dt=0.28 conda-forge::r-dtplyr=1.3.1 conda-forge::r-formattable=0.2.1 conda-forge::r-purrr=1.0.1 conda-forge::r-vegan=2.6_4 conda-forge::r-optparse=1.7.3 conda-forge::r-ggplot2=3.4.2 conda-forge::r-dplyr=1.1.2 conda-forge::r-data.table=1.14.8 conda-forge::r-patchwork=1.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    path(report_template)
    path(report_styles)
    path(report_logo)
    path(report_abstract)
    path(metadata)
    path(input_samplesheet)
    path(input_fasta)
    path(mqc_plots)
    path(cutadapt_summary)
    val(find_truncation_values)
    path(dada_filtntrim_args)
    path(dada_qual_stats)
    path(dada_pp_qual_stats)
    tuple val(meta), path(dada_err_svgs)
    path(dada_asv_table)
    path(dada_asv_fa)
    path(dada_tab)
    path(dada_stats)
    path(vsearch_cluster)
    path(barrnap_summary)
    path(filter_ssu_stats)
    path(filter_ssu_asv)
    path(filter_len_asv_stats)
    path(filter_len_asv_len_orig)
    path(filter_codons_fasta)
    path(filter_codons_stats)
    path(itsx_cutasv_summary)
    path(dada2_tax)
    tuple val(meta_ref), path(cut_dada_ref_taxonomy) // cutadapt log when params.cut_dada_ref_taxonomy
    path(sintax_tax)
    path(kraken2_tax)
    path(pplace_tax)
    tuple val(meta_pplace), path(pplace_heattree)
    path(qiime2_tax)
    val(run_qiime2)
    val(val_used_taxonomy)
    val(qiime2_filtertaxa) // <ASV count original +1>,<ASV count filtered +2>
    path(filter_stats_tsv)
    path(barplot)
    path(abundance_tables, stageAs: 'abundance_tables/*')
    val(alpha_rarefaction)
    path(diversity_indices)
    path(diversity_indices_alpha, stageAs: 'alpha_diversity/*') // prevent folder name collisons
    path(diversity_indices_beta, stageAs: 'beta_diversity/*') // prevent folder name collisons
    path(diversity_indices_adonis, stageAs: 'beta_diversity/adonis/*') // prevent folder name collisons
    path(ancom)
    path(picrust_pathways)
    path(sbdi, stageAs: 'sbdi/*')
    path(phyloseq, stageAs: 'phyloseq/*')

    output:
    path "*.svg"               , emit: svg, optional: true
    path "summary_report.html" , emit: report
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // make named R list (comma separated)
    // all non-boolean or non-numeric values must be encumbered by single quotes (')!
    // all elements must have a value, i.e. booleans also need to be set to TRUE
    def params_list_named  = [
        "css='$report_styles'",
        "report_logo='$report_logo'",
        "workflow_manifest_version='${workflow.manifest.version}'",
        "workflow_scriptid='${workflow.scriptId.substring(0,10)}'",
        params.report_title ? "report_title='$params.report_title'" : "",
        report_abstract ? "report_abstract='$params.report_abstract'" : "",
        meta.single_end ? "flag_single_end=TRUE" : "",
        metadata ? "metadata='$metadata'" : "",
        input_samplesheet ? "input_samplesheet='$input_samplesheet'" : "",
        input_fasta ? "input_fasta='$input_fasta'" : "",
        !input_fasta && !input_samplesheet ? "input_folder='$params.input_folder'" : "",
        mqc_plots ? "mqc_plot='${mqc_plots}/svg/fastqc_per_sequence_quality_scores_plot.svg'" : "",
        cutadapt_summary ?
            params.retain_untrimmed ? "flag_retain_untrimmed=TRUE,cutadapt_summary='$cutadapt_summary'" :
            "cutadapt_summary='$cutadapt_summary'" : "",
        find_truncation_values ? "trunc_qmin=$params.trunc_qmin,trunc_rmin=$params.trunc_rmin" : "",
        "trunclenf='$params.trunclenf'",
        "trunclenr='$params.trunclenr'",
        "max_ee=$params.max_ee",
        dada_qual_stats && meta.single_end ? "dada_qc_f_path='$dada_qual_stats',dada_pp_qc_f_path='$dada_pp_qual_stats'" :
            dada_qual_stats ? "dada_qc_f_path='FW_qual_stats.svg',dada_qc_r_path='RV_qual_stats.svg',dada_pp_qc_f_path='FW_preprocessed_qual_stats.svg',dada_pp_qc_r_path='RV_preprocessed_qual_stats.svg'" : "",
        dada_filtntrim_args ? "dada_filtntrim_args='$dada_filtntrim_args'" : "",
        "dada_sample_inference='$params.sample_inference'",
        dada_err_svgs && meta.run.size() == 1 && meta.single_end ?
            "dada_err_path='$dada_err_svgs',dada_err_run='"+meta.run+"'" :
            dada_err_svgs ? "dada_err_path='"+dada_err_svgs.join(',')+"',dada_err_run='"+meta.run.join(',')+"'" : "",
        dada_asv_table ? "asv_table_path='$dada_asv_table'" : "",
        dada_asv_fa ? "path_asv_fa='$dada_asv_fa'": "",
        dada_tab ? "path_dada2_tab='$dada_tab'" : "",
        dada_stats ? "dada_stats_path='$dada_stats'" : "",
        vsearch_cluster ? "vsearch_cluster='$vsearch_cluster',vsearch_cluster_id='$params.vsearch_cluster_id'" : "",
        params.skip_barrnap ? "" : "path_barrnap_sum='$barrnap_summary'",
        filter_ssu_stats ? "filter_ssu_stats='$filter_ssu_stats'" : "",
        filter_ssu_asv ? "filter_ssu_asv='$filter_ssu_asv',filter_ssu='$params.filter_ssu'" : "",
        filter_len_asv_stats ? "filter_len_asv='$filter_len_asv_stats'" : "",
        filter_len_asv_len_orig ? "filter_len_asv_len_orig='$filter_len_asv_len_orig'" : "",
        params.min_len_asv ? "min_len_asv=$params.min_len_asv" : "min_len_asv=0",
        params.max_len_asv ? "max_len_asv=$params.max_len_asv" : "max_len_asv=0",
        filter_codons_fasta ? "filter_codons_fasta='$filter_codons_fasta',stop_codons='$params.stop_codons'" : "",
        filter_codons_stats ? "filter_codons_stats='$filter_codons_stats'" : "",
        itsx_cutasv_summary ? "itsx_cutasv_summary='$itsx_cutasv_summary',cut_its='$params.cut_its'" : "",
        dada2_tax ? "dada2_taxonomy='$dada2_tax'" : "",
        dada2_tax && !params.dada_ref_tax_custom ? "dada2_ref_tax_title='${params.dada_ref_databases[params.dada_ref_taxonomy]["title"]}',dada2_ref_tax_file='${params.dada_ref_databases[params.dada_ref_taxonomy]["file"]}',dada2_ref_tax_citation='${params.dada_ref_databases[params.dada_ref_taxonomy]["citation"]}'" : "",
        cut_dada_ref_taxonomy ? "cut_dada_ref_taxonomy='$cut_dada_ref_taxonomy'" : "",
        sintax_tax ? "sintax_taxonomy='$sintax_tax',sintax_ref_tax_title='${params.sintax_ref_databases[params.sintax_ref_taxonomy]["title"]}',sintax_ref_tax_file='${params.sintax_ref_databases[params.sintax_ref_taxonomy]["file"]}',sintax_ref_tax_citation='${params.sintax_ref_databases[params.sintax_ref_taxonomy]["citation"]}'" : "",
        kraken2_tax ? "kraken2_taxonomy='$kraken2_tax',kraken2_confidence='$params.kraken2_confidence'" : "",
        kraken2_tax && !params.kraken2_ref_tax_custom ? "kraken2_ref_tax_title='${params.kraken2_ref_databases[params.kraken2_ref_taxonomy]["title"]}',kraken2_ref_tax_file='${params.kraken2_ref_databases[params.kraken2_ref_taxonomy]["file"]}',kraken2_ref_tax_citation='${params.kraken2_ref_databases[params.kraken2_ref_taxonomy]["citation"]}'" : "",
        pplace_tax ? "pplace_taxonomy='$pplace_tax',pplace_heattree='$pplace_heattree'" : "",
        qiime2_tax ? "qiime2_taxonomy='$qiime2_tax'" : "",
        qiime2_tax && params.qiime_ref_taxonomy ? "qiime2_ref_tax_title='${params.qiime_ref_databases[params.qiime_ref_taxonomy]["title"]}',qiime2_ref_tax_file='${params.qiime_ref_databases[params.qiime_ref_taxonomy]["file"]}',qiime2_ref_tax_citation='${params.qiime_ref_databases[params.qiime_ref_taxonomy]["citation"]}'" : "",
        run_qiime2 ? "val_used_taxonomy='$val_used_taxonomy'" : "",
        filter_stats_tsv ? "filter_stats_tsv='$filter_stats_tsv',qiime2_filtertaxa='$qiime2_filtertaxa',exclude_taxa='$params.exclude_taxa',min_frequency='$params.min_frequency',min_samples='$params.min_samples'" : "",
        barplot ? "barplot=TRUE" : "",
        barplot && params.metadata_category_barplot ? "metadata_category_barplot='$params.metadata_category_barplot'" : "",
        abundance_tables ? "abundance_tables=TRUE" : "",
        alpha_rarefaction ? "alpha_rarefaction=TRUE" : "",
        diversity_indices ? "diversity_indices_depth='$diversity_indices'": "",
        diversity_indices_alpha ? "diversity_indices_alpha=TRUE" : "",
        diversity_indices_beta ? "diversity_indices_beta='"+ diversity_indices_beta.join(",") +"'" : "",
        diversity_indices_adonis ? "diversity_indices_adonis='"+ diversity_indices_adonis.join(",") +"',qiime_adonis_formula='$params.qiime_adonis_formula'" : "",
        ancom ? "ancom='"+ ancom.join(",") +"'" : "",
        sbdi ? "sbdi='"+ sbdi.join(",") +"'" : "",
        phyloseq ? "phyloseq='"+ phyloseq.join(",") +"'" : "",
    ]
    // groovy list to R named list string; findAll removes empty entries
    params_list_named_string = params_list_named.findAll().join(',').trim()
    """
    #!/usr/bin/env Rscript
    library(rmarkdown)

    # Work around  https://github.com/rstudio/rmarkdown/issues/1508
    # If the symbolic link is not replaced by a physical file
    # output- and temporary files will be written to the original directory.
    file.copy("./${report_template}", "./template.Rmd", overwrite = TRUE)

    rmarkdown::render("template.Rmd", output_file = "summary_report.html", params = list($params_list_named_string), envir = new.env())

    writeLines(c("\\"${task.process}\\":",
        paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
        paste0("    rmarkdown: ", packageVersion("rmarkdown")),
        paste0("    knitr: ", packageVersion("knitr")) ),
        "versions.yml")
    """
}

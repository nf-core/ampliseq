process SUMMARY_REPORT  {

    label 'process_low'

    container 'docker.io/tillenglert/ampliseq_report:latest'
    /* this is from https://github.com/nf-core/modules/blob/master/modules/nf-core/rmarkdownnotebook/main.nf but doesnt work
    conda "conda-forge::r-base=4.1.0 conda-forge::r-rmarkdown=2.9 conda-forge::r-yaml=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' :
        'biocontainers/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' }"
    */

    input:
    path(report_template)
    path(report_styles)
    path(report_logo)
    path(mqc_plots)
    path(ca_summary)
    val(find_truncation_values)
    path(dada_filtntrim_args)
    path(dada_qual_stats)
    path(dada_pp_qual_stats)
    tuple val(meta), path(dada_err_svgs)
    path(dada_asv_table)
    path(dada_asv_fa)
    path(dada_tab)
    path(dada_stats)
    path(barrnap_summary)
    path(filter_ssu_stats)
    path(filter_ssu_asv)
    path(filter_len_asv_stats)
    path(filter_len_asv_len_orig)
    path(filter_codons_stats)
    path(itsx_cutasv_summary)
    path(dada2_tax)
    tuple val(meta_ref), path(cut_dada_ref_taxonomy) // cutadapt log when params.cut_dada_ref_taxonomy
    path(sintax_tax)
    path(pplace_tax)
    tuple val(meta_pplace), path(pplace_heattree)
    path(qiime2_tax)
    val(run_qiime2)
    val(val_used_taxonomy)
    val(qiime2_filtertaxa) // <ASV count original +1>,<ASV count filtered +2>
    path(filter_stats_tsv)
    path(barplot)
    val(abundance_tables)
    val(alpha_rarefaction)
    path(diversity_indices)
    path(diversity_indices_beta, stageAs: 'beta_diversity/*') // prevent folder name collisons
    path(diversity_indices_adonis, stageAs: 'beta_diversity/adonis/*') // prevent folder name collisons
    path(ancom)
    path(picrust_pathways)


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
        "logo='$report_logo'",
        "workflow_manifest_version='${workflow.manifest.version}'",
        "workflow_scriptid='${workflow.scriptId.substring(0,10)}'",
        meta.single_end ? "flag_single_end=TRUE" : "",
        mqc_plots ? "mqc_plot='${mqc_plots}/svg/mqc_fastqc_per_sequence_quality_scores_plot_1.svg'" : "",
        ca_summary ?
            params.retain_untrimmed ? "flag_retain_untrimmed=TRUE,ca_sum_path='$ca_summary'" :
            "ca_sum_path='$ca_summary'" : "",
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
        params.skip_barrnap ? "" : "path_barrnap_sum='$barrnap_summary'",
        filter_ssu_stats ? "filter_ssu_stats='$filter_ssu_stats',filter_ssu_asv='$filter_ssu_asv',filter_ssu='$params.filter_ssu'" : "",
        filter_len_asv_stats ? "filter_len_asv='$filter_len_asv_stats'" : "",
        filter_len_asv_len_orig ? "filter_len_asv_len_orig='$filter_len_asv_len_orig'" : "",
        params.min_len_asv ? "min_len_asv=$params.min_len_asv" : "min_len_asv=0",
        params.max_len_asv ? "max_len_asv=$params.max_len_asv" : "max_len_asv=0",
        filter_codons_stats ? "filter_codons='$filter_codons_stats',stop_codons='$params.stop_codons'" : "",
        itsx_cutasv_summary ? "itsx_cutasv_summary='$itsx_cutasv_summary',cut_its='$params.cut_its'" : "",
        !dada2_tax ? "" :
            params.dada_ref_tax_custom ? "dada2_taxonomy='$dada2_tax',flag_ref_tax_user=TRUE" :
            "dada2_taxonomy='$dada2_tax',dada2_ref_tax_title='${params.dada_ref_databases[params.dada_ref_taxonomy]["title"]}'",
        cut_dada_ref_taxonomy ? "cut_dada_ref_taxonomy='$cut_dada_ref_taxonomy'" : "",
        sintax_tax ? "sintax_taxonomy='$sintax_tax',sintax_ref_tax_title='${params.sintax_ref_databases[params.sintax_ref_taxonomy]["title"]}'" : "",
        pplace_tax ? "pplace_taxonomy='$pplace_tax',pplace_heattree='$pplace_heattree'" : "",
        qiime2_tax ? "qiime2_taxonomy='$qiime2_tax',qiime2_ref_tax_title='${params.qiime_ref_databases[params.qiime_ref_taxonomy]["title"]}'" : "",
        run_qiime2 ? "val_used_taxonomy='$val_used_taxonomy'" : "",
        filter_stats_tsv ? "filter_stats_tsv='$filter_stats_tsv',qiime2_filtertaxa='$qiime2_filtertaxa',exclude_taxa='$params.exclude_taxa',min_frequency='$params.min_frequency',min_samples='$params.min_samples'" : "",
        barplot ? "barplot=TRUE" : "",
        barplot && params.metadata_category_barplot ? "metadata_category_barplot='$params.metadata_category_barplot'" : "",
        abundance_tables ? "abundance_tables=TRUE" : "",
        alpha_rarefaction ? "alpha_rarefaction=TRUE" : "",
        diversity_indices ? "diversity_indices_depth='$diversity_indices',diversity_indices_beta='"+ diversity_indices_beta.join(",") +"'" : "",
        diversity_indices_adonis ? "diversity_indices_adonis='"+ diversity_indices_adonis.join(",") +"',qiime_adonis_formula='$params.qiime_adonis_formula'" : "",
        ancom ? "ancom='"+ ancom.join(",") +"'" : "",
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

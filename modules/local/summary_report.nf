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
    path(dada2_tax_reference)
    path(dada2_tax)
    path(sintax_tax)
    path(pplace_tax)
    path(qiime2_tax)
    val(run_qiime2)
    path(barplot)
    val(abundance_tables)
    val(alpha_rarefaction)
    val(diversity_indices)
    val(ancom)
    path(picrust_pathways)


    output:
    path "Summary_Report.html"      ,   emit: report
    //path "versions.yml"             ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def single_end = meta.single_end ? "--single_end" : ""
    def fastqc = params.skip_fastqc || params.skip_multiqc ? "--skip_fastqc" : "--mqc_plot ${mqc_plots}/svg/mqc_fastqc_per_sequence_quality_scores_plot_1.svg"
    def cutadapt = params.skip_cutadapt ? "--skip_cutadapt" :
        params.retain_untrimmed ? "--retain_untrimmed --ca_sum_path $ca_summary" :
        "--ca_sum_path $ca_summary"
    // Even when in "dada2_preprocessing.nf" is stated "qc_svg = ch_DADA2_QUALITY1_SVG.collect(sort:true)" the whole path, not only the file name, is used to sort. So FW cannot be guaranteed to be before RV!
    def dada_quality = params.skip_dada_quality ? "--skip_dada_quality" :
        meta.single_end ? "--dada_qc_f_path $dada_qual_stats --dada_pp_qc_f_path $dada_pp_qual_stats" :
        "--dada_qc_f_path 'FW_qual_stats.svg' --dada_qc_r_path 'RV_qual_stats.svg' --dada_pp_qc_f_path 'FW_preprocessed_qual_stats.svg' --dada_pp_qc_r_path 'RV_preprocessed_qual_stats.svg'"
    def find_truncation = find_truncation_values ? "--trunc_qmin $params.trunc_qmin --trunc_rmin $params.trunc_rmin" : ""
    // make comma separated list of error profile path when multiple sequencing runs were performed
    if ( meta.run.size() == 1 && meta.single_end ) {
        dada_err = "--dada_err_path $dada_err_svgs --dada_err_run " + meta.run
    } else {
        dada_err = "--dada_err_path " + dada_err_svgs.join(',') + " --dada_err_run " + meta.run.join(',')
    }
    def barrnap = params.skip_barrnap ? "--skip_barrnap" : "--path_barrnap_sum $barrnap_summary"
        barrnap += filter_ssu_stats ? " --filter_ssu_stats $filter_ssu_stats --filter_ssu_asv $filter_ssu_asv --filter_ssu $params.filter_ssu" : " --filter_ssu none"
    def filter_len_asv = filter_len_asv_stats ? "--filter_len_asv $filter_len_asv_stats --filter_len_asv_len_orig $filter_len_asv_len_orig" : ""
        filter_len_asv += params.min_len_asv ? " --min_len_asv $params.min_len_asv " : " --min_len_asv 0"
        filter_len_asv += params.max_len_asv ? " --max_len_asv $params.max_len_asv" : " --max_len_asv 0"
    def filter_codons = filter_codons_stats ? "--filter_codons $filter_codons_stats --stop_codons $params.stop_codons" : ""
    def itsx_cutasv = itsx_cutasv_summary ? "--itsx_cutasv_summary $itsx_cutasv_summary --cut_its $params.cut_its" : "--cut_its none"
    def dada2_taxonomy = !dada2_tax ? "" :
        params.dada_ref_tax_custom ? "--flag_dada2_taxonomy --dada2_taxonomy $dada2_tax --ref_tax_user" : "--flag_dada2_taxonomy --dada2_taxonomy $dada2_tax --ref_tax_path $dada2_tax_reference"
    def sintax_taxonomy = sintax_tax ? "--flag_sintax_taxonomy --sintax_taxonomy $sintax_tax" : ""
    def pplace_taxonomy = pplace_tax ? "--flag_pplace_taxonomy --pplace_taxonomy $pplace_tax" : ""
    def qiime2_taxonomy = qiime2_tax ? "--flag_qiime2_taxonomy --qiime2_taxonomy $qiime2_tax" : ""
    def qiime2 = run_qiime2 ? "" : "--flag_skip_qiime"
        qiime2 += barplot ? "" : " --flag_skip_barplot"
        qiime2 += abundance_tables ? "" : " --flag_skip_abundance_tables"
        qiime2 += alpha_rarefaction ? "" : " --flag_skip_alpha_rarefaction"
        qiime2 += diversity_indices ? "" : " --flag_skip_diversity_indices"
        qiime2 += ancom ? "" : " --flag_skip_ancom"
    def picrust = picrust_pathways ? "--picrust_pathways $picrust_pathways" : ""
    """
    generate_report.R   --report $report_template \\
                        --output "Summary_Report.html" \\
                        $fastqc \\
                        $cutadapt \\
                        $dada_quality \\
                        --asv_table_path $dada_asv_table \\
                        --path_asv_fa $dada_asv_fa \\
                        --path_dada2_tab $dada_tab \\
                        --dada_stats_path $dada_stats \\
                        --dada_filtntrim_args $dada_filtntrim_args \\
                        $dada_err \\
                        $barrnap \\
                        $single_end \\
                        $find_truncation \\
                        --trunclenf $params.trunclenf \\
                        --trunclenr $params.trunclenr \\
                        --max_ee $params.max_ee \\
                        $filter_len_asv \\
                        $filter_codons \\
                        $itsx_cutasv \\
                        $dada2_taxonomy \\
                        $sintax_taxonomy \\
                        $pplace_taxonomy \\
                        $qiime2_taxonomy \\
                        $qiime2 \\
                        $picrust
    """
    //--pl_results $results_dir \\
    //cat <<-END_VERSIONS > versions.yml
    //"${task.process}":
    //    R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    //END_VERSIONS
}

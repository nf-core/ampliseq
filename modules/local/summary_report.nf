process SUMMARY_REPORT  {

    label 'process_low'

    container 'docker.io/tillenglert/ampliseq_report:latest'
    //conda (params.enable_conda ? "bioconda:r-markdown==0.8--r3.4.1_1" : null)
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/r-markdown:0.8--r3.4.1_1' :
    //   'quay.io/biocontainers/r-markdown:0.8--r3.4.1_1' }"

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
    path(barrnap_gff)
    path(barrnap_summary)
    path(dada2_tax_reference)
    path(dada2_tax)
    path(sintax_tax)
    path(pplace_tax)
    path(qiime2_tax)


    output:
    path "Summary_Report.html"      ,   emit: report
    //path "versions.yml"             ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def single_end = meta.single_end ? "--single_end" : ""
    def fastqc = params.skip_fastqc ? "--skip_fastqc" : "--mqc_plot ${mqc_plots}/svg/mqc_fastqc_per_sequence_quality_scores_plot_1.svg"
    def cutadapt = params.skip_cutadapt ? "--skip_cutadapt" :
        params.retain_untrimmed ? "--retain_untrimmed --ca_sum_path $ca_summary" :
        "--ca_sum_path $ca_summary"
    // Even when in "dada2_preprocessing.nf" is stated "qc_svg = ch_DADA2_QUALITY1_SVG.collect(sort:true)" the whole path, not only the file name, is used to sort. So FW cannot be guaranteed to be before RV!
    def dada_quality = params.skip_dada_quality ? "--skip_dada_quality" :
        meta.single_end ? "--dada_qc_f_path $dada_qual_stats --dada_pp_qc_f_path $dada_pp_qual_stats" :
        "--dada_qc_f_path 'FW_qual_stats.svg' --dada_qc_r_path 'RV_qual_stats.svg' --dada_pp_qc_f_path 'FW_preprocessed_qual_stats.svg' --dada_pp_qc_r_path 'RV_preprocessed_qual_stats.svg'"
    def find_truncation = find_truncation_values ? "--trunc_qmin $params.trunc_qmin --trunc_rmin $params.trunc_rmin" : ""
    def dada_err = meta.single_end ? "--dada_1_err_path $dada_err_svgs" : "--dada_1_err_path ${dada_err_svgs[0]} --dada_2_err_path ${dada_err_svgs[1]}"
    def barrnap = params.skip_barrnap ? "--skip_barrnap" : "--path_rrna_arc ${barrnap_gff[0]} --path_rrna_bac ${barrnap_gff[1]} --path_rrna_euk ${barrnap_gff[2]} --path_rrna_mito ${barrnap_gff[3]} --path_barrnap_sum $barrnap_summary"
    def dada2_taxonomy = dada2_tax ? "--dada2_taxonomy $dada2_tax" : ""
    dada2_taxonomy += params.dada_ref_tax_custom ? " --ref_tax_user" : " --ref_tax_path $dada2_tax_reference"
    def sintax_taxonomy = sintax_tax ? "--sintax_taxonomy $sintax_tax" : ""
    def pplace_taxonomy = pplace_tax ? "--pplace_taxonomy $pplace_tax" : ""
    def qiime2_taxonomy = qiime2_tax ? "--qiime2_taxonomy $qiime2_tax" : ""
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
                        $dada2_taxonomy \\
                        $sintax_taxonomy \\
                        $pplace_taxonomy \\
                        $qiime2_taxonomy
    """
    //--pl_results $results_dir \\
    //cat <<-END_VERSIONS > versions.yml
    //"${task.process}":
    //    R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    //END_VERSIONS
}

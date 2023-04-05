process SUMMARY_REPORT  {

    label 'process_low'

    container 'tillenglert/ampliseq_report:latest'
    //conda (params.enable_conda ? "bioconda:r-markdown==0.8--r3.4.1_1" : null)
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/r-markdown:0.8--r3.4.1_1' :
    //   'quay.io/biocontainers/r-markdown:0.8--r3.4.1_1' }"

    input:
    path(report_template)
    path(report_styles)
    path(mqc_plots)
    path(ca_summary)
    path(dada_filtntrim_args)
    path(dada_fw_qual_stats)
    path(dada_rv_qual_stats)
    path(dada_pp_fw_qual_stats)
    path(dada_pp_rv_qual_stats)
    tuple val(meta), path(dada_err_svgs)
    path(dada_asv_table)
    path(dada_asv_fa)
    path(dada_tab)
    path(dada_stats)
    path(barrnap_gff)
    path(barrnap_summary)
    path(tax_reference)
    path(asv_tax)


    output:
    path "Summary_Report.html"      ,   emit: report
    //path "versions.yml"             ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fastqc = params.skip_fastqc ? "--skip_fastqc" : "--mqc_plot ${mqc_plots}/svg/mqc_fastqc_per_sequence_quality_scores_plot_1.svg"
    def cutadapt = params.skip_cutadapt ? "--skip_cutadapt" : "--ca_sum_path $ca_summary"
    def dada_quality = params.skip_dada_quality ? "--skip_dada_quality" :
        meta.single_end ? "--dada_qc_f_path $dada_fw_qual_stats --dada_pp_qc_f_path $dada_pp_fw_qual_stats" :
        "--dada_qc_f_path $dada_fw_qual_stats --dada_qc_r_path $dada_rv_qual_stats --dada_pp_qc_f_path $dada_pp_fw_qual_stats --dada_pp_qc_r_path $dada_pp_rv_qual_stats"
    def retain_untrimmed = params.retain_untrimmed ? "--retain_untrimmed" : ""
    def single_end = meta.single_end ? "--single_end" : ""
    def dada_err = meta.single_end ? "--dada_1_err_path $dada_err_svgs" : "--dada_1_err_path ${dada_err_svgs[0]} --dada_2_err_path ${dada_err_svgs[1]}"
    def barrnap = params.skip_barrnap ? "--skip_barrnap" : "--path_rrna_arc ${barrnap_gff[0]} --path_rrna_bac ${barrnap_gff[1]} --path_rrna_euk ${barrnap_gff[2]} --path_rrna_mito ${barrnap_gff[3]} --path_barrnap_sum $barrnap_summary"
    def taxonomy = params.skip_taxonomy ? "--skip_taxonomy" :
        params.dada_ref_tax_custom ? "--ref_tax_user --asv_tax_path $asv_tax" : "--ref_tax_path $tax_reference --asv_tax_path $asv_tax"
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
                        $taxonomy \\
                        $retain_untrimmed \\
                        $single_end \\
                        --trunclenf $params.trunclenf \\
                        --trunclenr $params.trunclenr \\
                        --trunc_qmin $params.trunc_qmin
    """
    //--pl_results $results_dir \\
    //cat <<-END_VERSIONS > versions.yml
    //"${task.process}":
    //    R: \$(R --version 2>&1 | sed -n 1p | sed 's/R version //' | sed 's/ (.*//')
    //END_VERSIONS
}

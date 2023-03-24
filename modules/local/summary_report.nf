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


    output:
    path "Summary_Report.html"      ,   emit: report
    //path "versions.yml"             ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def skip_fastqc = params.skip_fastqc ? "--skip_fastqc" : ""
    def cutadapt = params.skip_cutadapt ? "--skip_cutadapt" : "--ca_sum_path $ca_summary"
    def dada_quality = params.skip_dada_quality ? "--skip_dada_quality" : "--dada_qc_f_path $dada_fw_qual_stats --dada_qc_r_path $dada_rv_qual_stats --dada_pp_qc_f_path $dada_pp_fw_qual_stats --dada_pp_qc_r_path $dada_pp_rv_qual_stats"
    def skip_barrnap = params.skip_barrnap ? "--skip_barrnap" : ""
    def retain_untrimmed = params.retain_untrimmed ? "--retain_untrimmed" : ""
    """
    generate_report.R   --report $report_template \\
                        --output "Summary_Report.html" \\
                        --mqc_plot "${mqc_plots}/svg/mqc_fastqc_per_sequence_quality_scores_plot_1.svg" \\
                        --dada_filtntrim_args $dada_filtntrim_args \\
                        --dada_1_err_path ./1_1.err.svg \\
                        --dada_2_err_path ./1_2.err.svg \\
                        $skip_fastqc \\
                        $cutadapt \\
                        $dada_quality \\
                        $skip_barrnap \\
                        $retain_untrimmed \\
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

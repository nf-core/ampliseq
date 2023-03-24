#!/usr/bin/env Rscript

library(rmarkdown)
library(optparse)


option_list = list(
    make_option(c("-r", "--report"), type="character", default=NULL, help="report template file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="ampliseq_report.html", help="output file name", metavar="character"),
    make_option(c("--skip_fastqc"), action="store_true", default=FALSE, help="Trigger to skip fastqc reporting", metavar="logical"),
    make_option(c("--skip_cutadapt"), action="store_true", default=FALSE, help="Trigger to skip cutadapt filtering", metavar="logical"),
    make_option(c("--skip_dada_quality"), action="store_true", default=FALSE, help="Trigger to skip dada2 quality plotting", metavar="logical"),
    make_option(c("--skip_barrnap"), action="store_true", default=FALSE, help="Trigger to skip barrnap ASV filtering", metavar="logical"),
    make_option(c("--retain_untrimmed"), action="store_true", default=FALSE, help="Flag to retain the untrimmed sequences", metavar="logical"),
    make_option(c("--trunclenf"), type="numeric", default=-1, help="Flag to define truncation in forward strand", metavar="numeric"),
    make_option(c("--trunclenr"), type="numeric", default=-1, help="Flag to define truncation in reverse strand", metavar="numeric"),
    make_option(c("--trunc_qmin"), type="numeric", default=-1, help="Flag to define truncation via quality measure. Set to -1 if trunclen were given.", metavar="numeric"),
    make_option(c("--mqc_plot"), type="character", default=NULL, help="MultiQC plot per sequence quality", metavar="character"),
    make_option(c("--ca_sum_path"), type="character", default=NULL, help="cutadapt summary table", metavar="character"),
    make_option(c("--dada_filtntrim_args"), type="character", default=NULL, help="DADA2 arguments for filter and trim process", metavar="character"),
    make_option(c("--dada_qc_f_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--dada_qc_r_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--dada_pp_qc_f_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--dada_pp_qc_r_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--dada_1_err_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--dada_2_err_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--asv_table_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--path_asv_fa"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--path_dada2_tab"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--dada_stats_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--path_rrna_arc"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--path_rrna_bac"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--path_rrna_euk"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--path_rrna_mito"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--ref_tax_path"), type="character", default=NULL, help="MultiQC plots", metavar="character"),
    make_option(c("--asv_tax_path"), type="character", default=NULL, help="MultiQC plots", metavar="character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

rmarkdown::render(opt$report, output_file = opt$output,
                    params = list(
                        flag_skip_fastqc = opt$skip_fastqc,
                        flag_skip_cutadapt = opt$skip_cutadapt,
                        flag_skip_dada_quality = opt$skip_dada_quality,
                        flag_skip_barrnap = opt$skip_barrnap,
                        flag_retain_untrimmed = opt$retain_untrimmed,
                        flag_trunclenf = opt$trunclenf,
                        flag_trunclenr = opt$trunclenr,
                        flag_trunc_qmin = opt$trunc_qmin,
                        mqc_plot = opt$mqc_plot,
                        ca_sum_path = opt$ca_sum_path,
                        dada_filtntrim_args = opt$dada_filtntrim_args,
                        dada_qc_f_path = opt$dada_qc_f_path,
                        dada_qc_r_path = opt$dada_qc_r_path,
                        dada_pp_qc_f_path = opt$dada_pp_qc_f_path,
                        dada_pp_qc_r_path = opt$dada_pp_qc_r_path,
                        dada_1_err_path = opt$dada_1_err_path,
                        dada_2_err_path = opt$dada_2_err_path,
                        asv_table_path = opt$asv_table_path,
                        path_asv_fa = opt$path_asv_fa,
                        path_dada2_tab = opt$path_dada2_tab,
                        dada_stats_path = opt$dada_stats_path,
                        path_rrna_arc = opt$path_rrna_arc,
                        path_rrna_bac = opt$path_rrna_bac,
                        path_rrna_euk = opt$path_rrna_euk,
                        path_rrna_mito = opt$path_rrna_mito,
                        ref_tax_path = opt$ref_tax_path,
                        asv_tax_path = opt$asv_tax_path))

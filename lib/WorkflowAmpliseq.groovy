//
// This file holds several functions specific to the workflow/ampliseq.nf in the nf-core/ampliseq pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowAmpliseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if ( !params.input && !params.input_fasta && !params.input_folder ) {
            Nextflow.error("Missing input declaration: One of `--input`, `--input_fasta`, `--input_folder` is required.")
        }

        if ( !params.input_fasta && (!params.FW_primer || !params.RV_primer) && !params.skip_cutadapt ) {
            Nextflow.error("Incompatible parameters: `--FW_primer` and `--RV_primer` are required for primer trimming. If primer trimming is not needed, use `--skip_cutadapt`.")
        }

        if ( params.pacbio || params.iontorrent || params.single_end ) {
            if (params.trunclenr) { log.warn "Unused parameter: `--trunclenr` is ignored because the data is single end." }
        } else if (params.trunclenf && !params.trunclenr) {
            Nextflow.error("Invalid command: `--trunclenf` is set, but `--trunclenr` is not. Either both parameters `--trunclenf` and `--trunclenr` must be set or none.")
        } else if (!params.trunclenf && params.trunclenr) {
            Nextflow.error("Invalid command: `--trunclenr` is set, but `--trunclenf` is not. Either both parameters `--trunclenf` and `--trunclenr` must be set or none.")
        }

        if (!["pooled", "independent", "pseudo"].contains(params.sample_inference)) {
            def error_string = "Please set `--sample_inference` to one of the following:\n" +
                "\t-\"independent\" (lowest sensitivity and lowest resources),\n" +
                "\t-\"pseudo\" (balance between required resources and sensitivity),\n" +
                "\t-\"pooled\" (highest sensitivity and resources)."
            Nextflow.error(error_string)
        }

        if (params.double_primer && params.retain_untrimmed) {
            Nextflow.error("Incompatible parameters `--double_primer` and `--retain_untrimmed` cannot be set at the same time.")
        }

        if ( params.min_len_asv && params.max_len_asv && (params.min_len_asv > params.max_len_asv) ) {
            Nextflow.error("Incompatible parameters: `--min_len_asv` may not be greater than `--max_len_asv`.")
        }

        if ( params.skip_dada_quality && (params.trunclenf == null || params.trunclenr == null) ) {
            Nextflow.error("Incompatible parameters: `--skip_dada_quality` may not be used without setting `--trunclenf` and `--trunclenr`.")
        }

        if (params.tax_agglom_min > params.tax_agglom_max) {
            Nextflow.error("Incompatible parameters: `--tax_agglom_min` may not be greater than `--tax_agglom_max`.")
        }

        if (!params.dada_ref_tax_custom && params.dada_ref_tax_custom_sp) {
            Nextflow.error("Incompatible parameters: `--dada_ref_tax_custom_sp` requires `--dada_ref_tax_custom`.")
        }

        if (params.dada_ref_tax_custom && !params.dada_ref_tax_custom_sp && !params.skip_dada_addspecies) {
            Nextflow.error("Incompatible parameters: Either `--skip_dada_addspecies` or `--dada_ref_tax_custom_sp` is additionally required to `--dada_ref_tax_custom`.")
        }

        if (params.pplace_tree) {
            if (!params.pplace_aln) {
                Nextflow.error("Missing parameter: Phylogenetic placement requires in addition to `--pplace_tree` also `--pplace_aln`.")
            }
            if (!params.pplace_model) {
                Nextflow.error("Missing parameter: Phylogenetic placement requires in addition to `--pplace_tree` also `--pplace_model`.")
            }
        }

        if (params.dada_assign_taxlevels && params.sbdiexport && !params.sintax_ref_taxonomy) {
            Nextflow.error("Incompatible parameters: `--sbdiexport` expects specific taxonomics ranks (default) and therefore excludes modifying those using `--dada_assign_taxlevels`.")
        }

        if (params.skip_taxonomy && params.sbdiexport) {
            Nextflow.error("Incompatible parameters: `--sbdiexport` expects taxa annotation and therefore excludes `--skip_taxonomy`.")
        }

        if (params.skip_dada_taxonomy && params.sbdiexport) {
            if (!params.sintax_ref_taxonomy && (params.skip_qiime || (!params.qiime_ref_taxonomy && !params.qiime_ref_tax_custom))) {
                Nextflow.error("Incompatible parameters: `--sbdiexport` expects taxa annotation and therefore annotation with either DADA2, SINTAX, or QIIME2 is needed.")
            }
        }

        if ( (!params.FW_primer || !params.RV_primer) && (params.qiime_ref_taxonomy || params.qiime_ref_tax_custom) && !params.skip_qiime && !params.skip_taxonomy ) {
            Nextflow.error("Incompatible parameters: `--FW_primer` and `--RV_primer` are required for cutting the QIIME2 reference database to the amplicon sequences. Please specify primers or do not use `--qiime_ref_taxonomy`.")
        }

        if ( (!params.FW_primer || !params.RV_primer) && params.cut_dada_ref_taxonomy && !params.skip_taxonomy ) {
            Nextflow.error("Incompatible parameters: `--FW_primer` and `--RV_primer` are required for cutting the DADA2 reference database to the amplicon sequences. Please specify primers or do not use `--cut_dada_ref_taxonomy`.")
        }

        if ((params.qiime_ref_taxonomy || params.qiime_ref_tax_custom) && params.classifier) {
            Nextflow.error("Incompatible parameters: `--qiime_ref_taxonomy` and `--qiime_ref_tax_custom` will produce a classifier but `--classifier` points to a precomputed classifier, therefore, only use one of those.")
        }

        if (params.kraken2_ref_tax_custom && !params.kraken2_assign_taxlevels ) {
            Nextflow.error("Missing parameter: Taxonomic classification with a user provided database via `--kraken2_ref_tax_custom` requires `--kraken2_assign_taxlevels`")
        }

        if (params.filter_ssu && params.skip_barrnap) {
            Nextflow.error("Incompatible parameters: `--filter_ssu` cannot be used with `--skip_barrnap` because filtering for SSU's depends on barrnap.")
        }

        String[] sbdi_compatible_databases = ["coidb","coidb=221216","gtdb","gtdb=R08-RS214","gtdb=R07-RS207","gtdb=R06-RS202","gtdb=R05-RS95","midori2-co1","midori2-co1=gb250","pr2","pr2=5.0.0","pr2=4.14.0","pr2=4.13.0","rdp","rdp=18","sbdi-gtdb","sbdi-gtdb=R07-RS207-1","silva","silva=138","silva=132","unite-fungi","unite-fungi=9.0","unite-fungi=8.3","unite-fungi=8.2","unite-alleuk","unite-alleuk=9.0","unite-alleuk=8.3","unite-alleuk=8.2"]
        if (params.sbdiexport){
            if (params.sintax_ref_taxonomy ) {
                if (!Arrays.stream(sbdi_compatible_databases).anyMatch(entry -> params.sintax_ref_taxonomy.toString().equals(entry)) ) {
                    Nextflow.error("Incompatible parameters: `--sbdiexport` does not work with the chosen database of `--sintax_ref_taxonomy` because the expected taxonomic levels do not match.")
                }
            } else if (!Arrays.stream(sbdi_compatible_databases).anyMatch(entry -> params.dada_ref_taxonomy.toString().equals(entry)) ) {
                Nextflow.error("Incompatible parameters: `--sbdiexport` does not work with the chosen database of `--dada_ref_taxonomy` because the expected taxonomic levels do not match.")
            }
        }

        if (params.addsh && !params.dada_ref_databases[params.dada_ref_taxonomy]["shfile"]) {
            def validDBs = ""
            for (db in params.dada_ref_databases.keySet()) {
                if (params.dada_ref_databases[db]["shfile"]) {
                    validDBs += " " + db
                }
            }
            Nextflow.error("UNITE species hypothesis information is not available for the selected reference database, please use the option `--dada_ref_taxonomy` to select an appropriate database. Currently, the option `--addsh` can only be used together with the following UNITE reference databases:\n" + validDBs + ".")
        }

        if (params.addsh && params.cut_its == "none") {
            log.warn "Adding UNITE species hypothesis (SH) assignments is only feasible for ITS sequences. Please use option `--cut_its` to find ITS regions in the ASV sequences, unless the given sequences are already cut to the ITS region.\n"
        }

        // Error message for incompatible combination of --orf_start and --orf_end
        if ( params.orf_end && ( ( ( params.orf_end + 1 ) - params.orf_start ) % 3 != 0 ) ) {
            Nextflow.error("Incompatible parameters: The difference of  `--orf_end` and `--orf_start` must be a multiple of 3.")
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        // TODO nf-core: Optionally add in-text citation tools to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def citation_text = [
                "Tools used in the workflow included:",
                "FastQC (Andrews 2010),",
                "MultiQC (Ewels et al. 2016)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        // TODO Optionally add bibliographic entries to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        // TODO Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
        //meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        //meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }
}

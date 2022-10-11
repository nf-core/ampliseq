//
// This file holds several functions specific to the workflow/ampliseq.nf in the nf-core/ampliseq pipeline
//

import groovy.text.SimpleTemplateEngine

class WorkflowAmpliseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (params.enable_conda) { log.warn "Conda is enabled (`--enable_conda`), any steps involving QIIME2 are not available. Use a container engine instead of conda to enable all software." }

        if ( params.pacbio || params.iontorrent || params.single_end ) {
            if (params.trunclenr) { log.warn "Unused parameter: `--trunclenr` is ignored because the data is single end." }
        } else if (params.trunclenf && !params.trunclenr) {
            log.error "Invalid command: `--trunclenf` is set, but `--trunclenr` is not. Either both parameters `--trunclenf` and `--trunclenr` must be set or none."
            System.exit(1)
        } else if (!params.trunclenf && params.trunclenr) {
            log.error "Invalid command: `--trunclenr` is set, but `--trunclenf` is not. Either both parameters `--trunclenf` and `--trunclenr` must be set or none."
            System.exit(1)
        }

        if (!["pooled", "independent", "pseudo"].contains(params.sample_inference)) {
            log.error "Please set `--sample_inference` to one of the following:\n" +
                "\t-\"independent\" (lowest sensitivity and lowest resources),\n" +
                "\t-\"pseudo\" (balance between required resources and sensitivity),\n" +
                "\t-\"pooled\" (highest sensitivity and resources)."
            System.exit(1)
        }

        if (params.double_primer && params.retain_untrimmed) {
            log.error "Incompatible parameters `--double_primer` and `--retain_untrimmed` cannot be set at the same time."
            System.exit(1)
        }

        if ( params.min_len_asv && params.max_len_asv && (params.min_len_asv > params.max_len_asv) ) {
            log.error "Incompatible parameters: `--min_len_asv` may not be greater than `--max_len_asv`."
            System.exit(1)
        }

        if ( params.skip_dada_quality && (params.trunclenf == null || params.trunclenr == null) ) {
            log.error "Incompatible parameters: `--skip_dada_quality` may not be used without setting `--trunclenf` and `--trunclenr`."
            System.exit(1)
        }

        if (params.dada_tax_agglom_min > params.dada_tax_agglom_max) {
            log.error "Incompatible parameters: `--dada_tax_agglom_min` may not be greater than `--dada_tax_agglom_max`."
            System.exit(1)
        }

        if (params.qiime_tax_agglom_min > params.qiime_tax_agglom_max) {
            log.error "Incompatible parameters: `--qiime_tax_agglom_min` may not be greater than `--qiime_tax_agglom_max`."
            System.exit(1)
        }

        if (!params.dada_ref_tax_custom && params.dada_ref_tax_custom_sp) {
            log.error "Incompatible parameters: `--dada_ref_tax_custom_sp` requires `--dada_ref_tax_custom`."
            System.exit(1)
        }

        if (params.dada_ref_tax_custom && !params.dada_ref_tax_custom_sp && !params.skip_dada_addspecies) {
            log.error "Incompatible parameters: Either `--skip_dada_addspecies` or `--dada_ref_tax_custom_sp` is additionally required to `--dada_ref_tax_custom`."
            System.exit(1)
        }

        if (params.dada_assign_taxlevels && params.sbdiexport) {
            log.error "Incompatible parameters: `--sbdiexport` expects specific taxonomics ranks (default) and therefore excludes modifying those using `--dada_assign_taxlevels`."
            System.exit(1)
        }

        if (params.skip_dada_addspecies && params.sbdiexport) {
            log.error "Incompatible parameters: `--sbdiexport` expects species annotation and therefore excludes `--skip_dada_addspecies`."
            System.exit(1)
        }

        if (params.skip_taxonomy && params.sbdiexport) {
            log.error "Incompatible parameters: `--sbdiexport` expects taxa annotation and therefore excludes `--skip_taxonomy`."
            System.exit(1)
        }

        if ( (!params.FW_primer || !params.RV_primer) && params.qiime_ref_taxonomy && !params.skip_qiime && !params.skip_taxonomy ) {
            log.error "Incompatible parameters: `--FW_primer` and `--RV_primer` are required for cutting the QIIME2 reference database to the amplicon sequences. Please specify primers or do not use `--qiime_ref_taxonomy`."
            System.exit(1)
        }

        if ( (!params.FW_primer || !params.RV_primer) && params.cut_dada_ref_taxonomy && !params.skip_taxonomy ) {
            log.error "Incompatible parameters: `--FW_primer` and `--RV_primer` are required for cutting the DADA2 reference database to the amplicon sequences. Please specify primers or do not use `--cut_dada_ref_taxonomy`."
            System.exit(1)
        }

        if (params.qiime_ref_taxonomy && params.classifier) {
            log.error "Incompatible parameters: `--qiime_ref_taxonomy` will produce a classifier but `--classifier` points to a precomputed classifier, therefore, only use one of those."
            System.exit(1)
        }

        if (params.filter_ssu && params.skip_barrnap) {
            log.error "Incompatible parameters: `--filter_ssu` cannot be used with `--skip_barrnap` because filtering for SSU's depends on barrnap."
            System.exit(1)
        }

        String[] sbdi_incompatible_databases = ["midori2-co1=gb250","midori2-co1","rdp=18","rdp","sbdi-gtdb","sbdi-gtdb=R06-RS202-3","sbdi-gtdb=R06-RS202-1","silva=132","unite-fungi","unite-fungi=8.3","unite-fungi=8.2","unite-alleuk","unite-alleuk=8.3","unite-alleuk=8.2"]
        if ( params.sbdiexport && Arrays.stream(sbdi_incompatible_databases).anyMatch(entry -> params.dada_ref_taxonomy.toString().equals(entry)) ) {
            log.error "Incompatible parameters: `--sbdiexport` does not work with the chosen databse of `--dada_ref_taxonomy`, because the expected taxonomic levels do not match."
            System.exit(1)
        }

        if (params.addsh && !params.dada_ref_databases[params.dada_ref_taxonomy]["shfile"]) {
            def validDBs = ""
            for (db in params.dada_ref_databases.keySet()) {
                if (params.dada_ref_databases[db]["shfile"]) {
                    validDBs += " " + db
                }
            }
            log.error "UNITE species hypothesis information is not available for the selected reference database, please use the option `--dada_ref_taxonomy` to select an appropriate database. Currently, the option `--addsh` can only be used together with the following UNITE reference databases:\n" + validDBs + "."
            System.exit(1)
        }

        if (params.addsh && params.cut_its == "none") {
            log.warn "Adding UNITE species hypothesis (SH) assignments is only feasible for ITS sequences. Please use option `--cut_its` to find ITS regions in the ASV sequences, unless the given sequences are already cut to the ITS region.\n"
        }
    }

    //
    // Check string (String s) ends with one entry of an array of strings ("String[] extn")
    //
    public static boolean checkIfFileHasExtension(String s, String[] extn) {
        return Arrays.stream(extn).anyMatch(entry -> s.endsWith(entry));
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

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }
}

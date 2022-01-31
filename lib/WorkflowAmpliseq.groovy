//
// This file holds several functions specific to the workflow/ampliseq.nf in the nf-core/ampliseq pipeline
//

class WorkflowAmpliseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        if (params.enable_conda) { log.warn "Conda is enabled (`--enable_conda`), any steps involving QIIME2 are not available. Use a container engine instead of conda to enable all software." }

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

        if (params.dada_tax_agglom_min > params.dada_tax_agglom_max) {
            log.error "Incompatible parameters: `--dada_tax_agglom_min` may not be greater than `--dada_tax_agglom_max`."
            System.exit(1)
        }

        if (params.qiime_tax_agglom_min > params.qiime_tax_agglom_max) {
            log.error "Incompatible parameters: `--qiime_tax_agglom_min` may not be greater than `--qiime_tax_agglom_max`."
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
}

// Import generic module functions
//include { initOptions; saveFiles } from './functions'

process FORMAT_TAXONOMY {
    label 'process_low'
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    input:
    path(database)
    
    output:
    path( "*assignTaxonomy.fna" ), emit: assigntax
    path( "*addSpecies.fna"), emit: addspecies

    script:
    """
    ${params.genomes[params.dada_ref_taxonomy]["fmtscript"]}
    """
}

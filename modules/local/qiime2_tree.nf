// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options    = initOptions(params.options)

process QIIME2_TREE {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? { exit 1 "QIIME2 has no conda package" } : null)
    container "quay.io/qiime2/core:2021.2"

    input:
    path(repseq)

    output:
    path("rooted-tree.qza"), emit: qza
    path("tree.nwk")       , emit: nwk
    path "*.version.txt"   , emit: version

    script:
    def software     = getSoftwareName(task.process)
    """
    export XDG_CONFIG_HOME="\${PWD}/HOME"

    qiime alignment mafft \
        --i-sequences ${repseq} \
        --o-alignment aligned-rep-seqs.qza \
        --p-n-threads ${task.cpus}
    qiime alignment mask \
        --i-alignment aligned-rep-seqs.qza \
        --o-masked-alignment masked-aligned-rep-seqs.qza
    qiime phylogeny fasttree \
        --i-alignment masked-aligned-rep-seqs.qza \
        --p-n-threads ${task.cpus} \
        --o-tree unrooted-tree.qza
    qiime phylogeny midpoint-root \
        --i-tree unrooted-tree.qza \
        --o-rooted-tree rooted-tree.qza
    qiime tools export --input-path rooted-tree.qza  \
        --output-path phylogenetic_tree
    cp phylogenetic_tree/tree.nwk .

    echo \$(qiime --version | sed -e "s/q2cli version //g" | tr -d '`' | sed -e "s/Run qiime info for more version details.//g") > ${software}.version.txt
    """
}
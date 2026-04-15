process ITSXRUST_CUTASV {
    label 'process_medium'

    conda "bioconda::itsxrust=0.2.2 bioconda::hmmer=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/itsxrust:0.2.2--hdd79491_1' :
        'quay.io/biocontainers/itsxrust:0.2.2--hdd79491_1' }"

    input:
    path fasta
    val outfile

    output:
    path outfile         , emit: fasta
    path "ASV_ITS_seqs.summary.txt", emit: summary
    path "ASV_ITS_seqs.*fasta", emit: fastas
    path "versions.yml"  , emit: versions
    path "*.args.txt"    , emit: args

    script:
    def args = task.ext.args ?: ''
    def hmm_path = task.ext.hmm ?: "\${CONDA_PREFIX:-/usr/local}/share/itsxrust/hmm/F.hmm"
    """
    # Run ITSxRust with --region all so we always get full + ITS1 + ITS2 outputs
    # (mirrors ITSx behavior which always writes all region files)
    itsxrust extract \\
        --input $fasta \\
        --hmm $hmm_path \\
        $args \\
        --hmmer-cpu $task.cpus \\
        --output ASV_ITS_seqs \\
        --region all \\
        --output-format fasta \\
        --qc-json itsxrust_qc.json

    # ITSxRust outputs: ASV_ITS_seqs.its1.fasta, ASV_ITS_seqs.its2.fasta, ASV_ITS_seqs.full.fasta
    # The pipeline expects: ASV_ITS_seqs.ITS1.fasta, ASV_ITS_seqs.ITS2.fasta, ASV_ITS_seqs.full.fasta
    # Rename ITS1/ITS2 to match ITSx naming convention (uppercase)
    if [ -f ASV_ITS_seqs.its1.fasta ]; then
        cp ASV_ITS_seqs.its1.fasta ASV_ITS_seqs.ITS1.tmp && mv ASV_ITS_seqs.ITS1.tmp ASV_ITS_seqs.ITS1.fasta && rm -f ASV_ITS_seqs.its1.fasta
    fi
    if [ -f ASV_ITS_seqs.its2.fasta ]; then
        cp ASV_ITS_seqs.its2.fasta ASV_ITS_seqs.ITS2.tmp && mv ASV_ITS_seqs.ITS2.tmp ASV_ITS_seqs.ITS2.fasta && rm -f ASV_ITS_seqs.its2.fasta
    fi

    # Handle partial naming: if its_partial is used, the workflow expects
    # filenames like ASV_ITS_seqs.full_and_partial.fasta
    # This is handled by the outfile val — if outfile contains "full_and_partial",
    # we rename accordingly (ITSxRust's partial-chain fallback includes partials by default)
    if echo "$outfile" | grep -q "full_and_partial"; then
        if [ -f ASV_ITS_seqs.full.fasta ]; then
            cp ASV_ITS_seqs.full.fasta "$outfile"
        fi
        if [ -f ASV_ITS_seqs.ITS1.fasta ]; then
            cp ASV_ITS_seqs.ITS1.fasta ASV_ITS_seqs.ITS1.full_and_partial.fasta
        fi
        if [ -f ASV_ITS_seqs.ITS2.fasta ]; then
            cp ASV_ITS_seqs.ITS2.fasta ASV_ITS_seqs.ITS2.full_and_partial.fasta
        fi
    fi

    # Strip ITSxRust coordinate annotations from FASTA headers (for now)
    # e.g. ">seq1|full:47-433" becomes ">seq1"
    for f in ASV_ITS_seqs.*.fasta; do
        if [ -f "\$f" ]; then
            sed -i 's/|[a-z]*[0-9]*:[0-9]*-[0-9]*//' "\$f"
        fi
    done

    # Generate ITSx-compatible summary from ITSxRust QC JSON
    if [ -f itsxrust_qc.json ]; then
        total=\$(grep -o '"total_reads":[0-9]*' itsxrust_qc.json | grep -o '[0-9]*')
        kept=\$(grep -o '"kept":[0-9]*' itsxrust_qc.json | grep -o '[0-9]*')
        skipped=\$(grep -o '"skipped":[0-9]*' itsxrust_qc.json | grep -o '[0-9]*')
        echo "ITSxRust extraction summary" > ASV_ITS_seqs.summary.txt
        echo "Number of sequences in input file: \${total:-0}" >> ASV_ITS_seqs.summary.txt
        echo "Sequences detected as ITS by ITSx: \${kept:-0}" >> ASV_ITS_seqs.summary.txt
        echo "Number of sequences skipped: \${skipped:-0}" >> ASV_ITS_seqs.summary.txt
    else
        echo "ITSxRust: QC JSON not available" > ASV_ITS_seqs.summary.txt
    fi

    # Validate that the expected output file exists and is non-empty
    if [ ! -s "$outfile" ]; then
        echo "ERROR: No ITS regions found by ITSxRust. You might want to modify --cut_its, --its_partial, or --its_extractor options" >&2
        exit 1
    fi

    echo -e "ITSxRust\\t$args" > ITSxRust.args.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ITSxRust: \$( itsxrust --version 2>&1 | sed 's/itsxrust //' )
    END_VERSIONS
    """
}

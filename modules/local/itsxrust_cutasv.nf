process ITSXRUST_CUTASV {
    label 'process_medium'

    conda "bioconda::itsxrust=0.2.2 bioconda::hmmer=3.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/ayobi/itsxrust:0.2.2' :
        'ghcr.io/ayobi/itsxrust:0.2.2' }"

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
    def hmm_path = task.ext.hmm ?: '/usr/local/share/itsxrust/hmm/F.hmm'

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

    # Generate ITSx-compatible summary from ITSxRust QC JSON
    # The summary report parses lines like "Number of sequences in input: X"
    python3 <<'EOF'
import json, sys

try:
    with open("itsxrust_qc.json") as f:
        qc = json.load(f)
except (FileNotFoundError, json.JSONDecodeError):
    # Fallback: write minimal summary
    with open("ASV_ITS_seqs.summary.txt", "w") as out:
        out.write("ITSxRust: QC JSON not available\\n")
    sys.exit(0)

total = qc.get("total_reads", 0)
kept = qc.get("kept", 0)
skipped = qc.get("skipped", 0)

with open("ASV_ITS_seqs.summary.txt", "w") as out:
    out.write(f"ITSxRust extraction summary\\n")
    out.write(f"Number of sequences in input: {total}\\n")
    out.write(f"Number of sequences extracted: {kept}\\n")
    out.write(f"Number of sequences skipped: {skipped}\\n")
    # Confidence breakdown if available
    conf = qc.get("confidence", {})
    if conf:
        out.write(f"  Confident: {conf.get('confident', 0)}\\n")
        out.write(f"  Ambiguous: {conf.get('ambiguous', 0)}\\n")
        out.write(f"  Partial: {conf.get('partial', 0)}\\n")
    # Skip reasons if available
    reasons = qc.get("skip_reasons", {})
    if reasons:
        out.write("Skip reasons:\\n")
        for reason, count in reasons.items():
            out.write(f"  {reason}: {count}\\n")
EOF

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

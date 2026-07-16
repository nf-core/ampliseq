#!/bin/sh

# Reformat GloSED references for assignTaxonomy and addSpecies.
# Inputs (in current directory by default):
#   - GloSED__OTU_sequences.fasta.gz
#   - GloSED__Taxonomy.tsv.zip
# Outputs:
#   - assignTaxonomy.fna
#   - addSpecies.fna

set -eu

stream_seq() {
    case "$1" in
        *.zip) unzip -p "$1" ;;
        *) gzip -dc "$1" ;;
    esac
}

SEQ_GZ="${1:-GloSED__OTU_sequences.fasta.gz}"
TAX_ZIP="${2:-GloSED__Taxonomy.tsv.zip}"

if [ -z "${1:-}" ] && [ -z "${2:-}" ] && [ ! -f "$SEQ_GZ" ] && [ ! -f "$TAX_ZIP" ] && [ -f "GloSED__OTU_sequences_first1k.fasta.zip" ] && [ -f "GloSED_Taxonomy_first1k.tsv.zip" ]; then
    SEQ_GZ="GloSED__OTU_sequences_first1k.fasta.zip"
    TAX_ZIP="GloSED_Taxonomy_first1k.tsv.zip"
fi

[ -f "$SEQ_GZ" ] || { echo "ERROR: missing $SEQ_GZ" >&2; exit 1; }
[ -f "$TAX_ZIP" ] || { echo "ERROR: missing $TAX_ZIP" >&2; exit 1; }

TMPDIR_LOCAL="$(mktemp -d glosed_reformat.XXXXXX)"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT INT TERM

TAX_TSV="${TMPDIR_LOCAL}/glosed_taxonomy.tsv"
META_TSV="${TMPDIR_LOCAL}/glosed_meta.tsv"

unzip -p "$TAX_ZIP" > "$TAX_TSV"

awk -F '\t' 'BEGIN { OFS="\t" }
NR==1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "OTU") otu = i
        else if ($i == "Kingdom") kingdom = i
        else if ($i == "Phylum") phylum = i
        else if ($i == "Class") classcol = i
        else if ($i == "Order") ordercol = i
        else if ($i == "Family") family = i
        else if ($i == "Genus") genus = i
        else if ($i == "Species") species = i
    }
    next
}
{
    id = $otu
    k = $kingdom; p = $phylum; c = $classcol; o = $ordercol; f = $family; g = $genus; s = $species
    gsub(/ /, "_", k); gsub(/ /, "_", p); gsub(/ /, "_", c); gsub(/ /, "_", o); gsub(/ /, "_", f); gsub(/ /, "_", g); gsub(/ /, "_", s)

    tax = ""
    if (k != "." && k != "") {
        tax = k ";"
        if (p != "." && p != "") {
            tax = tax p ";"
            if (c != "." && c != "") {
                tax = tax c ";"
                if (o != "." && o != "") {
                    tax = tax o ";"
                    if (f != "." && f != "") {
                        tax = tax f ";"
                        if (g != "." && g != "") {
                            tax = tax g ";"
                            if (s != "." && s != "") {
                                tax = tax s ";"
                            }
                        }
                    }
                }
            }
        }
    }
    print id, tax, g, s
}
' "$TAX_TSV" > "$META_TSV"

stream_seq "$SEQ_GZ" | awk -F '\t' 'NR==FNR { tax[$1] = $2; next }
/^>/ { id = substr($0,2); if (id in tax && tax[id] != "") print ">" tax[id]; else print ">" id; next }
{ print }
' "$META_TSV" - > assignTaxonomy.fna

stream_seq "$SEQ_GZ" | awk -F '\t' 'NR==FNR { if ($3 != "." && $3 != "" && $3 != "NA" && $4 != "." && $4 != "" && $4 != "NA") addsp[$1] = $3 " " $4; next }
/^>/ { id = substr($0,2); keep = (id in addsp); if (keep) print ">" id " " addsp[id]; next }
{ if (keep) print }
' "$META_TSV" - > addSpecies.fna

echo "Created files: assignTaxonomy.fna addSpecies.fna"

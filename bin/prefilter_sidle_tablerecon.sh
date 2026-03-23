#!/usr/bin/env sh

set -eu

if [ "$#" -lt 2 ]; then
    echo "Usage: prefilter_sidle_tablerecon.sh <min_counts> <regional_table.qza...>" >&2
    exit 1
fi

min_counts="$1"
shift

# Export original regional tables as TSVs in the task work dir.
for region_table in "$@"; do
    region_table_base="$(basename "$region_table" .qza)"
    original_table_exported_folder="${region_table_base}_exported"
    original_table_tsv="${region_table_base}_feature-table.tsv"

    qiime tools export \
        --input-path "$region_table" \
        --output-path "${original_table_exported_folder}"
    biom convert \
        -i "${original_table_exported_folder}/feature-table.biom" \
        -o "${original_table_tsv}" \
        --to-tsv
done

# Determine total counts per sample across all regional tables, to be used as input for reconstruction.
{
    printf "sample-id\ttotal_count\n"
    awk '
    BEGIN { FS=OFS="\t" }

    $1 == "#OTU ID" {
        delete cols
        for (i=2; i<=NF; i++) cols[i] = $i
        next
    }

    /^#/ { next }

    {
        for (i=2; i<=NF; i++) sum[cols[i]] += $i
    }

    END {
        for (s in sum) {
            if (s != "") print s, sum[s]
        }
    }
    ' *_feature-table.tsv | sort -k1,1
} > total_counts.tsv

# Extract all samples with total counts >= min_counts into kept_samples.tsv.
awk -v min_counts="${min_counts}" '
BEGIN { FS=OFS="\t" }
NR==1 { print "sample-id"; next }
$2 >= min_counts { print $1 }
' total_counts.tsv > kept_samples.tsv

# Filter the regional tables to keep only samples meeting min_counts.
for region_table in "$@"; do
    region_table_base="$(basename "$region_table" .qza)"
    qiime feature-table filter-samples \
        --i-table "$region_table" \
        --m-metadata-file kept_samples.tsv \
        --o-filtered-table "${region_table_base}.filtered.qza"
done

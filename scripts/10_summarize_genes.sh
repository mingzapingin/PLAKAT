#!/usr/bin/env bash
set -euo pipefail
# Usage: 10_summarize_genes.sh <marker_annotations.tsv> <summary_out.tsv>

input="$1"   # results/core/<POS>_marker_annotations.tsv
out="$2"     # results/core/<POS>_marker_summary.tsv

echo -e "marker\tfeature_type\tgene\tlocus_tag" > "$out"

# Notes:
# - $6  = feature_type
# - $8  = gene
# - $9  = locus_tag
# - $23 = attributes (fallback for locus_tag if needed)
tail -n +2 "$input" | awk -F'\t' '
BEGIN { OFS = FS }
{
  marker       = $1
  feature_type = $6

  # normalize gene
  gene = $8
  if (gene == "" || gene == "." || gene == "NA") gene = "unnamed"

  # prefer explicit locus_tag column ($9); fallback to attributes ($23)
  locus = $9
  if (locus == "" || locus == "." || locus == "NA") {
    n = split($23, A, ";")
    locus = "NA"
    for (i=1; i<=n; i++) {
      if (A[i] ~ /^locus_tag=/) {
        split(A[i], kv, "=")
        locus = kv[2]
        break
      }
    }
  }

  print marker, feature_type, gene, locus
}' >> "$out"

echo "✅ Gene summary written to: $out"

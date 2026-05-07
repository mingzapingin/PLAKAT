#!/usr/bin/env bash
set -euo pipefail
# 09_annotate_with_gff.sh
# Usage: 09_annotate_with_gff.sh <marker.bed> <gff_dir> <out_tsv> [MIN_FRAC_MARKER=0] [MIN_BP=0]
# Example: 09_annotate_with_gff.sh results/core/M.marinum_markers_on_ref.bed data/ref/gff results/core/M.marinum_marker_annotations.tsv 0.20 100

bed="$1"             # >=4-column BED (chrom start end name)
gff_dir="$2"         # directory containing reference_combined.gff
out="$3"
min_frac_marker="${4:-0}"
min_bp="${5:-0}"

command -v bedtools >/dev/null 2>&1 || { echo "❌ bedtools not found in PATH"; exit 127; }

gff="$gff_dir/reference_combined.gff"
[[ -s "$bed" ]] || { echo "❌ BED file missing/empty: $bed"; exit 1; }
[[ -s "$gff" ]] || { echo "❌ GFF file missing/empty: $gff"; exit 1; }

# Temps
gff_norm="$(mktemp)"
bed_norm="$(mktemp)"
tmp="$(mktemp)"
bed_chroms="$(mktemp)"
gff_chroms="$(mktemp)"
both="$(mktemp)"
cleanup() { rm -f "$gff_norm" "$bed_norm" "$tmp" "$bed_chroms" "$gff_chroms" "$both"; }
trap cleanup EXIT

# Normalize seqids (strip trailing ".<ver>") on BOTH sides
awk 'BEGIN{FS=OFS="\t"} /^#/{print; next} {sub(/\.[0-9]+$/,"",$1); print}' "$gff" > "$gff_norm"
awk 'BEGIN{FS=OFS="\t"} {sub(/\.[0-9]+$/,"",$1); print}' "$bed" > "$bed_norm"

# Preflight: ensure shared chromosome IDs exist
cut -f1 "$bed_norm" | sort -u > "$bed_chroms"
awk 'BEGIN{FS=OFS="\t"} !/^#/ {print $1}' "$gff_norm" | sort -u > "$gff_chroms"
comm -12 "$bed_chroms" "$gff_chroms" > "$both"
if [[ ! -s "$both" ]]; then
  echo "❌ No shared chrom IDs between BED and GFF after normalization."
  echo "•   Example BED chrom : $(head -1 "$bed_chroms")"
  echo "•   Example GFF chrom : $(head -1 "$gff_chroms")"
  echo "Hint: Make sure the BED was produced against the SAME FASTA used to create $gff"
  exit 2
fi

# Intersect (adds overlap length as last column; non-overlaps appear with 0)
bedtools intersect -a "$bed_norm" -b "$gff_norm" -wao > "$tmp"

# Header
echo -e "marker_id\tchrom\tstart\tend\tmarker_len\tfeature_type\tfeature_source\tgene\tlocus_tag\tname\tproduct\tfeat_id\tparent\tfeat_chrom\tfeat_start\tfeat_end\tfeat_len\tstrand\tphase\toverlap_bp\tfrac_marker\tfrac_feat\tattributes" > "$out"

# Parse & filter
awk -v MINF="$min_frac_marker" -v MINBP="$min_bp" 'BEGIN{OFS="\t"}
function attrget(s, key,   n,i,a){
  n = split(s, a, ";")
  for (i=1; i<=n; i++) {
    gsub(/^[ \t]+|[ \t]+$/, "", a[i])
    if (a[i] ~ "^" key "=") return substr(a[i], length(key)+2)
  }
  return ""
}
{
  # BED side (expects at least 4 cols)
  bed_chrom=$1; bed_start=$2+0; bed_end=$3+0; marker_id=$4
  marker_len = bed_end - bed_start

  # Intersect padding (from -wao)
  ovl_bp = $NF + 0

  # GFF side (index from tail)
  gff_attr  = $(NF-1)
  gff_phase = $(NF-2)
  gff_strand= $(NF-3)
  gff_end   = $(NF-5) + 0
  gff_start = $(NF-6) + 0
  gff_type  = $(NF-7)
  gff_src   = $(NF-8)
  gff_seq   = $(NF-9)

  # Skip no-overlap rows (from -wao) and malformed lengths
  if (ovl_bp < MINBP) next
  if (gff_type == "." || marker_len <= 0) next

  feat_len = gff_end - gff_start + 1
  if (feat_len <= 0) next

  frac_marker = ovl_bp / marker_len
  frac_feat   = ovl_bp / feat_len

  if (frac_marker < MINF) next

  gene      = attrget(gff_attr, "gene")
  locus_tag = attrget(gff_attr, "locus_tag")
  name      = attrget(gff_attr, "Name")
  product   = attrget(gff_attr, "product")
  feat_id   = attrget(gff_attr, "ID")
  parent    = attrget(gff_attr, "Parent")
  if (gene == "" && name != "") gene = name

  print marker_id, bed_chrom, bed_start, bed_end, marker_len,
        gff_type, gff_src,
        (gene==""?"NA":gene), (locus_tag==""?"NA":locus_tag),
        (name==""?"NA":name), (product==""?"NA":product),
        (feat_id==""?"NA":feat_id), (parent==""?"NA":parent),
        gff_seq, gff_start, gff_end, feat_len, gff_strand, gff_phase,
        ovl_bp, frac_marker, frac_feat, gff_attr
}' "$tmp" >> "$out"

# Summary
bed_rows=$(wc -l < "$bed" | awk '{print $1+0}')
annot_rows=$(( $(wc -l < "$out" | awk '{print $1+0}') - 1 ))

printf "▶ Annotated markers written: %s\n" "$out"
printf "•   BED rows (input) : %8d\n" "$bed_rows"
printf "•   annotated rows   : %8d\n" "$annot_rows"

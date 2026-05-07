#!/usr/bin/env bash
set -euo pipefail
# Usage: 07_prepare_reference.sh  <ref_list.txt>  <out_dir>

ref_list="$1"           # e.g.  M.marinum_ref.txt
ref_dir="$2"            # e.g.  data/ref
out_dir="$3"            # e.g.  data/genomes/

gbk_dir="$ref_dir/gbk"
gff_dir="$ref_dir/gff"
fa_dir="$ref_dir/fasta"
cache_dir="$out_dir/cache"   

mkdir -p "$gff_dir" "$fa_dir" "$cache_dir" "$gbk_dir"

while read -r line; do
  # skip blanks/comments
  [[ -z "$line" || "$line" =~ ^# ]] && continue

  # strip trailing " |REF" if present
  line="${line%%|REF*}"
  line="$(echo "$line" | sed 's/[[:space:]]\+$//')"

  # support both "ACC" and "chr:ACC@STRAIN" formats
  acc="$line"
  if [[ "$line" == *:*@* ]]; then
      acc="${line#*:}"      # drop "chr:" or "plasmid:"
      acc="${acc%@*}"       # drop "@STRAIN"
  fi

  gbk="$gbk_dir/${acc}.gbk"
  gff="$gff_dir/${acc}.gff"
  clean="$gff_dir/${acc}.clean.gff"
  fa="$fa_dir/${acc}.fasta"
  cached_fa="$cache_dir/${acc}.fasta"

  # ── download .gbk if missing ────────────────────────────────────────
  if [[ ! -s "$gbk" ]]; then
      echo "▶ Fetching $acc ..."
      efetch -db nuccore -id "$acc" -format gbwithparts -mode text > "$gbk"
  else
      echo "▶ Using cached $gbk"
  fi

  # ── convert to GFF3 (+FASTA) via EMBOSS seqret ─────────────────────
  if [[ ! -s "$gff" ]]; then
      seqret -sequence "$gbk" -feature -osformat gff3 -outseq "$gff"
  fi

  # strip FASTA section & comments (9-col only) → clean GFF
  if [[ ! -s "$clean" ]]; then
      awk '/^##FASTA/{exit} {print}' "$gff" | grep -v '^#' | awk 'NF==9' > "$clean"
  fi

  # ── Use cached FASTA if available, otherwise extract and cache ──────
  if [[ -s "$cached_fa" ]]; then
      echo "▶ Using cached FASTA for $acc"
      cp "$cached_fa" "$fa"
  elif [[ ! -s "$fa" ]]; then
      echo "▶ Extracting FASTA from $gbk"
      seqret -sequence "$gbk" -osformat fasta -auto -stdout > "$fa"
      cp "$fa" "$cached_fa"
  fi
done < "$ref_list"

rm -f "$fa_dir/reference_combined.fasta" "$gff_dir/reference_combined.gff"
cat "$fa_dir"/*.fasta > "$fa_dir/reference_combined.fasta"
cat "$gff_dir"/*.clean.gff > "$gff_dir/reference_combined.gff"

echo "▶ Reference prepared:"
echo "   - $fa_dir/reference_combined.fasta  (BLAST input)"
echo "   - $gff_dir/reference_combined.gff     (combined annotation)"

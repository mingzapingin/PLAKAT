#!/usr/bin/env bash
# 01_slice_genome.sh
# Usage: 01_slice_genome.sh  out_dir  input_fasta  window=500  step=250  [chrom_only=0]  [ids_file]
#   ids_file: file with one accession per line (exact match to the first token of FASTA header)

set -euo pipefail

outdir="$1"
input="$2"
size="${3:-500}"
step="${4:-250}"
chrom_only="${5:-0}"
ids_file="${6:-}"

tmp_ids=""
tmp_hdr=""
cleanup() {
  set +e
  if [ -n "${tmp_ids:-}" ] && [ -f "$tmp_ids" ]; then rm -f "$tmp_ids"; fi
  if [ -n "${tmp_hdr:-}" ] && [ -f "$tmp_hdr" ]; then rm -f "$tmp_hdr"; fi
  return 0
}
trap cleanup EXIT

command -v seqkit >/dev/null 2>&1 || { echo "❌ seqkit not found in PATH"; exit 127; }

mkdir -p "$outdir"
base=$(basename "$input" .fasta)
out_path="$outdir/${base}_frags_${size}bp.fasta"
rm -f "$out_path"

# Prefer explicit IDs list if provided (exact accession match by first header token)
if [[ -n "$ids_file" && -s "$ids_file" ]]; then
  tmp_ids=$(mktemp)
  awk -v idfile="$ids_file" '
    BEGIN{
      while ((getline < idfile) > 0) {
        gsub(/\r/,""); gsub(/^[ \t]+|[ \t]+$/,"");
        if (length($0) > 0) want[$0]=1;
      }
      close(idfile); keep=0;
    }
    /^>/{
      header=$0; sub(/^>/,"",header);
      split(header, a, /[[:space:]]+/);
      id=a[1]; keep = (id in want);
    }
    { if (keep) print }
  ' "$input" > "$tmp_ids"

  grep -q '^>' "$tmp_ids" || { echo "❌ No sequences matched IDs in: $ids_file" >&2; exit 2; }
  input="$tmp_ids"

elif [[ "$chrom_only" == "1" ]]; then
  tmp_hdr=$(mktemp)
  awk '/^>/{p = ($0 !~ /plasmid/i)} p' "$input" > "$tmp_hdr"
  grep -q '^>' "$tmp_hdr" || { echo "❌ chrom_only filter removed all sequences" >&2; exit 2; }
  input="$tmp_hdr"
fi

seqkit sliding -W "$size" -s "$step" "$input" -o "$out_path"
echo "▶ Wrote: $out_path"

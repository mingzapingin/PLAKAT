#!/usr/bin/env bash
set -euo pipefail
# 08_blast_vs_ref.sh
# Usage: 08_blast_vs_ref.sh <markers.fa> <ref_concat.fa> <out_bed> [threads=4] [task=megablast] [min_pident=95] [min_qcov=80]
#
# Output:
#   - <out_bed>                : 4-col BED (chrom, start, end, marker_id), 0-based, best hit per marker
#   - <out_bed>.unmapped.txt   : markers with no qualifying hit
#   - <out_bed>.multi.txt      : markers that had >1 qualifying HSP (flagged; best one used in BED)

markers="$1"
ref_fa="$2"
out_bed="$3"
threads="${4:-4}"
task="${5:-megablast}"
min_pid="${6:-95}"
min_qcov="${7:-80}"

# ---- sanity --------------------------------------------------------------
[[ -s "$markers" ]] || { echo "❌ Missing/empty markers FASTA: $markers"; exit 2; }
[[ -s "$ref_fa"  ]] || { echo "❌ Missing/empty reference FASTA: $ref_fa"; exit 2; }
case "$out_bed" in *.bed) : ;; *) echo "❌ Output must be .bed (got: $out_bed)"; exit 2;; esac
mkdir -p "$(dirname "$out_bed")"

# ---- lightweight stats ---------------------------------------------------
mcount=$(grep -c '^>' "$markers" || echo 0)
rseqs=$(grep -c '^>' "$ref_fa" || echo 0)
rbp="NA"
if command -v seqkit >/dev/null 2>&1; then
  rbp=$(seqkit stats -T "$ref_fa" 2>/dev/null | awk 'NR==2{print $5}')
fi

db_dir="data/db"; mkdir -p "$db_dir"
ref_stem=$(basename "$ref_fa" .fasta)
db_path="$db_dir/${ref_stem}"

echo "▶ Reference DB prep"
echo "•   DB path         : $db_path"
echo "•   DB input FASTA  : $ref_fa"
echo "•   DB sequences    : $rseqs"
echo "•   DB total bp     : $rbp"
echo "•   markers (queries): $mcount"

# ---- (re)build database if missing --------------------------------------
if [[ ! -f ${db_path}.nin && ! -f ${db_path}.00.nin ]]; then
  makeblastdb -in "$ref_fa" -dbtype nucl -out "$db_path"
fi

# ---- BLAST mapping -------------------------------------------------------
echo "▶ BLAST mapping markers → reference"
echo "•   thresholds      : pident≥${min_pid}%, qcov≥${min_qcov}%, task=${task}"
echo "•   masking         : dust=no, soft_masking=false"
echo "•   threads         : ${threads}"

start=$(date +%s)
tmp_tsv="$(mktemp)"
trap 'rm -f "$tmp_tsv" "$out_bed.ids" "$out_bed.unmapped.txt" "$out_bed.multi.txt.tmp"' EXIT

blastn \
  -query "$markers" \
  -db "$db_path" \
  -task "$task" \
  -dust no -soft_masking false \
  -perc_identity "$min_pid" \
  -qcov_hsp_perc "$min_qcov" \
  -num_threads "$threads" \
  -outfmt '6 qseqid sseqid pident length qlen sstart send evalue bitscore qcovhsp' \
  > "$tmp_tsv"

# ---- choose best HSP per marker; flag multi-hit markers ------------------
# Keep the record with the highest bitscore; tie-break by pident then length.
# Build a set of markers that had >1 qualifying HSP overall.
awk -v OFS="\t" '
{
  q=$1; s=$2; pid=$3+0; len=$4+0; qlen=$5+0; ss=$6+0; se=$7+0; e=$8; bit=$9+0; qc=$10+0;

  nhsp[q]++
  if (!(q in best) || bit > b_bit[q] || (bit==b_bit[q] && (pid > b_pid[q] || (pid==b_pid[q] && len > b_len[q])))) {
    best[q]=$0; b_bit[q]=bit; b_pid[q]=pid; b_len[q]=len;
  }
}
END{
  for (q in nhsp) if (nhsp[q]>1) print q > "'"$out_bed"'.multi.txt.tmp"
  for (q in best) print best[q]
}
' "$tmp_tsv" > "$out_bed.ids"

# ---- convert best hits to BED (0-based), strip version suffix ------------

awk -v OFS="\t" '
function stripver(x){ sub(/\.[0-9]+$/,"",x); return x }
{
  q=$1; s=$2; ss=$6+0; se=$7+0;
  chr=stripver(s);
  start = (ss<se ? ss : se) - 1;
  end   = (ss<se ? se : ss);
  if (start<0) start=0;
  print chr, start, end, q
}
' "$out_bed.ids" > "$out_bed"

# ---- unmapped + multi counts --------------------------------------------
# list unmapped queries (present in markers but absent among best-hit ids)
awk '/^>/{sub(/^>/,""); print $1}' "$markers" | sort -u > "$out_bed".allq.tmp
cut -f1 "$out_bed.ids" | sort -u > "$out_bed".mappedq.tmp
comm -23 "$out_bed".allq.tmp "$out_bed".mappedq.tmp > "$out_bed".unmapped.txt || true

multi_n=0
if [[ -f "$out_bed.multi.txt.tmp" ]]; then
  sort -u "$out_bed.multi.txt.tmp" > "$out_bed".multi.txt
  multi_n=$(wc -l < "$out_bed".multi.txt | awk "{print \$1}")
  rm -f "$out_bed.multi.txt.tmp"
fi

mapped_n=$(wc -l < "$out_bed" 2>/dev/null | awk '{print $1+0}')
unmap_n=$(wc -l < "$out_bed".unmapped.txt 2>/dev/null | awk '{print $1+0}')

dur=$(( $(date +%s) - start ))

printf "▶ Wrote BED: %s\n" "$out_bed"
printf "•   markers total  : %8d\n" "$mcount"
printf "•   mapped (best)  : %8d\n" "$mapped_n"
printf "•   unmapped       : %8d\n" "$unmap_n"
printf "•   multi-hit flag : %8d\n" "$multi_n"
printf "•   runtime        : %dm%02ds\n" "$((dur/60))" "$((dur%60))"

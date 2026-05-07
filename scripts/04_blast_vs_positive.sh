#!/usr/bin/env bash
set -euo pipefail
# Usage: 04_blast_vs_positive.sh <force_db 0|1> <centroids.fa> <pos_fasta> <out_tsv> [threads=4] [task=megablast]

force_db="$1"; shift
centroids="$1"; shift
pos_fasta="$1"; shift
out_tsv="$1"; shift

# optional args with guarding shifts
if [[ $# -ge 1 ]]; then
  threads="$1"; shift
else
  threads="4"
fi
if [[ $# -ge 1 ]]; then
  task="$1"; shift
else
  task="megablast"
fi

# sanity
[[ -s "$centroids" ]] || { echo "❌ Missing/empty centroids FASTA: $centroids"; exit 2; }
[[ -s "$pos_fasta"  ]] || { echo "❌ Missing/empty positives FASTA: $pos_fasta"; exit 2; }
case "$out_tsv" in *.tsv) ;; *) echo "❌ Output must be a .tsv (got: $out_tsv)"; exit 2 ;; esac
if [[ "$(basename "$pos_fasta")" == "NEGATIVE_TEMP.fasta" ]]; then
  echo "❌ You passed the NEGATIVE DB as positives DB: $pos_fasta"; exit 2
fi

# validate task for blastn
case "$task" in
  megablast|dc-megablast|blastn|blastn-short|rmblastn) ;;
  *) echo "❌ Invalid task: '$task' (allowed: megablast, dc-megablast, blastn, blastn-short, rmblastn)"; exit 2 ;;
esac

mkdir -p "$(dirname "$out_tsv")"

# lightweight stats
qcount=$(grep -c '^>' "$centroids" || echo 0)
db_seqs=$(grep -c '^>' "$pos_fasta" || echo 0)
db_bp="NA"
if command -v seqkit >/dev/null 2>&1; then
  db_bp=$(seqkit stats -T "$pos_fasta" 2>/dev/null | awk 'NR==2{print $5}')
fi

db_dir="data/db"; mkdir -p "$db_dir"
db_stem=$(basename "$pos_fasta" .fasta)
db_path="$db_dir/${db_stem}"

echo "▶ Positives DB prep"
echo "•   DB path          : $db_path"
echo "•   DB input FASTA   : $pos_fasta"
echo "•   DB sequences     : $db_seqs"
echo "•   DB total bp      : $db_bp"
echo "•   Query centroids  : $qcount"

# (re)build DB
if [[ "$force_db" == "1" ]]; then rm -f "${db_path}".{nin,nhr,nsq}; fi
if [[ ! -f ${db_path}.nin ]]; then makeblastdb -in "$pos_fasta" -dbtype nucl -out "$db_path"; fi

echo "▶ BLAST centroids → positives"
echo "•   thresholds      : pident≥98%, qcov≥80%, task=${task}"
echo "•   masking         : dust=no, soft_masking=false"
echo "•   threads         : ${threads}"

start=$(date +%s)
blastn \
  -query "$centroids" \
  -db "$db_path" \
  -out "$out_tsv" \
  -outfmt '6 qseqid sseqid pident length qlen qstart qend slen sstart send evalue bitscore qcovhsp stitle' \
  -task "$task" \
  -dust no -soft_masking false \
  -perc_identity 98 \
  -qcov_hsp_perc 80 \
  -num_threads "$threads"
dur=$(( $(date +%s) - start ))
hcount=$(wc -l < "$out_tsv" 2>/dev/null | awk '{print $1+0}')

printf "▶ Centroids→positives BLAST written: %s\n" "$out_tsv"
printf "•   queries         : %d\n" "$qcount"
printf "•   hits (rows)     : %d\n" "$hcount"
printf "•   runtime         : %dm%02ds\n" "$((dur/60))" "$((dur%60))"

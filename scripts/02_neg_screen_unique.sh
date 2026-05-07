#!/usr/bin/env bash
# 02_neg_screen_unique.sh
# Usage: 02_neg_screen_unique.sh <force_db 0|1> <windows_fasta> <negative_fasta> <unique_out_fasta> <out_tsv> [threads=4]
# Does:
#   1) Build (or reuse) a BLAST nucleotide DB for the negatives
#   2) BLAST windows vs negatives, output only qseqid (one per hit)
#   3) Compute set-difference (all windows - hit windows)
#   4) Extract no-hit (unique) windows to FASTA

set -euo pipefail

force_db="$1"; shift
windows="$1"; shift
neg_fasta="$1"; shift
unique_out="$1"; shift
out_tsv="$1"; shift
threads="${1:-4}"

command -v makeblastdb >/dev/null 2>&1 || { echo "❌ makeblastdb not found"; exit 127; }
command -v blastn      >/dev/null 2>&1 || { echo "❌ blastn not found"; exit 127; }
command -v seqkit      >/dev/null 2>&1 || { echo "❌ seqkit not found"; exit 127; }

# lightweight stats
qcount=$(grep -c '^>' "$windows" || echo 0)
db_seqs=$(grep -c '^>' "$neg_fasta" || echo 0)
db_bp="NA"
if command -v seqkit >/dev/null 2>&1; then
  db_bp=$(seqkit stats -T "$neg_fasta" 2>/dev/null | awk 'NR==2{print $5}')
fi

db_dir="data/db"
mkdir -p "$db_dir"
stem=$(basename "$neg_fasta" .fasta)
db_path="$db_dir/${stem}"

echo "▶ Negative DB prep"
echo "•   DB path          : $db_path"
echo "•   DB input FASTA   : $neg_fasta"
echo "•   DB sequences     : $db_seqs"
echo "•   DB total bp      : $db_bp"
echo "•   Query windows    : $qcount"

# (Re)build DB
if [[ "$force_db" == "1" ]]; then
  echo "▶ Rebuilding DB for $stem"
  rm -f "${db_path}".*
fi
if [[ ! -f ${db_path}.nhr && ! -f ${db_path}.00.nhr ]]; then
  makeblastdb -in "$neg_fasta" -dbtype nucl -out "$db_path"
else
  echo "▶ Using cached DB $db_path"
fi

mkdir -p "$(dirname "$unique_out")"
tmp_all=$(mktemp)
tmp_hit=$(mktemp)
tmp_ids=$(mktemp)
tmp_pat=$(mktemp)

# all window IDs (names)
seqkit fx2tab -n "$windows" | cut -f1 | sort -u > "$tmp_all"

# BLAST all windows vs negatives; capture a RICH TSV for ALL_HITS
OUTFMT='6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp slen stitle'
blastn -task megablast \
  -query "$windows" \
  -db "$db_path" \
  -evalue 1e-10 \
  -perc_identity 85 \
  -num_threads "$threads" \
  -dust yes \
  -outfmt "$OUTFMT" \
  2> >(grep -v 'Examining 5 or more matches is recommended' >&2) \
| tee "$out_tsv" \
| cut -f1 \
| sort -u > "$tmp_hit"

# IDs with no hits in negatives
comm -23 "$tmp_all" "$tmp_hit" > "$tmp_ids"

# Extract exact-ID matches to FASTA (anchor the patterns)
awk '{print "^"$0"$"}' "$tmp_ids" > "$tmp_pat"
seqkit grep -r -n -f "$tmp_pat" "$windows" > "$unique_out"

# --- pretty summary ---
total=$(wc -l < "$tmp_all" | awk '{print $1}')
hits=$(wc -l < "$tmp_hit" | awk '{print $1}')
unique=$(( total - hits ))

pct_hit=$(awk -v h="$hits" -v t="$total" 'BEGIN{printf("%.2f%%", (t?100*h/t:0))}')
pct_unique=$(awk -v u="$unique" -v t="$total" 'BEGIN{printf("%.2f%%", (t?100*u/t:0))}')

printf "▶ Unique windows written: %s\n" "$unique_out"
printf "•   total  : %8d (100%%)\n" "$total"
printf "•   hits   : %8d (%s)\n"   "$hits"   "$pct_hit"
printf "•   unique : %8d (%s)\n"   "$unique" "$pct_unique"

rm -f "$tmp_all" "$tmp_hit" "$tmp_ids" "$tmp_pat"

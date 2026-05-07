#!/usr/bin/env bash
set -euo pipefail
# 03_dedup_cluster.sh  <unique_fasta>  <out_dir>  [threads=4]  [id=0.95]
# Clusters 500bp windows: seqkit rmdup → cd-hit-est

in_fa="$1"
outdir="$2"
threads="${3:-4}"
id="${4:-0.95}"          # allow override (e.g. 0.9)

mkdir -p "$outdir"
command -v seqkit     >/dev/null 2>&1 || { echo "❌ seqkit not found"; exit 127; }
command -v cd-hit-est >/dev/null 2>&1 || { echo "❌ cd-hit-est not found"; exit 127; }

# File stems/paths
stem="$(basename "$in_fa" .fasta)"
dedup_fa="$outdir/${stem}_dedup.fasta"
clust_fa="$outdir/${stem}_95.fasta"     # keep old name for compatibility
log_fp="$outdir/${stem}_cdhit.log"

# Helper to count FASTA records safely
count_fa() { grep -c '^>' "$1" 2>/dev/null || echo 0; }

# Choose a valid word length for cd-hit-est based on identity
# cd-hit-est rules of thumb
# >=0.95 → n=10; 0.90–<0.95 → n=10; 0.88–<0.90 → n=9; 0.85–<0.88 → n=8; 0.80–<0.85 → n=7; 0.75–<0.80 → n=6; else n=5
pick_n() {
  awk -v c="$1" 'BEGIN{
    n = (c>=0.95)?10 : (c>=0.90)?10 : (c>=0.88)?9 : (c>=0.85)?8 : (c>=0.80)?7 : (c>=0.75)?6 : 5;
    print n;
  }'
}
nword="$(pick_n "$id")"

echo "▶ seqkit rmdup  →  $dedup_fa"
seqkit rmdup -s "$in_fa" -o "$dedup_fa"

in_n="$(count_fa "$in_fa")"
dedup_n="$(count_fa "$dedup_fa")"
if [[ "$dedup_n" -eq 0 ]]; then
  echo "❌ No sequences left after rmdup. Input had $in_n records. Aborting."
  exit 2
fi

# Clean stale outputs so cd-hit-est can overwrite freely
rm -f "$clust_fa" "${clust_fa}.clstr" "$log_fp"

echo "▶ CD-HIT-EST $(awk -v c="$id" 'BEGIN{printf("%.2f",c)*100}') %  →  $clust_fa"
set +e
cd-hit-est -i "$dedup_fa" -o "$clust_fa" \
           -c "$id" -n "$nword" \
           -T "$threads" -M 0 \
           -d 0 2>"$log_fp"
status=$?
set -e

if [[ $status -ne 0 ]]; then
  echo "❌ cd-hit-est failed (exit $status). First lines of log:"
  head -n 30 "$log_fp" || true
  # A couple of common auto-fixes:
  if grep -qi "word length" "$log_fp"; then
    echo "↻ Retrying with relaxed word length (-n 9) keeping -c=$id"
    set +e
    cd-hit-est -i "$dedup_fa" -o "$clust_fa" -c "$id" -n 9 -T "$threads" -M 0 -d 0 2>>"$log_fp"
    status=$?
    set -e
  fi
  if [[ $status -ne 0 ]]; then
    echo "Log saved to: $log_fp"
    exit 2
  fi
fi

final_n="$(count_fa "$clust_fa")"
dup_removed=$(( in_n - dedup_n ))
merged_by_cdhit=$(( dedup_n - final_n ))

pct() { awk -v a="$1" -v b="$2" 'BEGIN{printf("%.2f%%", (b?100*a/b:0))}'; }
pct_dup_removed="$(pct "$dup_removed" "$in_n")"
pct_after_dedup="$(pct "$dedup_n" "$in_n")"
pct_final="$(pct "$final_n" "$in_n")"
pct_merged_cdhit="$(pct "$merged_by_cdhit" "$dedup_n")"

printf "▶ Dedup & cluster summary for %s\n" "$stem"
printf "•   input windows       : %8d (100%%)\n"     "$in_n"
printf "•   after rmdup         : %8d (%s)\n"        "$dedup_n" "$pct_after_dedup"
printf "•     └─ duplicates rm  : %8d (%s of input)\n" "$dup_removed" "$pct_dup_removed"
printf "•   after CD-HIT (%.0f%%) : %8d (%s of input)\n" "$(awk -v c="$id" 'BEGIN{print c*100}')" "$final_n" "$pct_final"
printf "•     └─ merged by cd-hit : %8d (%s of rmdup)\n"  "$merged_by_cdhit" "$pct_merged_cdhit"

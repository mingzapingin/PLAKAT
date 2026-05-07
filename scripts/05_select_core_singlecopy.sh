#!/usr/bin/env bash
# 05_select_core_singlecopy.sh  (robust)
# Usage:
#   05_select_core_singlecopy.sh <blast.tsv> <centroids.fa> <core_out.fa> <chr_ids.txt> [min_pident=98] [min_qcov=80] [min_chr(K)=ALL] [min_frac=]
# Notes:
#   - blast.tsv outfmt 6 must include: qseqid sseqid pident length qlen (in that order)
#   - chr_ids.txt: one allowed sseqid per line (exact string match)
set -euo pipefail

blast_tsv="$1"
centroids_fa="$2"
core_out="$3"
chr_ids="$4"
min_pid="${5:-98}"
min_qcov="${6:-80}"
want_k="${7:-}"          # integer K-of-N (wins if given)
want_frac="${8:-}"       # fraction of N (used if K not given)

[[ -s "$blast_tsv"    ]] || { echo "❌ Missing BLAST table: $blast_tsv"; exit 2; }
[[ -s "$centroids_fa" ]] || { echo "❌ Missing centroids FASTA: $centroids_fa"; exit 2; }
[[ -s "$chr_ids"      ]] || { echo "❌ Missing chr IDs list: $chr_ids"; exit 2; }


# Count required chromosomes (N)
N=$(awk 'NF>0{n++} END{print (n+0)}' "$chr_ids")

# Derive K (required hits)
if [[ -n "${want_k:-}" ]]; then
  K="$want_k"
elif [[ -n "${want_frac:-}" ]]; then
  K=$(awk -v n="$N" -v f="$want_frac" 'BEGIN{
         if (f<0) f=0; if (f>1) f=1;
         k = int(f*n + 0.999999); if (k<1) k=1; if (k>n) k=n; print k
       }')
else
  K="$N"
fi
# Clamp
if [[ "$K" -lt 1 ]]; then K=1; fi
if [[ "$K" -gt "$N" ]]; then K="$N"; fi

ids_tmp="$(mktemp)"; trap 'rm -f "$ids_tmp"' EXIT

# Quick sanity: do sseqids in BLAST overlap chr_ids?
have_overlap=$(
  awk 'NR==FNR{a[$0]=1; next} ($2 in a){ok=1} END{print (ok?1:0)}' "$chr_ids" "$blast_tsv"
)
if [[ "$have_overlap" -eq 0 ]]; then
  ex_bed=$(awk 'NR==1{print $2; exit}' "$blast_tsv")
  ex_chr=$(awk 'NF>0{print; exit}' "$chr_ids")
  echo "❌ No overlap between BLAST sseqid and chr_ids." >&2
  echo "• example BLAST sseqid: $ex_bed" >&2
  echo "• example chr_ids item : $ex_chr" >&2
  exit 2
fi

# Pick queries that:
#   • hit allowed sseqids only
#   • meet thresholds (pident ≥ min_pid AND qcov ≥ min_qcov)
#   • are single-copy per subject (≤1 HSP per q,s)
#   • are present on ≥K distinct required chromosomes
awk -v MINPID="$min_pid" -v MINQC="$min_qcov" -v CHR="$chr_ids" -v K="$K" '
BEGIN{
  while ((getline < CHR) > 0) { gsub(/\r/,""); gsub(/^[ \t]+|[ \t]+$/,""); if (length($0)) allow[$0]=1; nchr++ }
  close(CHR)
}
{
  q=$1; s=$2; pid=$3+0; hlen=$4+0; qlen=$5+0
  if (!(s in allow)) { miss_s++ ; next }
  qcov = (qlen>0 ? 100.0*hlen/qlen : 0)
  if (pid < MINPID || qcov < MINQC) { fail_thr++ ; next }
  nhsp[q SUBSEP s]++        # HSPs per (q,s) for single-copy check
  seen[q,s]=1               # distinct subject coverage per q
  pass_thr=1
}
END{
  # single-copy filter
  for (ks in nhsp) { split(ks,t,SUBSEP); q=t[1]; if (nhsp[ks]>1) bad[q]=1 }
  # presence count
  for (qs in seen) { split(qs,t,SUBSEP); q=t[1]; hitcnt[q]++ }
  # output keepers
  kept=0
  for (q in hitcnt) if (hitcnt[q]>=K && !bad[q]) { print q; kept++ }
  # diagnostics to stderr
  printf("▶ Core selection diagnostics:\n") > "/dev/stderr"
  printf("•   allowed chr (N)   : %d\n", nchr) > "/dev/stderr"
  printf("•   threshold fails   : %d\n", (fail_thr+0)) > "/dev/stderr"
  printf("•   wrong chr sseqid  : %d\n", (miss_s+0)) > "/dev/stderr"
  sc_multi=0; for (ks in nhsp) if (nhsp[ks]>1) sc_multi++
  printf("•   multi-HSP (q,s)   : %d\n", sc_multi) > "/dev/stderr"
  printf("•   kept IDs          : %d\n", kept) > "/dev/stderr"
}
' "$blast_tsv" > "$core_out.ids"

# If no IDs, create empty FASTA and still summarize nicely
if [[ ! -s "$core_out.ids" ]]; then
  : > "$core_out"
else
  # Extract the kept centroids
  seqkit grep -f "$core_out.ids" "$centroids_fa" > "$core_out"
fi

# Counts (sanitize to integers)
kept=$(awk 'BEGIN{c=0} /^>/{c++} END{print c+0}' "$core_out")
tot=$(awk 'BEGIN{c=0} /^>/{c++} END{print c+0}' "$centroids_fa")

printf "▶ Selected core single-copy markers → %s\n" "$core_out"
printf "•   required chr   : %s of %s\n" "$K" "$N"
printf "•   thresholds     : pident≥%s, qcov≥%s\n" "$min_pid" "$min_qcov"
printf "•   kept markers   : %s (of %s)\n" "$kept" "$tot"

rm -f "$core_out.ids"

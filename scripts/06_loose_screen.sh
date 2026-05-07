#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   06_loose_screen.sh CORE_FASTA NEG_FASTA_OR_DB FINAL_FASTA_OUT NEG_TSV_OUT THREADS [PID=80] [QCOV=50]
#
# Example:
#   06_loose_screen.sh \
#     results/core/M.marinum_core_singlecopy.fasta \
#     data/NEGATIVE_TEMP.fasta \
#     results/core/M.marinum_final_markers.fasta \
#     results/blast/M.marinum_centroids_vs_neg.tsv \
#     4 80 50

CORE_FASTA="$1"
NEG_SRC="$2"          # FASTA or BLAST DB prefix
FINAL_FASTA_OUT="$3"
NEG_TSV_OUT="$4"
THREADS="$5"
PID="${6:-80}"
QCOV="${7:-50}"

# Tools we rely on:
for exe in blastn makeblastdb seqkit awk sort uniq tee; do
  command -v "$exe" >/dev/null 2>&1 || { echo "ERROR: Missing required tool: $exe" >&2; exit 3; }
done

[[ -s "$CORE_FASTA" ]] || { echo "ERROR: CORE_FASTA not found or empty: $CORE_FASTA" >&2; exit 4; }

# Ensure output dirs exist
mkdir -p "$(dirname "$FINAL_FASTA_OUT")" "$(dirname "$NEG_TSV_OUT")"

# Decide if NEG_SRC is a DB prefix or a FASTA we need to index
NEG_DB_PREFIX="$NEG_SRC"
if [[ -f "$NEG_SRC" ]]; then
  # NEG_SRC is a FASTA; make a DB alongside it (idempotent)
  NEG_DB_PREFIX="${NEG_SRC}.blastdb"
  if [[ ! -f "${NEG_DB_PREFIX}.nsq" ]]; then
    echo "▶ Building BLAST DB for negatives at: ${NEG_DB_PREFIX} (from $NEG_SRC)"
    makeblastdb -in "$NEG_SRC" -dbtype nucl -parse_seqids -out "$NEG_DB_PREFIX" 1>&2
  fi
else
  # Assume it's a DB prefix; sanity check
  if [[ ! -f "${NEG_DB_PREFIX}.nsq" && ! -f "${NEG_DB_PREFIX}.00.nsq" ]]; then
    echo "ERROR: NEG_SRC does not look like a BLAST DB prefix or a FASTA: $NEG_SRC" >&2
    exit 5
  fi
fi

# Temp files
TMP_DIR="$(mktemp -d -t step6_XXXXXX)"
trap 'rm -rf "$TMP_DIR"' EXIT
BAD_IDS="$TMP_DIR/drop_ids.txt"

echo "▶ Running loose NEG screen with thresholds: pident>=${PID}, qcovhsp>=${QCOV}"
echo "   core: $CORE_FASTA"
echo "   db  : $NEG_DB_PREFIX"
echo "   tsv : $NEG_TSV_OUT"
echo "   out : $FINAL_FASTA_OUT"
echo "   thr : $THREADS threads"

# Rich outfmt to power ALL_HITS assembly later; qseqid first so cut -f1 works for drop list
OUTFMT='6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp slen stitle'

# BLAST with thresholds applied so TSV == exact evidence used for dropping
# (No max_target_seqs cap; we want all qualifying hits.)
blastn \
  -db "$NEG_DB_PREFIX" \
  -query "$CORE_FASTA" \
  -task blastn \
  -evalue 1e-6 \
  -perc_identity "$PID" \
  -qcov_hsp_perc "$QCOV" \
  -num_threads "$THREADS" \
  -dust yes \
  -outfmt "$OUTFMT" \
| tee "$NEG_TSV_OUT" \
| cut -f1 \
| sort -u > "$BAD_IDS"

# Remove bad IDs from CORE_FASTA → FINAL_FASTA_OUT
if [[ -s "$BAD_IDS" ]]; then
  seqkit grep -v -f "$BAD_IDS" "$CORE_FASTA" > "$FINAL_FASTA_OUT"
else
  # No drops; copy through
  cp "$CORE_FASTA" "$FINAL_FASTA_OUT"
fi

# Quick summary
TOTAL_Q=$(grep -c '^>' "$CORE_FASTA" || true)
DROPPED=$(wc -l < "$BAD_IDS" 2>/dev/null || echo 0)
KEPT=$(grep -c '^>' "$FINAL_FASTA_OUT" || true)

echo "▶ Loose NEG screen summary:"
echo "   total_markers : $TOTAL_Q"
echo "   dropped       : $DROPPED"
echo "   kept          : $KEPT"
echo "   evidence_tsv  : $NEG_TSV_OUT"

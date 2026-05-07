#!/usr/bin/env bash
# 11_design_primers_from_markers.sh
# Verbose wrapper around primer3_core to design N primer pairs per marker
# and aggregate results into a single TSV.

set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  design_primers_from_markers.sh <final_markers.fasta> <out.tsv> [options]

Required:
  final_markers.fasta   Combined FASTA of all marker sequences (e.g., 500 bp each)
  out.tsv               Destination TSV file

Options (all optional; defaults shown):
  --num-return N              [default: 5]
  --product-size MIN-MAX      [default: 200-450]
  --tm MIN-MAX                [default: 65-75]        # °C, nearest-neighbor model
  --gc MIN-MAX                [default: 40-60]        # percent
  --len MIN,OPT,MAX           [default: 16,18,24]     # nt
  --gc-clamp N                [default: 1]
  --max-poly-x N              [default: 4]
  --dna-conc N                [default: 250]          # nM
  --salt-mono N               [default: 50]           # mM
  --salt-di N                 [default: 1.5]          # mM
  --dntp N                    [default: 0.6]          # mM
  --thermo-path DIR           # primer3_config directory (auto-detected if possible)
  --quiet                     # reduce verbosity
  -h | --help

Notes:
- Script prints progress and chosen parameters for each marker unless --quiet.
- Output TSV columns:
    marker_id, marker_seq, pair_index, left_seq, right_seq, product_size,
    left_len, right_len, left_tm, right_tm, left_gc, right_gc,
    left_self_any_th, left_self_end_th, right_self_any_th, right_self_end_th,
    pair_compl_any_th, pair_compl_end_th, pair_penalty
USAGE
}

# ---- Defaults ----
NUM_RETURN=5
PRODUCT_RANGE="200-450"
TM_MIN=65
TM_MAX=75
GC_MIN=40
GC_MAX=60
LEN_MIN=16
LEN_OPT=18
LEN_MAX=24
GC_CLAMP=1
MAX_POLY_X=4
DNA_CONC=250
SALT_MONO=50
SALT_DI=1.5
DNTP=0.6
THERMO_PATH=""
QUIET=0

log() { [ "$QUIET" -eq 0 ] && echo "[INFO] $*" >&2 || true; }

# ---- Parse args ----
if [ $# -lt 2 ]; then usage; exit 1; fi
FASTA=$1; shift
OUT_TSV=$1; shift

while [ $# -gt 0 ]; do
  case "$1" in
    --num-return) NUM_RETURN=$2; shift 2;;
    --product-size) PRODUCT_RANGE=$2; shift 2;;
    --tm) TM_MIN=$(echo "$2" | awk -F- '{print $1}'); TM_MAX=$(echo "$2" | awk -F- '{print $2}'); shift 2;;
    --gc) GC_MIN=$(echo "$2" | awk -F- '{print $1}'); GC_MAX=$(echo "$2" | awk -F- '{print $2}'); shift 2;;
    --len) LEN_MIN=$(echo "$2" | awk -F, '{print $1}');
           LEN_OPT=$(echo "$2" | awk -F, '{print $2}');
           LEN_MAX=$(echo "$2" | awk -F, '{print $3}'); shift 2;;
    --gc-clamp) GC_CLAMP=$2; shift 2;;
    --max-poly-x) MAX_POLY_X=$2; shift 2;;
    --dna-conc) DNA_CONC=$2; shift 2;;
    --salt-mono) SALT_MONO=$2; shift 2;;
    --salt-di) SALT_DI=$2; shift 2;;
    --dntp) DNTP=$2; shift 2;;
    --thermo-path) THERMO_PATH=$2; shift 2;;
    --quiet) QUIET=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage; exit 1;;
  esac
done

# ---- Checks ----
if [ ! -s "$FASTA" ]; then echo "[ERROR] FASTA not found: $FASTA" >&2; exit 1; fi
if ! command -v primer3_core >/dev/null 2>&1; then
  echo "[ERROR] primer3_core not found in PATH." >&2
  echo "        Install Primer3 and ensure 'primer3_core' is available." >&2
  exit 1
fi

# Try to auto-detect primer3 thermodynamic parameters path if not given
if [ -z "$THERMO_PATH" ]; then
  P3C="${PRIMER3_CONFIG:-}"  # safe even when PRIMER3_CONFIG is unset
  for d in \
    "$P3C" \
    "/opt/homebrew/share/primer3/primer3_config" \
    "/usr/local/share/primer3/primer3_config" \
    "/usr/share/primer3/primer3_config"
  do
    if [ -n "$d" ] && [ -d "$d" ]; then
      THERMO_PATH="$d"
      break
    fi
  done
fi
if [ -n "$THERMO_PATH" ] && [ -d "$THERMO_PATH" ]; then
  log "Using thermodynamic parameters at: $THERMO_PATH"
else
  log "Thermo parameter directory not found; relying on primer3 defaults (if compiled in)."
  THERMO_PATH=""
fi

# Compute OPT values from min/max if needed
TM_OPT=$(awk -v a="$TM_MIN" -v b="$TM_MAX" 'BEGIN{printf("%.1f",(a+b)/2.0)}')
GC_OPT=$(awk -v a="$GC_MIN" -v b="$GC_MAX" 'BEGIN{printf("%.1f",(a+b)/2.0)}')

# Prepare output
TMPDIR=$(mktemp -d -t p3wrap.XXXXXX)
trap 'rm -rf "$TMPDIR"' EXIT

echo -e "marker_id\tmarker_seq\tproduct_seq\tpair_index\tleft_seq\tright_seq\tproduct_size\tleft_len\tright_len\tleft_tm\tright_tm\tleft_gc\tright_gc\tleft_self_any_th\tleft_self_end_th\tright_self_any_th\tright_self_end_th\tpair_compl_any_th\tpair_compl_end_th\tpair_penalty" > "$OUT_TSV"

# Count markers
N_MARKERS=$(grep -c "^>" "$FASTA" || true)
echo "▶Found $N_MARKERS markers in $FASTA"
echo "▶Parameters: num_return=$NUM_RETURN, product_size=$PRODUCT_RANGE, Tm=$TM_MIN-$TM_MAX (opt $TM_OPT), GC=$GC_MIN-$GC_MAX (opt $GC_OPT), len=$LEN_MIN,$LEN_OPT,$LEN_MAX"

# Convert FASTA to (id \t seq) lines and iterate
awk '
  BEGIN{ORS=""; id=""; seq=""}
  /^>/{
    if(length(seq)>0){print id "\t" seq "\n"}
    id=substr($0,2); seq="";
    next
  }
  {
    gsub(/[ \t\r]/,"");
    if(length($0)>0){seq=seq toupper($0)}
  }
  END{
    if(length(seq)>0){print id "\t" seq "\n"}
  }
' "$FASTA" | while IFS=$'\t' read -r MARKER_ID MARKER_SEQ; do
  [ -z "$MARKER_ID" ] && continue

  # sanitize marker id for filenames (avoid ':' on macOS)
  SAFE_ID=$(printf "%s" "$MARKER_ID" | tr '/: ' '___')

  P3IN="$TMPDIR/$SAFE_ID.p3in"
  P3OUT="$TMPDIR/$SAFE_ID.p3out"
  echo "▶ Designing primers for: $MARKER_ID (len=${#MARKER_SEQ})"


  {
    echo "SEQUENCE_ID=$MARKER_ID"
    echo "SEQUENCE_TEMPLATE=$MARKER_SEQ"
    echo "PRIMER_TASK=pick_pcr_primers"
    echo "PRIMER_NUM_RETURN=$NUM_RETURN"
    echo "PRIMER_PRODUCT_SIZE_RANGE=$PRODUCT_RANGE"
    echo "PRIMER_MIN_SIZE=$LEN_MIN"
    echo "PRIMER_OPT_SIZE=$LEN_OPT"
    echo "PRIMER_MAX_SIZE=$LEN_MAX"
    echo "PRIMER_MIN_TM=$TM_MIN"
    echo "PRIMER_OPT_TM=$TM_OPT"
    echo "PRIMER_MAX_TM=$TM_MAX"
    echo "PRIMER_MIN_GC=$GC_MIN"
    echo "PRIMER_OPT_GC_PERCENT=$GC_OPT"
    echo "PRIMER_MAX_GC=$GC_MAX"
    echo "PRIMER_GC_CLAMP=$GC_CLAMP"
    echo "PRIMER_MAX_POLY_X=$MAX_POLY_X"
    echo "PRIMER_DNA_CONC=$DNA_CONC"
    echo "PRIMER_SALT_MONOVALENT=$SALT_MONO"
    echo "PRIMER_SALT_DIVALENT=$SALT_DI"
    echo "PRIMER_DNTP_CONC=$DNTP"
    echo "PRIMER_PICK_INTERNAL_OLIGO=0"
    echo "PRIMER_EXPLAIN_FLAG=1"
    if [ -n "$THERMO_PATH" ]; then
      echo "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=$THERMO_PATH"
      # Reasonable dimer/hairpin thresholds (thermo-based)
      echo "PRIMER_MAX_SELF_ANY_TH=45.0"
      echo "PRIMER_MAX_SELF_END_TH=35.0"
      echo "PRIMER_PAIR_MAX_COMPL_ANY_TH=45.0"
      echo "PRIMER_PAIR_MAX_COMPL_END_TH=35.0"
    fi
    echo "="
  } > "$P3IN"

  primer3_core -output="$P3OUT" "$P3IN" >/dev/null

  # Parse Boulder-IO to TSV rows
  awk -v id="$MARKER_ID" -v seq="$MARKER_SEQ" -v want="$NUM_RETURN" '
  BEGIN{ FS="=" }
  {
    key=$1; val=$2
    gsub(/\r$/,"",val)
    if (key ~ /^PRIMER_LEFT_[0-9]+_SEQUENCE$/)       { idx=key; sub(/^PRIMER_LEFT_/,"",idx); sub(/_SEQUENCE$/,"",idx); left_seq[idx]=val }
    else if (key ~ /^PRIMER_RIGHT_[0-9]+_SEQUENCE$/) { idx=key; sub(/^PRIMER_RIGHT_/,"",idx); sub(/_SEQUENCE$/,"",idx); right_seq[idx]=val }
    else if (key ~ /^PRIMER_LEFT_[0-9]+$/)           { idx=key; sub(/^PRIMER_LEFT_/,"",idx);        split(val,a,","); left_start[idx]=a[1]; left_len[idx]=a[2] }
    else if (key ~ /^PRIMER_RIGHT_[0-9]+$/)          { idx=key; sub(/^PRIMER_RIGHT_/,"",idx);       split(val,a,","); right_end[idx]=a[1]; right_len[idx]=a[2] }
    else if (key ~ /^PRIMER_LEFT_[0-9]+_TM$/)        { idx=key; sub(/^PRIMER_LEFT_/,"",idx); sub(/_TM$/,"",idx); left_tm[idx]=val }
    else if (key ~ /^PRIMER_RIGHT_[0-9]+_TM$/)       { idx=key; sub(/^PRIMER_RIGHT_/,"",idx); sub(/_TM$/,"",idx); right_tm[idx]=val }
    else if (key ~ /^PRIMER_LEFT_[0-9]+_GC_PERCENT$/){ idx=key; sub(/^PRIMER_LEFT_/,"",idx); sub(/_GC_PERCENT$/,"",idx); left_gc[idx]=val }
    else if (key ~ /^PRIMER_RIGHT_[0-9]+_GC_PERCENT$/){idx=key; sub(/^PRIMER_RIGHT_/,"",idx); sub(/_GC_PERCENT$/,"",idx); right_gc[idx]=val }
    else if (key ~ /^PRIMER_LEFT_[0-9]+_SELF_ANY_TH$/){idx=key; sub(/^PRIMER_LEFT_/,"",idx); sub(/_SELF_ANY_TH$/,"",idx); l_any[idx]=val }
    else if (key ~ /^PRIMER_LEFT_[0-9]+_SELF_END_TH$/){idx=key; sub(/^PRIMER_LEFT_/,"",idx); sub(/_SELF_END_TH$/,"",idx); l_end[idx]=val }
    else if (key ~ /^PRIMER_RIGHT_[0-9]+_SELF_ANY_TH$/){idx=key; sub(/^PRIMER_RIGHT_/,"",idx); sub(/_SELF_ANY_TH$/,"",idx); r_any[idx]=val }
    else if (key ~ /^PRIMER_RIGHT_[0-9]+_SELF_END_TH$/){idx=key; sub(/^PRIMER_RIGHT_/,"",idx); sub(/_SELF_END_TH$/,"",idx); r_end[idx]=val }
    else if (key ~ /^PRIMER_PAIR_[0-9]+_COMPL_ANY_TH$/){ idx=key; sub(/^PRIMER_PAIR_/,"",idx); sub(/_COMPL_ANY_TH$/,"",idx); p_any[idx]=val }
    else if (key ~ /^PRIMER_PAIR_[0-9]+_COMPL_END_TH$/){ idx=key; sub(/^PRIMER_PAIR_/,"",idx); sub(/_COMPL_END_TH$/,"",idx); p_end[idx]=val }
    else if (key ~ /^PRIMER_PAIR_[0-9]+_PRODUCT_SIZE$/){ idx=key; sub(/^PRIMER_PAIR_/,"",idx); sub(/_PRODUCT_SIZE$/,"",idx); p_size[idx]=val }
    else if (key ~ /^PRIMER_PAIR_[0-9]+_PENALTY$/)     { idx=key; sub(/^PRIMER_PAIR_/,"",idx); sub(/_PENALTY$/,"",idx); p_pen[idx]=val }
  }
  END{
    for (i=0; i<want; i++) {
      if ((i in left_seq) && (i in right_seq) && (i in p_size)) {
        # Precompute values
        ls=left_seq[i];  rs=right_seq[i];  ps=p_size[i]
        ll=(i in left_len ? left_len[i] : "");   rl=(i in right_len ? right_len[i] : "")
        ltm=(i in left_tm ? left_tm[i] : "");    rtm=(i in right_tm ? right_tm[i] : "")
        lgc=(i in left_gc ? left_gc[i] : "");    rgc=(i in right_gc ? right_gc[i] : "")
        lany=(i in l_any ? l_any[i] : "");       lend=(i in l_end ? l_end[i] : "")
        rany=(i in r_any ? r_any[i] : "");       rend=(i in r_end ? r_end[i] : "")
        pany=(i in p_any ? p_any[i] : "");       pend=(i in p_end ? p_end[i] : "")
        ppen=(i in p_pen ? p_pen[i] : "")
        # Derive product sequence from 0-based coords: left_start (first base), right_end (last base)
        # substr() is 1-based in awk, so +1 on the start, and add +1 to length for inclusive end
        prseq=""
        if ((i in left_start) && (i in right_end) && right_end[i] >= left_start[i]) {
          prseq = substr(seq, left_start[i] + 1, right_end[i] - left_start[i] + 1)
        }
        printf("%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
               id, seq, prseq, i, ls, rs, ps, ll, rl, ltm, rtm, lgc, rgc, lany, lend, rany, rend, pany, pend, ppen)
      }
    }
  }
' "$P3OUT" >> "$OUT_TSV"

  # Count how many rows we just wrote for this marker (no non-zero exit codes)
  COUNT=$(awk -v id="$MARKER_ID" -F'\t' '($1==id){n++} END{printf("%d", (n? n:0))}' "$OUT_TSV")

  if [ "$COUNT" = "0" ]; then
    # Surface Primer3 reasons (requires PRIMER_EXPLAIN_FLAG=1, which you already set)
    EXPLAIN_PAIR=$(grep -m1 '^PRIMER_PAIR_EXPLAIN='  "$P3OUT" | cut -d= -f2- || true)
    EXPLAIN_LEFT=$(grep -m1 '^PRIMER_LEFT_EXPLAIN='  "$P3OUT" | cut -d= -f2- || true)
    EXPLAIN_RIGHT=$(grep -m1 '^PRIMER_RIGHT_EXPLAIN=' "$P3OUT" | cut -d= -f2- || true)
    echo "⚠ No primer pairs for: $MARKER_ID" >&2
    [ -n "$EXPLAIN_PAIR"  ] && echo "   • PAIRS : $EXPLAIN_PAIR"  >&2
    [ -n "$EXPLAIN_LEFT"  ] && echo "   • LEFT  : $EXPLAIN_LEFT"  >&2
    [ -n "$EXPLAIN_RIGHT" ] && echo "   • RIGHT : $EXPLAIN_RIGHT" >&2
  fi

  echo "▶ Finished $MARKER_ID"
done

echo "▶ Done. Wrote TSV → $OUT_TSV"

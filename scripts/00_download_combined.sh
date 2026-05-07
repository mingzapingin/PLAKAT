#!/usr/bin/env bash
set -euo pipefail

force_dl="$1"; shift        # 1 = ALWAYS download, 0 = use cache if present          
outdir="$1"; shift          # e.g. data/genomes
lists=("$@")

mkdir -p "$outdir"
cache_dir="$outdir/cache"          # ← 1️⃣  new cache folder
mkdir -p "$cache_dir"

for list in "${lists[@]}"; do
    base=$(basename "$list" .txt)
    combined="$outdir/${base}_combined.fasta"
    echo "▶ Building $combined"
    : > "$combined"

    # read the list
    while IFS= read -r line; do
        # skip blank or commented lines
        [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue

        # strip trailing " |REF" if present
        line="${line%|REF}"
        line="$(echo "$line" | sed 's/[[:space:]]*$//')"

        # ── parse Option-B: role:ACC@strain  OR bare ACC ─────────────
        seqacc=""
        if [[ "$line" =~ ^(chr|plasmid):([A-Za-z0-9._]+) ]]; then
            seqacc="${BASH_REMATCH[2]}"
        else
            # fallback: take the first token up to space or @
            seqacc="${line%%[@[:space:]]*}"
        fi

        # sanity check
        if [[ -z "$seqacc" ]]; then
            echo "  ⚠️  Cannot parse accession from: $line" >&2
            continue
        fi

        # use parsed accession for cache filename (safe)
        single="$cache_dir/${seqacc}.fasta"

        if [[ "$force_dl" == "1" || ! -s "$single" ]]; then
            echo "  • $seqacc  (downloading …)"
            if ! efetch -db nuccore -id "$seqacc" -format fasta > "$single".tmp 2>/dev/null; then
                echo "    ❌ efetch failed for $seqacc" >&2
                rm -f "$single".tmp
                continue
            fi
            # basic integrity check: FASTA header present
            if ! grep -q '^>' "$single".tmp; then
                echo "    ❌ no FASTA header for $seqacc (skipping)" >&2
                rm -f "$single".tmp
                continue
            fi
            mv "$single".tmp "$single"
        else
            echo "  • $seqacc  (cached)"
        fi

        cat "$single" >> "$combined"
    done < "$list"
done
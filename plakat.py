#!/usr/bin/env python3
"""
PLAKAT  ─  Precised Local Alignment for marKers Analysis Tool
Master controller script.

First milestone:
    • Parse user-supplied “accession-list” text files
    • Build three in-memory lists:
        1. positive_genomes   (all M. marinum strains except the reference)
        2. reference_genome   (single M. marinum reference strain)
        3. negative_genomes   (non-target species: M. ulcerans, M. chelonae, …)

Future milestones will call the modular bash scripts for each pipeline step.
"""

import argparse
import pathlib
import sys
from typing import List, Dict
import subprocess, shutil
import os, stat 
import datetime, shlex, re
import json

# ────────────────────────────────────────────────────────────────────────────────
# Helpers
# ────────────────────────────────────────────────────────────────────────────────
def count_chr_plasmid(list_path: pathlib.Path):
    lines = [ln.strip() for ln in open(list_path)
             if ln.strip() and not ln.lstrip().startswith("#")]
    # strip trailing '|REF' if present
    lines = [ln.replace("|REF", "").strip() for ln in lines]
    n_chr = sum(ln.startswith("chr:") for ln in lines)
    n_plm = sum(ln.startswith("plasmid:") for ln in lines)
    return n_chr, n_plm

def open_run_log():
    run_id = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    log_dir = pathlib.Path("results/logs") / run_id
    log_dir.mkdir(parents=True, exist_ok=True)
    return log_dir, (log_dir / "plakat.log").open("a", encoding="utf-8")

def run_logged(cmd, logfh, echo_all=False):
    """Stream stdout/stderr to log; optionally filter what appears in console."""
    # Lines shown in console when echo_all=False
    echo_pattern = r'^(▶|•|❌|⚠️|SUMMARY|Step|Finished)'
    pat = None if echo_all else re.compile(echo_pattern)

    print(" $", " ".join(map(shlex.quote, cmd)))
    logfh.write("\n$ " + " ".join(map(shlex.quote, cmd)) + "\n"); logfh.flush()

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in p.stdout:
        logfh.write(line)
        if pat is None or pat.match(line):
            print(line, end="")
    p.wait(); logfh.flush()
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd)

# ── Accession normalization & list readers for params ──────────────────────────
def _norm_acc_token(s: str, strip_version: bool = True) -> str:
    # remove "|REF" marker and any "chr:" / "plasmid:" prefixes
    s = s.strip().replace("|REF", "")
    if ":" in s:
        # e.g., "chr:NZ_CP012345.1" or "plasmid:NZ_ABCD0001.1"
        s = s.split(":", 1)[1].strip()
    # strip db|ACC| wrappers
    m = re.match(r"^[A-Za-z]+\|([^|]+)\|?$", s)
    if m:
        s = m.group(1)
    # optionally strip trailing .version
    if strip_version and "." in s:
        head, tail = s.rsplit(".", 1)
        if tail.isdigit() and len(tail) <= 2:
            s = head
    return s

# ── Accession normalization & list readers for params ──────────────────────────
def _norm_acc_token(s: str, strip_version: bool = True) -> str:
    # remove "|REF" marker and any "chr:" / "plasmid:" prefixes
    s = s.strip().replace("|REF", "")
    if ":" in s:
        # e.g., "chr:NZ_CP012345.1" or "plasmid:NZ_ABCD0001.1"
        s = s.split(":", 1)[1].strip()
    # strip db|ACC| wrappers
    m = re.match(r"^[A-Za-z]+\|([^|]+)\|?$", s)
    if m:
        s = m.group(1)
    # optionally strip trailing .version
    if strip_version and "." in s:
        head, tail = s.rsplit(".", 1)
        if tail.isdigit() and len(tail) <= 2:
            s = head
    return s

def _read_accessions_from_txt(list_path: pathlib.Path, strip_version: bool = True) -> list[str]:
    accs = []
    with open(list_path, "r", encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.lstrip().startswith("#"):
                continue
            
            ln = ln.split("@", 1)[0].strip()
            accs.append(_norm_acc_token(ln, strip_version=strip_version))
    return accs

def _collect_accessions_from_lists(paths: list[pathlib.Path], strip_version: bool = True) -> list[str]:
    s = set()
    for p in paths:
        for acc in _read_accessions_from_txt(p, strip_version=strip_version):
            if acc:
                s.add(acc)
    return sorted(s)

def _norm_tag_from_list(path_like) -> str:
    """
    Derive a stable 'POS tag' from a list filename (e.g., 'M.marinum_pos.txt' -> 'M.marinum_pos').
    """
    stem = pathlib.Path(path_like).stem
    # keep letters, digits, dot, underscore, hyphen; replace others with underscore
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem)

def _params_init(params_dir: pathlib.Path, pos_tag: str) -> pathlib.Path:
    """
    Create an empty params file for this run (JSON text, valid YAML).
    Returns the path to the file.
    """
    params_dir.mkdir(parents=True, exist_ok=True)
    out = params_dir / f"{pos_tag}.yml"  # JSON is valid YAML
    skeleton = {
        "created_at": datetime.datetime.now().isoformat(),
        "pos_tag": pos_tag,
        # you will fill these gradually across steps:
        "genome_mbp": None,
        "pos_accessions": [],
        "neg_accessions_big_run": [],
        "pident_min": None,
        "qcov_min": None,
        "win_bp": None,
        "step_bp": None,
        "mash_matrix": None,
        "manifests_dir": None,
        "blast_inputs": {},          # e.g., {unique_vs_neg, centroids_vs_pos, centroids_vs_neg}
        "normalization": {
            "strip_accession_version": True,
            "strip_db_prefix": True
        },
    }
    # write JSON (YAML-compatible)
    out.write_text(json.dumps(skeleton, indent=2) + "\n", encoding="utf-8")
    return out

def _params_update(params_path: pathlib.Path, **kvs) -> None:
    """
    Update the params file. If a provided value is a dict and the existing
    value is also a dict, perform a shallow merge; else overwrite.
    """
    import json
    try:
        current = json.loads(params_path.read_text(encoding="utf-8"))
    except Exception:
        current = {}

    for k, v in kvs.items():
        if isinstance(v, dict) and isinstance(current.get(k), dict):
            current[k].update(v)  # shallow merge into existing map
        else:
            current[k] = v        # overwrite

    tmp = params_path.with_suffix(".yml.tmp")
    tmp.write_text(json.dumps(current, indent=2) + "\n", encoding="utf-8")
    tmp.replace(params_path)


def _compute_genome_mbp_from_fasta(fa: pathlib.Path) -> float:
    # Use seqkit if available (fast), else pure-Python fallback
    try:
        import shutil, subprocess
        if shutil.which("seqkit"):
            cmd = ["seqkit", "fx2tab", "-nl", str(fa)]
            out = subprocess.check_output(cmd, text=True)
            total_bp = 0
            for ln in out.splitlines():
                parts = ln.strip().split("\t")
                if len(parts) >= 2:
                    total_bp += int(parts[1])
            return round(total_bp / 1e6, 6)
    except Exception:
        pass
    # fallback: python parser
    total = 0
    with open(fa, "r", encoding="utf-8") as fh:
        seq = []
        for ln in fh:
            if ln.startswith(">"):
                if seq:
                    total += sum(len(s.strip()) for s in seq)
                    seq = []
            else:
                seq.append(ln)
        if seq:
            total += sum(len(s.strip()) for s in seq)
    return round(total / 1e6, 6)

# ────────────────────────────────────────────────────────────────────────────────
# CLI
# ────────────────────────────────────────────────────────────────────────────────
def make_cli() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="plakat",
        description="PLAKAT – Precised Local Alignment for marKers Analysis Tool"
    )

    p.add_argument(
        "--positive-list", "-p",
        required=True,
        type=pathlib.Path,
        help="Text file of accession numbers for the TARGET species (positive set)"
    )
    p.add_argument(
        "--reference-list", "-r",
        required=True,
        type=pathlib.Path,
        help="Text file (usually 1 entry) for the REFERENCE strain used for coordinate mapping"
    )
    p.add_argument(
        "--negative-list", "-n",
        required=True,
        nargs="+",
        type=pathlib.Path,
        help="One or more text files of accession numbers for NON-target species"
    )
    p.add_argument(
        "--dry",
        action="store_true",
        help="Parse input lists and exit (no workflow steps)"
    )
    p.add_argument("--force-download", "-F",
        action="store_true",
        help="Always re-download each accession (ignore FASTA cache)")
    p.add_argument("--rebuild-db",
        action="store_true",
        help="Force BLAST DB rebuild but keep FASTA cache")
    p.add_argument(
    "--verbose", "-v",
    action="store_true",
    help="Show all tool stdout in console (everything is always saved to the log)."
    )
    p.add_argument(
        "--params",
        action="store_true",
        help="If set, create and maintain a run params file under results/params/"
    )
    p.add_argument(
        "--params-dir",
        type=pathlib.Path,
        default=pathlib.Path("results/params"),
        help="Directory to store the params file (default: results/params)"
    )
    return p


def main(argv=None) -> None:
    cli = make_cli()
    args = cli.parse_args(argv)

    log_dir, logfh = open_run_log()
    logfh.write(f"# PLAKAT run {datetime.datetime.now().isoformat()}\n"); logfh.flush()

    # ── Validate files exist ───────────────────────────────────────────────
    for path in [args.positive_list, args.reference_list, *args.negative_list]:
        if not path.is_file():
            cli.error(f"▶ List file not found: {path}")

    # ── Echo a quick summary ───────────────────────────────────────────────
    def count(txt): return sum(1 for _ in open(txt) if _.strip() and not _.startswith("#"))
    print(f"\nSUMMARY")
    print(f"▶  Positive : {args.positive_list.name}  ({count(args.positive_list)} accessions)")
    print(f"▶  Reference: {args.reference_list.name}  ({count(args.reference_list)} accessions)")
    for neg in args.negative_list:
        print(f"▶  Negative : {neg.name:<20} ({count(neg)} accessions)")
    print()

    # ── Normalize lists so you can reuse the same file for pos+ref ────────────
    lists_tmp = pathlib.Path("data/lists_tmp")
    lists_tmp.mkdir(parents=True, exist_ok=True)

    def _read_lines(p: pathlib.Path):
        return [ln.strip() for ln in open(p) if ln.strip() and not ln.lstrip().startswith("#")]

    # Build REF-only list (keep both chr and plasmid for now)
    ref_src_lines = _read_lines(args.reference_list)
    ref_only = [ln for ln in ref_src_lines if "|REF" in ln]
    if not ref_only:
        # fallback: pick chr: lines if present, else first non-comment line
        ref_only = [ln for ln in ref_src_lines if ln.startswith("chr:")] or ref_src_lines[:1]

    # Build POS list excluding REF entries (avoid double-counting the reference)
    pos_src_lines = _read_lines(args.positive_list)
    pos_no_ref = [ln for ln in pos_src_lines if "|REF" not in ln]

    # Write effective lists with controlled basenames
    eff_pos = lists_tmp / args.positive_list.name              # keep same basename → downstream names stable
    eff_ref = lists_tmp / f"{args.reference_list.stem}_ref.txt"  # distinct basename for REF

    eff_pos.write_text("\n".join(pos_no_ref) + ("\n" if pos_no_ref else ""), encoding="utf-8")
    eff_ref.write_text("\n".join(ref_only) + ("\n" if ref_only else ""), encoding="utf-8")

    # (Optional) print a quick delta so user sees what’s happening
    print(f"▶  Using positives from: {eff_pos}  ({len(pos_no_ref)}/{len(pos_src_lines)} entries after removing |REF)")
    print(f"▶  Using reference from: {eff_ref}  ({len(ref_only)} REF entries)")

    num_chr, num_plasmid = count_chr_plasmid(eff_ref)
    print(f"▶ Reference has {num_chr} chromosome(s) and {num_plasmid} plasmid(s)")

    params_path = None
    if args.params:
        pos_tag = _norm_tag_from_list(args.positive_list)  # derive tag from POS list filename
        params_path = _params_init(args.params_dir, pos_tag)
        print(f"▶ Params file initialized at: {params_path}")
        # seed a few known paths now if you want:
        _params_update(params_path,
            manifests_dir=str(pathlib.Path(args.positive_list).parent.resolve()),
        )
    
    if args.dry:
        sys.exit(0)

    # ── Step 00: download one combined FASTA per list ──────────────────────
    scripts_dir   = pathlib.Path(__file__).parent / "scripts"

    # ── ensure every bash module is executable ─────────────────────────────
    for sh in scripts_dir.glob("*.sh"):
        mode = os.stat(sh).st_mode
        os.chmod(sh, mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

    download_scr  = scripts_dir / "00_download_combined.sh"      
    outdir        = pathlib.Path("data/genomes")
    outdir.mkdir(parents=True, exist_ok=True)

    force_dl = "1" if args.force_download else "0"   # 0 = download if missing
    force_db  = "1" if args.rebuild_db else "0"      # 0 = build db if missing

    cmd = [
    str(download_scr),
    force_dl,
    str(outdir),
    str(eff_pos),         # ← positives without |REF
    str(eff_ref),         # ← REF-only (chr+plasmid for now)
    *map(str, args.negative_list),
    ]   

    print("▶ Running step 0 [Download]:\n   ", " ".join(cmd), "\n")
    run_logged(cmd, logfh, echo_all=args.verbose)

    if params_path:
        # include REF in the POS universe (useful for POS_presence_rate later)
        pos_accessions = _collect_accessions_from_lists([eff_pos, eff_ref], strip_version=True)
        # NEG universe across all provided NEG list files
        neg_accessions = _collect_accessions_from_lists(list(args.negative_list), strip_version=True)

        _params_update(params_path,
            pos_accessions=pos_accessions,
            neg_accessions_big_run=neg_accessions,
            # optionally record the list file paths you actually used
            lists={
                "pos": str(eff_pos.resolve()),
                "ref": str(eff_ref.resolve()),
                "neg": [str(p.resolve()) for p in args.negative_list],
            },
        )
    
    # Build chr IDs file (as you already have)
    ref_chr_ids = lists_tmp / f"{args.reference_list.stem}_ref_chr_ids.txt"
    with open(ref_chr_ids, "w", encoding="utf-8") as fh:
        for ln in ref_only:
            ln = ln.replace("|REF", "").strip()
            if ln.startswith("chr:"):
                acc = ln.split(":", 1)[1].split("@", 1)[0]
                fh.write(acc + "\n")

    # ── Step 01: sliding-window slicing ──────────────────────────────────
    slice_script = scripts_dir / "01_slice_genome.sh"
    sliced_dir   = pathlib.Path("results/sliced")
    sliced_dir.mkdir(parents=True, exist_ok=True)

    #Use the reference combined FASTA (matches eff_ref basename)
    ref_stem     = eff_ref.stem                           # e.g. "M.marinum_ref"
    ref_combined = outdir / f"{ref_stem}_combined.fasta"  # data/genomes/M.marinum_ref_combined.fasta

    size = 500
    overlap = 250

    cmd_slice = [
        str(slice_script),
        str(sliced_dir),
        str(ref_combined),    # ← slice the REF combined
        str(size),
        str(overlap),
        "1",                  # chrom_only safety
        str(ref_chr_ids),     # explicit chr accessions
    ]
    print("▶ Running slicing step 1:\n   ", " ".join(cmd_slice), "\n")
    run_logged(cmd_slice, logfh, echo_all=args.verbose)

    if params_path:
        mbp = _compute_genome_mbp_from_fasta(ref_combined)
        _params_update(params_path, genome_mbp=mbp,win_bp=size,step_bp=overlap)

    # reference-based frags path
    frag_fasta = sliced_dir / f"{ref_combined.stem}_frags_{size}bp.fasta"

    # ── Step 02: negative screen + unique FASTA (02+03 combined) ─────────
    neg_script = scripts_dir / "02_neg_screen_unique.sh"
    unique_fasta = pathlib.Path("results/filtered") / f"{ref_stem}_unique_{size}bp.fasta"
    unique_fasta.parent.mkdir(parents=True, exist_ok=True)
    blast_unique_out = pathlib.Path("results/blast") / f"{ref_stem}_unique_vs_neg.tsv"
    blast_unique_out.parent.mkdir(parents=True, exist_ok=True)

    # build a temporary “negative fasta” = concat all negatives into one
    neg_fasta_tmp = outdir / "NEGATIVE_TEMP.fasta"
    with open(neg_fasta_tmp, "w") as handle:
        for neg in args.negative_list:
            stem = neg.stem
            fasta = outdir / f"{stem}_combined.fasta"
            handle.write(open(fasta).read())

    cmd_neg = [
        str(neg_script),
        force_db,               # "1" = rebuild DB, "0" = reuse
        str(frag_fasta),        # windows from Step 01 (reference-based)
        str(neg_fasta_tmp),     # concatenated negatives
        str(unique_fasta),      # output FASTA of no-hit (unique) windows
        str(blast_unique_out),
        "4",                    # threads
    ]

    print("▶ Running negative screen & unique selection (step 2):\n   ", " ".join(cmd_neg), "\n")
    run_logged(cmd_neg, logfh, echo_all=args.verbose)

    if params_path:
        _params_update(params_path,
            blast_inputs={"unique_vs_neg": str(blast_unique_out.resolve())}
    )

    # ── Step 03: deduplicate + cluster  ───────────────────────────
    dedup_script = scripts_dir / "03_dedup_cluster.sh"
    dedup_dir = pathlib.Path("results/dedup")
    dedup_dir.mkdir(parents=True, exist_ok=True)

    cmd_dedup = [
        str(dedup_script),
        str(unique_fasta),   # ← output from combined step
        str(dedup_dir),
        "4"
    ]
    print("▶ Running de-duplication & clustering (step 3):\n   ", " ".join(cmd_dedup), "\n")
    run_logged(cmd_dedup, logfh, echo_all=args.verbose)
    
    # Step 4 — BLAST centroids vs POSITIVES (not negatives)
    blast_pos_script = scripts_dir / "04_blast_vs_positive.sh"

    pos_stem     = args.positive_list.stem
    pos_combined = outdir / f"{pos_stem}_combined.fasta"   # built in Step 00 from eff_pos
    unique_stem  = unique_fasta.stem                       # e.g., M.marinum_ref_unique_500bp
    clustered_fa = dedup_dir / f"{unique_stem}_95.fasta"

    blast_pos_out = pathlib.Path("results/blast") / f"{pos_stem}_centroids_vs_pos.tsv"
    blast_pos_out.parent.mkdir(parents=True, exist_ok=True)

    cmd_blast_pos = [
        str(blast_pos_script),
        force_db,                # "1" = rebuild DB, "0" = reuse
        str(clustered_fa),       # queries = CD-HIT centroids
        str(pos_combined),       # DB = positives (no REF)
        str(blast_pos_out),      # output TSV (must end .tsv)
        "4",                     # threads
        "megablast",             # (optional) faster for same-species
    ]
    print("▶ Running centroid→positive BLAST (step 4):\n   ", " ".join(map(str, cmd_blast_pos)), "\n")
    run_logged(cmd_blast_pos, logfh, echo_all=args.verbose)
    
    # Build POSITIVE chromosome sseqid list (accessions only)
    pos_chr_ids = lists_tmp / f"{args.positive_list.stem}_pos_chr_ids.txt"
    with open(pos_chr_ids, "w", encoding="utf-8") as fh:
        for ln in pos_no_ref:
            ln = ln.replace("|REF", "").strip()
            if ln.startswith("chr:"):
                acc = ln.split(":", 1)[1].split("@", 1)[0]
                fh.write(acc + "\n")

    # ── Step 05: select core single-copy markers across ALL positive chromosomes ──
    core_script = scripts_dir / "05_select_core_singlecopy.sh"
    core_out    = pathlib.Path("results/core") / f"{pos_stem}_core_singlecopy.fasta"
    core_out.parent.mkdir(parents=True, exist_ok=True)

    cmd_core = [
        str(core_script),
        str(blast_pos_out),   # 1: TSV (centroids vs positives)
        str(clustered_fa),    # 2: centroids FASTA (dedup + 95%)
        str(core_out),        # 3: OUTPUT core markers FASTA   
        str(pos_chr_ids),     # 4: chromosome sseqid list (one per line)
        "98",                 # 5: min pident (optional)
        "80",                 # 6: min qcov (optional)
    ]

    print("▶ Selecting core single-copy markers (step 5):\n   ", " ".join(map(str, cmd_core)), "\n")
    run_logged(cmd_core, logfh, echo_all=args.verbose)
    
    if params_path:
        _params_update(params_path,
            blast_inputs={"centroids_vs_pos": str(blast_pos_out.resolve())}
    )
        
    # ── Step 06: Loose BLAST screen vs negatives ───────────────────────────────
    loose_script   = scripts_dir / "06_loose_screen.sh"
    final_markers  = pathlib.Path("results/core") / f"{pos_stem}_final_markers.fasta"
    final_markers.parent.mkdir(parents=True, exist_ok=True)
    blast_neg_out = pathlib.Path("results/blast") / f"{pos_stem}_centroids_vs_neg.tsv"
    pident_min = 80
    qcov_min = 50
    cmd_loose = [
        str(loose_script),
        str(core_out),        # from Step 5 (core single-copy)
        str(neg_fasta_tmp),   # the same NEGATIVE_TEMP.fasta built earlier
        str(final_markers),
        str(blast_neg_out),
        "4",                   # threads
        str(pident_min),
        str(qcov_min)
    ]

    print("▶ Running loose BLAST screen (step 6):\n   ", " ".join(map(str, cmd_loose)), "\n")
    run_logged(cmd_loose, logfh, echo_all=args.verbose)
    
    if params_path:
        _params_update(params_path,
            pident_min = pident_min,
            qcov_min = qcov_min,
            blast_inputs={"centroids_vs_neg": str(blast_neg_out.resolve())}
    )
        
    # ── Step 07: Prepare reference (download .gbk ➜ clean GFF, concat FASTA) ──
    prep_ref_script = scripts_dir / "07_prepare_reference.sh"   # script file unchanged
    ref_dir         = pathlib.Path("data/ref")
    ref_dir.mkdir(parents=True, exist_ok=True)

    cmd_prep = [
        str(prep_ref_script),
        str(eff_ref),     # ← use the REF-only temp list you wrote earlier
        str(ref_dir),
        str(outdir)
    ]
    print("▶ Preparing reference (step 7):\n   ", " ".join(map(str, cmd_prep)), "\n")
    run_logged(cmd_prep, logfh, echo_all=args.verbose)

    ref_combined = pathlib.Path("data/ref/fasta") / "reference_combined.fasta"
    
    # ── Step 08: Map final markers to reference ─────────────────────────────
    map_script   = scripts_dir / "08_blast_vs_ref.sh"   # ← was 09_
    marker_bed   = pathlib.Path("results/annotation") / f"{pos_stem}_markers_on_ref.bed"
    marker_bed.parent.mkdir(parents=True, exist_ok=True)

    cmd_map = [
        str(map_script),
        str(final_markers),   # from step 7
        str(ref_combined),    # make sure this matches the file your Step 7 wrote
        str(marker_bed),
        "4"                   # threads
    ]

    print("▶ Mapping markers to reference (step 8):\n   ", " ".join(cmd_map), "\n")
    run_logged(cmd_map, logfh, echo_all=args.verbose)
    
    # ── Step 09: Annotate marker BED with GFF genes ─────────────────────────
    annotate_script = scripts_dir / "09_annotate_with_gff.sh"  # ← was 10_
    gff_dir         = pathlib.Path("data/ref/gff")
    anno_out        = pathlib.Path("results/annotation") / f"{pos_stem}_marker_annotations.tsv"
    anno_out.parent.mkdir(parents=True, exist_ok=True)

    cmd_anno = [
        str(annotate_script),
        str(marker_bed),   # from step 8
        str(gff_dir),
        str(anno_out),
        # optional: min frac overlap of marker, min bp overlap
        # e.g. "0.20", "100"
    ]
    print("▶ Annotating markers with reference GFF (step 9):\n   ", " ".join(cmd_anno), "\n")
    run_logged(cmd_anno, logfh, echo_all=args.verbose)
    
    # ── Step 10: collapse to gene summary ───────────────────────────────────
    summ_script  = scripts_dir / "10_summarize_genes.sh"
    anno_in      = pathlib.Path("results/annotation") / f"{pos_stem}_marker_annotations.tsv"
    summary_out  = pathlib.Path("results/annotation") / f"{pos_stem}_marker_summary.tsv"
    summary_out.parent.mkdir(parents=True, exist_ok=True)

    cmd_summ = [
        str(summ_script),
        str(anno_in),
        str(summary_out)
    ]
    print("▶ Running gene-summary step (10):\n   ", " ".join(cmd_summ), "\n")
    run_logged(cmd_summ, logfh, echo_all=args.verbose)

    # ── Step 11: generate primer pairs from final markers ───────────────────────────────────
    primer_script  = scripts_dir / "11_design_primers_from_markers.sh"
    primer_out  = pathlib.Path("results/primer") / f"{pos_stem}_marker_summary.tsv"
    primer_out.parent.mkdir(parents=True, exist_ok=True)

    cmd_prim = [
        str(primer_script),
        str(final_markers),
        str(primer_out)
    ]
    print("▶ Running generate primer step (11):\n   ", " ".join(cmd_summ), "\n")
    run_logged(cmd_prim, logfh, echo_all=args.verbose)

    print("▶ Finish !!")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
import argparse, json, re, shutil, subprocess, sys
from pathlib import Path
from typing import Dict, List, Tuple

# ───────────────── helpers ─────────────────

def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def run(cmd: List[str]) -> str:
    try:
        out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
        return out.decode("utf-8", "replace")
    except subprocess.CalledProcessError as e:
        print(f"❌ Command failed: {' '.join(cmd)}\n{e.output.decode('utf-8','replace')}", file=sys.stderr)
        sys.exit(1)

def datasets_jsonl(taxon: str, report: str, levels: str, source: str, exact: bool) -> List[dict]:
    """
    Call NCBI datasets and return a list of JSON objects (one per line).
    """
    if not have("datasets"):
        print("❌ 'datasets' CLI not found. Install ncbi-datasets-cli.", file=sys.stderr)
        sys.exit(1)

    cmd = ["datasets", "summary", "genome", "taxon", str(taxon),
           "--assembly-version", "latest", "--as-json-lines"]
    if levels:
        cmd += ["--assembly-level", levels]
    if report == "sequence":
        cmd += ["--report", "sequence"]
    if source in ("refseq", "genbank"):
        cmd += ["--assembly-source", source]
    if exact:
        cmd += ["--tax-exact-match"]

    txt = run(cmd).strip()
    if not txt:
        return []
    out = []
    for line in txt.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            out.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return out

def sanitize(s: str) -> str:
    s = (s or "").strip()
    s = re.sub(r"[ \t/;:|]+", "_", s)
    return s or "NA"

def infer_strain_from_assembly_obj(obj: dict) -> str:
    org = obj.get("organism") or {}
    inf = org.get("infraspecific_names") or {}
    strain = inf.get("strain") or inf.get("isolate") or ""
    if not strain:
        strain = obj.get("strain") or obj.get("isolate") or ""
    return strain

def get_assembly_accession(obj: dict) -> str:
    return (obj.get("accession")
            or obj.get("assembly_accession")
            or (obj.get("assembly", {}) or {}).get("assembly_accession", ""))

def get_organism_name(obj: dict) -> str:
    org = obj.get("organism") or {}
    return org.get("organism_name") or obj.get("organism_name") or ""

def get_assembly_source(obj: dict) -> str:
    return (obj.get("assembly_info") or {}).get("assembly_source") or obj.get("assembly_source") or ""

def get_refseq_category(obj: dict) -> str:
    return ((obj.get("assembly_info") or {}).get("refseq_category")
            or obj.get("refseq_category") or "").lower()

def choose_seq_acc(seq_obj: dict, prefer_refseq: bool) -> Tuple[str, str]:
    # accept snake_case and camelCase
    ref = (seq_obj.get("refseq_accession") or seq_obj.get("refseqAccession")
           or seq_obj.get("refseq_acc") or "")
    gen = (seq_obj.get("genbank_accession") or seq_obj.get("genbankAccession")
           or seq_obj.get("genbank_acc") or "")
    if prefer_refseq and ref:
        return ref, "refseq"
    if gen:
        return gen, "genbank"
    if ref:
        return ref, "refseq"
    return "", ""

def detect_role_chr_plasmid(seq_obj: dict) -> str:
    loc = (seq_obj.get("assigned_molecule_location_type")
           or seq_obj.get("assignedMoleculeLocationType") or "")
    role = seq_obj.get("role") or ""
    s = f"{loc} {role}".lower()
    return "plasmid" if "plasmid" in s else "chr"

def extract_seq_name(seq_obj: dict) -> str:
    return (seq_obj.get("sequence_name") or seq_obj.get("sequenceName")
            or seq_obj.get("name") or seq_obj.get("chr_name") or "")

def build_from_json(asm_jsonl: List[dict], seq_jsonl: List[dict],
                    prefer_refseq: bool, sanitize_strain_flag: bool
                    ) -> Tuple[List[str], List[List[str]]]:
    """
    Returns:
      optionB_lines: ["chr:ACC@STRAIN", ...]  (not used for writing; kept for debugging)
      manifest_rows: [[assembly, seq_acc, role, strain, organism, source, sequence_name], ...]
    """
    # Map assembly accession -> (strain, organism, source)
    asm_map: Dict[str, Tuple[str,str,str]] = {}
    for obj in asm_jsonl:
        acc = get_assembly_accession(obj)
        if not acc:
            continue
        strain = infer_strain_from_assembly_obj(obj) or "NA"
        if sanitize_strain_flag:
            strain = sanitize(strain)
        org = get_organism_name(obj) or ""
        src = get_assembly_source(obj) or ""
        asm_map[acc] = (strain, org, src)

    optionB_lines: List[str] = []
    manifest_rows: List[List[str]] = []

    for obj in seq_jsonl:
        asm_acc = obj.get("assembly_accession") or obj.get("accession") or ""
        if not asm_acc:
            continue
        seq_list = obj["sequences"] if isinstance(obj.get("sequences"), list) else [obj]
        for s in seq_list:
            seq_acc, src_tag = choose_seq_acc(s, prefer_refseq)
            if not seq_acc:
                continue
            role = detect_role_chr_plasmid(s)
            seqname = extract_seq_name(s)
            strain, organism, asm_src = asm_map.get(asm_acc, ("NA", "", ""))

            optionB_lines.append(f"{role}:{seq_acc}@{strain}")
            manifest_rows.append([asm_acc, seq_acc, role, strain, organism,
                                  (src_tag or asm_src or ""), seqname])
    return optionB_lines, manifest_rows
def _is_rep(o):
    cat = ((o.get("assembly_info") or {}).get("refseq_category") or o.get("refseq_category") or "").lower()
    return cat == "representative genome"

def _is_type_material(o):
    # True if assembly or biosample flags type material
    if (o.get("type_material") or (o.get("assembly_info") or {}).get("type_material")):
        return True
    bs = (o.get("assembly_info") or {}).get("biosample") or o.get("biosample") or {}
    for a in bs.get("attributes", []):
        if str(a.get("name","")).lower() == "type-material":
            return True
    return False

def _assembly_level(o):
    return ((o.get("assembly_info") or {}).get("assembly_level") or "").lower()

def _num_contigs(o):
    try:
        return int((o.get("assembly_stats") or {}).get("number_of_contigs", 10**9))
    except Exception:
        return 10**9

def _total_len(o):
    try:
        return int((o.get("assembly_stats") or {}).get("total_sequence_length", 0))
    except Exception:
        return 0

def _is_refseq_src(o):
    srcdb = (o.get("source_database") or "").upper()
    src   = ((o.get("assembly_info") or {}).get("assembly_source") or "").lower()
    return "REFSEQ" in srcdb or src == "refseq"

def _release_date_key(o):
    d = ((o.get("assembly_info") or {}).get("release_date") or "").replace("Z","")
    # newer is better → use negative timestamp surrogate by string sort fallback
    return d

# ───────────────── main ─────────────────

def main():
    ap = argparse.ArgumentParser(description="Create Option-B species list and manifest.tsv from NCBI Datasets by taxon.")
    ap.add_argument("--taxon", required=True, help="NCBI Taxon ID (e.g., 1781)")
    ap.add_argument("--levels", default="chromosome,complete", help="Assembly levels (comma-separated)")
    ap.add_argument("--source", default="all", choices=["refseq","genbank","all"], help="Assembly source")
    ap.add_argument("--exact-match", action="store_true", help="Exact taxon node only")
    ap.add_argument("--prefer-refseq", type=lambda x: str(x).lower()!="false", default=True)
    ap.add_argument("--sanitize-strain", type=lambda x: str(x).lower()!="false", default=True)

    ap.add_argument("--out-dir", default=".", help="Output directory")
    ap.add_argument("--out", default="species.txt", help="Option-B list filename")
    ap.add_argument("--manifest", default="manifest.tsv", help="Manifest filename")

    # Reference tagging (optional overrides)
    ap.add_argument("--ref-assembly", default="", help="Mark this assembly (e.g., GCF_003391395.1) as |REF")
    ap.add_argument("--ref-strain", default="", help="Mark this strain (sanitized) as |REF")

    args = ap.parse_args()
    prefer_refseq = bool(args.prefer_refseq)
    sanitize_strain_flag = bool(args.sanitize_strain)
    src = args.source if args.source in ("refseq","genbank") else ""

    print(f"▶ Fetching assemblies for taxon {args.taxon} (levels={args.levels or 'ANY'}, source={args.source}, exact={args.exact_match})", file=sys.stderr)
    asm_jsonl = datasets_jsonl(args.taxon, report="assembly",
                               levels=args.levels, source=src, exact=args.exact_match)
    if not asm_jsonl:
        print("⚠️ No assemblies under current level filter; retrying without --assembly-level …", file=sys.stderr)
        asm_jsonl = datasets_jsonl(args.taxon, report="assembly",
                                   levels=None, source=src, exact=args.exact_match)
        if not asm_jsonl:
            print("⚠️ Still no assemblies found.", file=sys.stderr)
            outdir = Path(args.out_dir); outdir.mkdir(parents=True, exist_ok=True)
            (outdir / args.out).write_text("", encoding="utf-8")
            with open(outdir / args.manifest, "w", encoding="utf-8") as fh:
                fh.write("assembly_accession\tseq_accession\trole\tstrain\torganism\tsource\tsequence_name\n")
            sys.exit(0)

    # Identify NCBI RefSeq “reference genome” assemblies
    ref_assemblies = {
        (get_assembly_accession(o)) for o in asm_jsonl
        if get_refseq_category(o) == "reference genome"
    }

    # 1) If no “reference genome”, try “representative genome”
    if not ref_assemblies:
        rep = { get_assembly_accession(o) for o in asm_jsonl if _is_rep(o) }
        if rep:
            ref_assemblies = rep

    # 2) If still empty, pick a single best assembly deterministically
    if not ref_assemblies:
        def _score(o):
            # lower is better
            return (
                0 if _is_refseq_src(o) else 1,
                0 if _is_type_material(o) else 1,
                0 if _assembly_level(o) == "complete genome" else 1,
                _num_contigs(o),
                -_total_len(o),
                # sort by date descending (newer preferred): reverse via minus by string won’t work,
                # so we’ll sort ascending on string then invert later via tuple placement; keep at the end
                _release_date_key(o)
            )
        best = sorted(asm_jsonl, key=_score)[0]
        ref_assemblies = { get_assembly_accession(best) }

    seq_jsonl = datasets_jsonl(args.taxon, report="sequence",
                               levels=args.levels, source=src, exact=args.exact_match)
    if not seq_jsonl:
        print("⚠️ No sequences under current level filter; retrying without --assembly-level …", file=sys.stderr)
        seq_jsonl = datasets_jsonl(args.taxon, report="sequence",
                                   levels=None, source=src, exact=args.exact_match)
        if not seq_jsonl:
            print("⚠️ Still no sequences found.", file=sys.stderr)
            outdir = Path(args.out_dir); outdir.mkdir(parents=True, exist_ok=True)
            (outdir / args.out).write_text("", encoding="utf-8")
            with open(outdir / args.manifest, "w", encoding="utf-8") as fh:
                fh.write("assembly_accession\tseq_accession\trole\tstrain\torganism\tsource\tsequence_name\n")
            sys.exit(0)

    # Build tables
    optionB_lines, manifest_rows = build_from_json(asm_jsonl, seq_jsonl, prefer_refseq, sanitize_strain_flag)

    # Dedup the manifest (set of tuples) and sort for stability
    manifest_rows = sorted({tuple(r) for r in manifest_rows})

    # Compose one .txt line per (role, seq_acc, strain), tagging |REF:
    #  - Tag if assembly is in ref_assemblies
    #  - OR if matches manual overrides (--ref-assembly / --ref-strain)
    best: Dict[Tuple[str,str,str], bool] = {}
    ref_asm_arg = args.ref_assembly.strip()
    ref_str_arg = args.ref_strain.strip()

    for asm_acc, seq_acc, role, strain, organism, source, seqname in manifest_rows:
        key = (role, seq_acc, strain)
        is_ref = (asm_acc in ref_assemblies) or (ref_asm_arg and asm_acc == ref_asm_arg) or (ref_str_arg and strain == ref_str_arg)
        best[key] = best.get(key, False) or is_ref

    lines = [f"{r}:{s}@{t}" + (" |REF" if flag else "") for (r,s,t), flag in sorted(best.items())]

    # Write outputs
    outdir = Path(args.out_dir); outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / args.out
    man_path = outdir / args.manifest

    out_path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")

    with open(man_path, "w", encoding="utf-8") as fh:
        fh.write("assembly_accession\tseq_accession\trole\tstrain\torganism\tsource\tsequence_name\n")
        for r in manifest_rows:
            fh.write("\t".join(r) + "\n")


    print(f"▶ Wrote: {out_path}  ({len(lines)} lines)", file=sys.stderr)
    print(f"▶ Wrote: {man_path}  ({len(manifest_rows)} rows)", file=sys.stderr)

if __name__ == "__main__":
    main()

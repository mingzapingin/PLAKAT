# PLAKAT

**PLAKAT** — **Precised Local Alignment marKers Analysis Tool** — is a genome-guided workflow for discovering short, species-specific DNA marker regions from bacterial genomes.

The pipeline is designed for closely related bacterial clades where traditional marker genes such as 16S rRNA, `rpoB`, `hsp65`, or MLSA panels may not provide enough resolution for species-level or subspecies-level identification.

PLAKAT searches across the whole reference genome, keeps candidate regions that are conserved in the target species, removes regions shared with close relatives, checks for single-copy behavior, maps retained markers back to the reference genome, annotates marker locations, and designs primers for downstream wet-lab validation.

---

## Project goals

PLAKAT aims to:

- discover short marker regions that are specific to a target bacterial species;
- screen candidate markers against close relatives and other non-target genomes;
- retain markers that are conserved in the positive group and absent or divergent in the negative group;
- prioritize single-copy, reference-mappable regions;
- annotate marker positions using reference genome features;
- design PCR/Sanger/qPCR-compatible primers for retained markers.

---

## Repository structure

```text
PLAKAT/
├── plakat.py
├── make_manifest.py
├── scripts/
│   ├── 00_download_combined.sh
│   ├── 01_slice_genome.sh
│   ├── 02_neg_screen_unique.sh
│   ├── 03_dedup_cluster.sh
│   ├── 04_blast_vs_positive.sh
│   ├── 05_select_core_singlecopy.sh
│   ├── 06_loose_screen.sh
│   ├── 07_prepare_reference.sh
│   ├── 08_blast_vs_ref.sh
│   ├── 09_annotate_with_gff.sh
│   ├── 10_summarize_genes.sh
│   └── 11_design_primers_from_markers.sh
├── README.md
└── .gitignore
```

Large downloaded genomes, BLAST databases, intermediate files, logs, and final results are intentionally excluded from this repository.

---

## Main scripts

### `make_manifest.py`

Builds standardized genome accession lists and manifest tables from NCBI taxon IDs using the NCBI Datasets CLI.

It can produce:

- an accession list in PLAKAT input format;
- a `manifest.tsv` file mapping assemblies, sequence accessions, roles, strains, organism names, sources, and sequence names;
- automatic `|REF` tagging for a selected reference assembly.

Example:

```bash
python make_manifest.py \
  --taxon 1781 \
  --levels chromosome,complete \
  --source refseq \
  --out-dir manifests/M.marinum \
  --out M.marinum.txt \
  --manifest M.marinum_manifest.tsv
```

### `plakat.py`

Main controller script for the PLAKAT workflow.

It accepts:

- a positive genome list;
- a reference genome list;
- one or more negative genome lists.

Example:

```bash
python plakat.py \
  --positive-list manifests/M.marinum/M.marinum.txt \
  --reference-list manifests/M.marinum/M.marinum.txt \
  --negative-list manifests/M.ulcerans/M.ulcerans.txt manifests/M.chelonae/M.chelonae.txt \
  --params
```

---

## Input list format

PLAKAT accepts accession lists with one sequence per line.

Recommended format:

```text
chr:NC_000000.1@Strain_name
plasmid:NC_000001.1@Strain_name
chr:NC_000002.1@Reference_strain|REF
```

Notes:

- `chr:` marks chromosome sequences.
- `plasmid:` marks plasmid sequences.
- `@Strain_name` stores the strain label.
- `|REF` marks the reference sequence used for coordinate mapping.
- Blank lines and lines beginning with `#` are ignored.

The same species list can be used as both the positive list and reference list if the reference entry is marked with `|REF`.

---

## Workflow summary

PLAKAT currently performs the following major steps:

1. **Download combined FASTAs**

   Downloads or reuses cached FASTA files for positive, reference, and negative accessions.

2. **Slice reference genome**

   Splits the reference genome into overlapping windows, typically 500 bp windows with 250 bp step size.

3. **Screen against negative genomes**

   Removes reference windows that show similarity to non-target genomes.

4. **Deduplicate and cluster**

   Removes duplicate windows and clusters near-identical candidates.

5. **BLAST against positive genomes**

   Checks whether candidate markers are present and conserved across target genomes.

6. **Select core single-copy candidates**

   Keeps markers that satisfy positive-group conservation and chromosome-only single-copy criteria.

7. **Loose negative screen**

   Applies a relaxed out-group screen to remove borderline cross-reactive markers.

8. **Prepare reference bundle**

   Builds a reference FASTA and GFF bundle for coordinate mapping and annotation.

9. **Map markers back to reference**

   BLASTs final markers against the reference genome and creates BED coordinates.

10. **Annotate markers**

   Intersects marker coordinates with GFF features to identify nearby or overlapping genes.

11. **Summarize marker genes**

   Creates a compact marker summary table.

12. **Design primers**

   Uses Primer3 to design primer pairs for retained markers.

---

## Main outputs

Typical generated outputs include:

```text
results/sliced/
results/filtered/
results/dedup/
results/blast/
results/core/
results/annotation/
results/primer/
results/logs/
results/params/
```

Important final outputs may include:

```text
results/core/<target>_final_markers.fasta
results/core/<target>_marker_annotations.tsv
results/core/<target>_marker_summary.tsv
results/primer/<target>_primers.tsv
results/logs/<run_id>/plakat.log
results/params/<target>.yml
```

These files are generated during analysis and are not tracked by Git.

---

## Software requirements

PLAKAT is intended to run in a Unix-like environment such as Linux, macOS, or Windows WSL2.

Core requirements include:

- Python 3
- Bash
- NCBI Datasets CLI
- Entrez Direct / `efetch`
- BLAST+
- SeqKit
- CD-HIT
- BEDTools
- Primer3
- EMBOSS `seqret` for some reference-preparation steps

Some optional downstream analyses may require additional tools.

---

## Installation notes

Clone the repository:

```bash
git clone https://github.com/mingzapingin/PLAKAT.git
cd PLAKAT
```

Make sure scripts are executable:

```bash
chmod +x scripts/*.sh
```

Check that required tools are available:

```bash
which datasets
which efetch
which blastn
which makeblastdb
which seqkit
which cd-hit-est
which bedtools
which primer3_core
```

---

## Example dry run

Use `--dry` to check input parsing without running the full workflow:

```bash
python plakat.py \
  --positive-list path/to/positive.txt \
  --reference-list path/to/reference.txt \
  --negative-list path/to/negative1.txt path/to/negative2.txt \
  --dry
```

---

## Example full run

```bash
python plakat.py \
  --positive-list path/to/M.marinum.txt \
  --reference-list path/to/M.marinum.txt \
  --negative-list path/to/M.ulcerans.txt path/to/M.chelonae.txt path/to/M.shottsii.txt \
  --params
```

Use `--verbose` to print full tool output to the console:

```bash
python plakat.py \
  -p path/to/positive.txt \
  -r path/to/reference.txt \
  -n path/to/negative1.txt path/to/negative2.txt \
  --params \
  --verbose
```

---

## Git tracking policy

This repository tracks only source code and documentation.

The following are intentionally ignored:

- downloaded genomes;
- BLAST databases;
- cached NCBI files;
- intermediate FASTA files;
- result tables;
- logs;
- temporary input lists;
- large analysis outputs.

This keeps the GitHub repository lightweight and reproducible.

---


## Author

**Kimhun Tuntikawinwong**

PLAKAT was developed for genome-guided discovery of species-specific marker regions in high-similarity bacterial clades, with an initial focus on *Mycobacterium* species.

---

## License

License information has not yet been specified.

Before reuse or redistribution, please contact the author.

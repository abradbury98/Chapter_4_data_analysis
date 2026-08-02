# scripts/python

Python utilities for processing ABRicate AMR outputs, normalising Kraken2 taxonomic classification results, and generating MetaPhlAn visualisations.

---

## Contents

| Script | Purpose |
|--------|---------|
| [`combine_abricate.py`](#combine_abricatepy) | Merge per-sample ABRicate Excel results into combined workbooks |
| [`compare_abricate.py`](#compare_abricatepy) | Compare AMR gene presence/absence across samples or conditions |
| [`hamronize_and_dedup.py`](#hamronize_and_deduppy) | Standardise and deduplicate ABRicate outputs via hAMRonization |
| [`kraken2_bpm_normalised_reads.py`](#kraken2_bpm_normalised_readspy) | BPM normalisation of Kraken2 read-level classification |
| [`kraken2_bpm_normalise_contigs.py`](#kraken2_bpm_normalise_contigspy) | BPM normalisation of Kraken2 contig-level classification |
| [`metaphlan_visualisations.py`](#metaphlan_visualisationspy) | Visualise MetaPhlAn4 community profiles |
| [`normalize_abricate.py`](#normalize_abricatepy) | Normalise ABRicate results by sequencing depth |
| [`normalize_dedup_only.py`](#normalize_dedup_onlypy) | Lightweight deduplication without full hAMRonization |
| [`process_abricate.py`](#process_abricatepy) | Filter and reformat raw ABRicate outputs |

---

## Global Prerequisites

Install core dependencies before running any script:

```bash
pip install pandas numpy matplotlib seaborn openpyxl
```

Scripts that use hAMRonization require an additional install:

```bash
pip install hamronization
```

---

## Scripts

### `combine_abricate.py`

Merges ABRicate results from multiple samples into group-level combined Excel workbooks with formatted output sheets.

#### Inputs

Excel (`.xlsx`) files in `~/Desktop/ABRicate_dereplicated/`, each with a sheet named `ABRicate Results` containing columns: `GENE`, `PRODUCT`, `Type`, `DATABASE`, `%COVERAGE`, `%IDENTITY`, `RESISTANCE`, `ACCESSION`. Paths and sample group definitions are hardcoded — update them before running.

#### Outputs

One `.xlsx` file per sample group (e.g., `COMBINED_ASC_CAT.xlsx`) saved to the same input directory. Each workbook contains:
- Combined results table (one row per unique gene, with `Count` and `Samples` columns)
- Metadata info sheet

#### Dependencies

```bash
pip install pandas openpyxl
```

#### Running

```bash
python scripts/python/combine_abricate.py
```

> **Note:** No command-line arguments — all configuration is set inside the script.

---

### `compare_abricate.py`

Compares AMR gene presence/absence across samples or enrichment conditions using ABRicate outputs.

#### Inputs

ABRicate TSV files (two or more samples/conditions).

#### Outputs

Comparison table or summary of shared and unique gene hits.

#### Dependencies

```bash
pip install pandas
```

#### Running

```bash
python scripts/python/compare_abricate.py --help
```

---

### `hamronize_and_dedup.py`

Standardises ABRicate outputs via hAMRonization and deduplicates overlapping gene hits across databases using coordinate-based single-linkage clustering.

#### Inputs

Folders matching `abricate_reads_*` on the Desktop, each containing barcode subdirectories (e.g., `all_barcode01_porechop_nanofilt/`) with TSV files from multiple ABRicate database scans.

Expected structure:

```
~/Desktop/
└── abricate_reads_<enrichment>/
    ├── all_barcode01_porechop_nanofilt/
    │   ├── results_card.tsv
    │   ├── results_resfinder.tsv
    │   └── ...
    └── all_barcode03_porechop_nanofilt/
        └── ...
```

#### Outputs

Written to `~/Desktop/hamronized/<ENRICHMENT>/`:

| File | Description |
|------|-------------|
| `CH01_enrichment_deduplicated.tsv` | Per-sample deduplicated AMR report |
| `*.json` | hAMRonization intermediate JSON files |
| `combined_enrichment_report.tsv` | All samples merged into a single file |

#### Dependencies

```bash
pip install pandas hamronization
```

#### Running

```bash
python scripts/python/hamronize_and_dedup.py
```

> **Note:** Paths are hardcoded — update the Desktop directory at the top of the script before running.

---

### `kraken2_bpm_normalised_reads.py`

Calculates Bases Per Million (BPM) taxonomic abundance from Kraken2 read-level classification, traversing the NCBI taxonomy tree to aggregate bases at genus, family, and species level.

**Formula:** `BPM = (bases assigned to taxon clade / total bases in sample) × 1,000,000`

#### Inputs

| File | Description |
|------|-------------|
| Kraken2 `.out` files | One per sample; standard Kraken2 per-read output format |
| `nodes.dmp` | NCBI taxonomy node hierarchy |
| `names.dmp` | NCBI taxonomy scientific names |

Download NCBI taxonomy files from: `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`

#### Outputs

Tab-separated file with 12 columns per row: `group`, `sample`, `taxonomic_level`, `taxid`, `taxon_name`, `rank`, `base_count`, `BPM`, `relative_abundance`, `superkingdom`, and summary statistics printed to console.

#### Dependencies

```bash
pip install pandas numpy
```

#### Running

```bash
python scripts/python/kraken2_bpm_normalised_reads.py --help
```

---

### `kraken2_bpm_normalise_contigs.py`

Applies the same BPM normalisation as `kraken2_bpm_normalised_reads.py` but on Kraken2 **contig-level** classification outputs.

#### Inputs

| File | Description |
|------|-------------|
| Kraken2 `.out` files | Contig-level classification (from `kraken2_nt_assembly_contigs.sl`) |
| `nodes.dmp` | NCBI taxonomy node hierarchy |
| `names.dmp` | NCBI taxonomy scientific names |

#### Outputs

BPM-normalised abundance table at contig level (same column format as the read-level script).

#### Dependencies

```bash
pip install pandas numpy
```

#### Running

```bash
python scripts/python/kraken2_bpm_normalise_contigs.py --help
```

---

### `metaphlan_visualisations.py`

Generates taxonomic composition figures from a MetaPhlAn4 profile file.

#### Inputs

| File | Description |
|------|-------------|
| `barcode01_profile.txt` | MetaPhlAn4 output — clade names and relative abundances |

#### Outputs

| File | Description |
|------|-------------|
| `barcode01_phylum_bar.png` | Horizontal bar chart — phylum-level abundance |
| `barcode01_species_bar.png` | Horizontal bar chart — species-level abundance, colour-coded by phylum |
| `barcode01_phylum_pie.png` | Pie chart — phylum composition |
| `barcode01_cladogram.png` | Circular cladogram via GraPhlAn *(optional)* |

#### Dependencies

```bash
pip install matplotlib pandas seaborn
# Optional — for cladogram output:
pip install graphlan export2graphlan
```

#### Running

```bash
python scripts/python/metaphlan_visualisations.py --help
```

---

### `normalize_abricate.py`

Normalises ABRicate AMR results by sequencing depth or sample size to allow cross-sample comparisons.

#### Inputs

ABRicate TSV file(s) and sequencing depth or read count metadata.

#### Outputs

Normalised AMR abundance table (TSV).

#### Dependencies

```bash
pip install pandas numpy
```

#### Running

```bash
python scripts/python/normalize_abricate.py --help
```

---

### `normalize_dedup_only.py`

Lightweight deduplication of ABRicate outputs without the full hAMRonization standardisation step. Use when cross-database standardisation is not required.

#### Inputs

ABRicate TSV file(s).

#### Outputs

Deduplicated TSV file(s).

#### Dependencies

```bash
pip install pandas
```

#### Running

```bash
python scripts/python/normalize_dedup_only.py --help
```

---

### `process_abricate.py`

Filters ABRicate outputs by coverage/identity thresholds, reformats columns, and prepares data for downstream analysis scripts.

#### Inputs

ABRicate TSV file(s).

#### Outputs

Filtered and reformatted TSV file.

#### Dependencies

```bash
pip install pandas
```

#### Running

```bash
python scripts/python/process_abricate.py --help
```

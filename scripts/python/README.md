# scripts/python

Python utilities for AMR result processing, normalization, and taxonomy-based summarization/visualization.

## Files

- `combine_abricate.py`  
  Combines ABRicate output files into consolidated tables.

- `compare_abricate.py`  
  Compares ABRicate results across inputs/conditions.

- `hamronize_and_dedup.py`  
  Harmonizes AMR outputs and removes duplicates.

- `kraken2_bpm_normalise_contigs.py`  
  Normalizes Kraken2 contig-level results (BPM-style normalization).

- `kraken2_bpm_normalised_reads.py`  
  Normalizes Kraken2 read-level results (BPM-style normalization).

- `metaphlan_visualisations.py`  
  Generates visualizations from MetaPhlAn outputs.

- `normalize_abricate.py`  
  Normalizes ABRicate result counts/metrics.

- `normalize_dedup_only.py`  
  Performs deduplication-focused normalization.

- `process_abricate.py`  
  General processing/clean-up pipeline for ABRicate outputs.

## Typical usage

```bash
python scripts/python/<script_name>.py --help
```

If `--help` is not implemented, inspect the script’s argument parser or top-level variables.

## Dependencies

Likely includes `pandas`, `numpy`, and plotting/IO libraries; create an environment and install script-specific requirements.

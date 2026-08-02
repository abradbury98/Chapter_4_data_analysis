"""
Calculate BPM (Bases Per Million) from a Kraken2 .out file.

BPM = (bases assigned to taxon clade / total bases in sample) x 1,000,000

Kraken2 assigns each read to the most specific taxon it can resolve. Reads
assigned to a species (e.g. Comamonas kerstersii) are NOT re-counted at the
genus node unless explicitly rolled up. This script walks each read's taxid up
the NCBI taxonomy tree and aggregates bases at genus level, so genus BPM
reflects ALL reads within that genus clade (not just reads left unresolved at
genus level).

Kraken2 .out columns:
  1. C/U  - classified or unclassified
  2. read_id
  3. taxon name + "(taxid N)"
  4. read_length (actual bases)
  5. kmer hits string

Total bases = sum of ALL read lengths (classified + unclassified).
"""

import re
import sys
import time
from collections import defaultdict

KRAKEN_DIRS = {
    "ASC":      "/Users/alicebradbury/Desktop/kraken.out_report_ASC",
    "ASC_RVS":  "/Users/alicebradbury/Desktop/kraken.out_report_ASC_RVS",
    "ASC_BFBB": "/Users/alicebradbury/Desktop/kraken.out_reports_ASC_BFBB",
    "EVI":      "/Users/alicebradbury/Desktop/kraken.out_report_EVI",
    "EVI_RVS":  "/Users/alicebradbury/Desktop/kraken.out_report_EVI_RVS",
    "EVI_BFBB": "/Users/alicebradbury/Desktop/kraken.out_report_EVI_BFBB",
}
NODES_DMP   = "/Users/alicebradbury/Desktop/AMRplusplus/data/kraken_db/k2_minusb_20250714/nodes.dmp"
NAMES_DMP   = "/Users/alicebradbury/Desktop/AMRplusplus/data/kraken_db/k2_minusb_20250714/names.dmp"
OUTPUT_FILE = "/Users/alicebradbury/Desktop/all_samples_bpm.tsv"

TAXID_RE = re.compile(r"\(taxid (\d+)\)")


def load_taxonomy(nodes_path, names_path):
    """
    Parse nodes.dmp and names.dmp.
    Returns:
        ranks  : taxid (str) -> rank (str)
        parents: taxid (str) -> parent_taxid (str)
        names  : taxid (str) -> scientific name (str)
    """
    ranks, parents = {}, {}
    with open(nodes_path) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 3:
                continue
            taxid  = parts[0].strip()
            parent = parts[1].strip()
            rank   = parts[2].strip()
            ranks[taxid]   = rank
            parents[taxid] = parent

    names = {}
    with open(names_path) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 4:
                continue
            taxid     = parts[0].strip()
            name      = parts[1].strip()
            name_class = parts[3].rstrip("\t|\n").strip()
            if name_class == "scientific name":
                names[taxid] = name

    return ranks, parents, names


def rank_ancestor(taxid, target_rank, ranks, parents):
    """
    Walk up the taxonomy tree from taxid and return the first ancestor
    at target_rank (e.g. 'genus', 'family', 'species').
    Returns None if root is reached without finding one.
    """
    visited = set()
    current = taxid
    while current and current not in visited:
        if ranks.get(current) == target_rank:
            return current
        parent = parents.get(current)
        if parent is None or parent == current:
            break
        visited.add(current)
        current = parent
    return None


# Fixed NCBI taxids for each domain — stable across all taxonomy versions.
# Rank names ("superkingdom", "domain") differ between taxonomy releases,
# but these taxids never change.
DOMAIN_TAXIDS = {
    "2":     "Bacteria",
    "2157":  "Archaea",
    "2759":  "Eukaryota",
    "10239": "Viruses",
}

_superkingdom_cache: dict = {}

def get_superkingdom(taxid, ranks, parents, tax_names):
    """
    Walk up the taxonomy tree until one of the known domain taxids is reached.
    Returns 'Bacteria', 'Archaea', 'Eukaryota', 'Viruses', or 'Unknown'.
    Uses fixed taxid anchors rather than rank names, which vary between
    taxonomy releases (older = 'superkingdom', newer = 'domain').
    Results are cached to avoid repeated tree walks for the same taxid.
    """
    if taxid in _superkingdom_cache:
        return _superkingdom_cache[taxid]
    visited = set()
    current = taxid
    while current and current not in visited:
        if current in DOMAIN_TAXIDS:
            result = DOMAIN_TAXIDS[current]
            _superkingdom_cache[taxid] = result
            return result
        parent = parents.get(current)
        if parent is None or parent == current:
            break
        visited.add(current)
        current = parent
    _superkingdom_cache[taxid] = "Unknown"
    return "Unknown"


def parse_taxon_field(field):
    m = TAXID_RE.search(field)
    taxid = m.group(1) if m else "0"
    return taxid


TAX_LEVELS = ["species", "family", "genus"]


def _read_kraken_file(path, max_retries=3, retry_delay=10):
    """
    Read a Kraken2 .out file into (raw_bases, total_bases).
    Retries up to max_retries times on OSError / TimeoutError before giving up.
    """
    last_exc = None
    for attempt in range(1, max_retries + 1):
        raw_bases   = defaultdict(int)
        total_bases = 0
        try:
            with open(path) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    taxon_field = parts[2]
                    try:
                        read_length = int(parts[3])
                    except ValueError:
                        continue
                    taxid = parse_taxon_field(taxon_field)
                    raw_bases[taxid] += read_length
                    total_bases      += read_length
            return raw_bases, total_bases          # success
        except OSError as exc:
            last_exc = exc
            if attempt < max_retries:
                print(f"  WARNING: attempt {attempt} failed reading {path}: {exc}",
                      file=sys.stderr)
                print(f"  Retrying in {retry_delay}s ...", file=sys.stderr)
                time.sleep(retry_delay)
            else:
                print(f"  ERROR: all {max_retries} attempts failed for {path}: {exc}",
                      file=sys.stderr)
    raise last_exc


def process_sample(kraken_out_path, ranks, parents, tax_names):
    """Parse one .out file and return rows for species, family and genus levels."""
    try:
        raw_bases, total_bases = _read_kraken_file(kraken_out_path)
    except OSError as exc:
        print(f"  SKIPPING {kraken_out_path}: {exc}", file=sys.stderr)
        return []

    if total_bases == 0:
        return []

    rows = []
    for level in TAX_LEVELS:
        level_bases = defaultdict(int)
        for taxid, bases in raw_bases.items():
            if taxid == "0":
                continue
            ancestor = rank_ancestor(taxid, level, ranks, parents)
            if ancestor:
                level_bases[ancestor] += bases

        classified_bases = sum(level_bases.values())
        if classified_bases == 0:
            continue

        for anc_taxid, bases in level_bases.items():
            bpm       = (bases / total_bases) * 1_000_000
            rel_abund = (bases / classified_bases) * 100
            rows.append({
                "tax_level":        level,
                "taxid":            anc_taxid,
                "taxon":            tax_names.get(anc_taxid, f"taxid:{anc_taxid}"),
                "rank":             level,
                "bases":            bases,
                "total_bases":      total_bases,
                "classified_bases": classified_bases,
                "bpm":              round(bpm, 4),
                "rel_abund":        round(rel_abund, 4),
                "superkingdom":     get_superkingdom(anc_taxid, ranks, parents, tax_names),
            })

    return rows


def main():
    import glob
    import os

    print(f"Loading taxonomy: {NODES_DMP}")
    ranks, parents, tax_names = load_taxonomy(NODES_DMP, NAMES_DMP)

    all_rows = []
    total_samples = 0
    for group, kraken_dir in KRAKEN_DIRS.items():
        out_files = sorted(glob.glob(os.path.join(kraken_dir, "*.out")))
        if not out_files:
            print(f"WARNING: no .out files found in {kraken_dir}", file=sys.stderr)
            continue
        print(f"\n--- Group: {group} ---")
        for path in out_files:
            sample = os.path.basename(path).replace("_reads.out", "")
            print(f"Processing: {sample}")
            rows = process_sample(path, ranks, parents, tax_names)
            for r in rows:
                r["sample"] = sample
                r["group"]  = group
            all_rows.extend(rows)
            total_samples += 1
            if rows:
                tb = rows[0]["total_bases"]
                cb = rows[0]["classified_bases"]
                print(f"  total bases: {tb:,}  classified: {cb:,} ({cb/tb*100:.1f}%)  genera: {len(rows)}")

    # Write combined TSV
    header = "group\tsample\ttax_level\ttaxid\ttaxon\trank\tbases\ttotal_bases\tclassified_bases\tbpm\trel_abund\tsuperkingdom"
    with open(OUTPUT_FILE, "w") as fh:
        fh.write(header + "\n")
        for r in sorted(all_rows, key=lambda r: (r["group"], r["sample"], r["tax_level"], -r["rel_abund"])):
            fh.write(
                f"{r['group']}\t{r['sample']}\t{r['tax_level']}\t{r['taxid']}\t{r['taxon']}\t{r['rank']}\t"
                f"{r['bases']}\t{r['total_bases']}\t{r['classified_bases']}\t"
                f"{r['bpm']}\t{r['rel_abund']}\t{r['superkingdom']}\n"
            )

    print(f"\nOutput written to: {OUTPUT_FILE}")
    print(f"Total samples: {total_samples}, total rows: {len(all_rows)}")

    # Sanity check: show superkingdom breakdown
    from collections import Counter
    sk_counts = Counter(r["superkingdom"] for r in all_rows if r["tax_level"] == "genus")
    print("\nSuperkingdom breakdown (genus-level rows):")
    for sk, n in sorted(sk_counts.items(), key=lambda x: -x[1]):
        print(f"  {sk}: {n:,}")


if __name__ == "__main__":
    main()

"""
Calculate BPM (Bases Per Million) and raw base counts from Kraken2 contig .out files.

BPM = (bases in contigs assigned to taxon clade / total contig bases in sample) x 1,000,000

Contigs are assembled sequences; column 4 is the contig length in bases.
The script rolls up bases to genus/family/species level by walking the NCBI
taxonomy tree, so a genus total reflects all contigs classified anywhere within
that genus clade.

Kraken2 .out columns:
  1. C/U  - classified or unclassified
  2. contig_id
  3. taxon name + "(taxid N)"
  4. contig_length (bases)
  5. kmer hits string

Output TSV includes both normalised (BPM, rel_abund) and unnormalised
(raw_bases, contig_count) columns so the R plotting script can use either.
"""

import re
import sys
import time
from collections import defaultdict

KRAKEN_DIRS = {
    "ASC_CAT": "/Users/alicebradbury/Desktop/Kraken2_contigs/ASC_CAT",
    "EVI_CAT": "/Users/alicebradbury/Desktop/Kraken2_contigs/EVI_CAT",
    "ASC_RVS": "/Users/alicebradbury/Desktop/Kraken2_contigs/ASC_RVS",
    "EVI_RVS": "/Users/alicebradbury/Desktop/Kraken2_contigs/EVI_RVS",
}
NODES_DMP   = "/Users/alicebradbury/Desktop/AMRplusplus/data/kraken_db/k2_minusb_20250714/nodes.dmp"
NAMES_DMP   = "/Users/alicebradbury/Desktop/AMRplusplus/data/kraken_db/k2_minusb_20250714/names.dmp"
OUTPUT_FILE = "/Users/alicebradbury/Desktop/all_samples_contig_bpm.tsv"

TAXID_RE = re.compile(r"\(taxid (\d+)\)")


def load_taxonomy(nodes_path, names_path):
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
            taxid      = parts[0].strip()
            name       = parts[1].strip()
            name_class = parts[3].rstrip("\t|\n").strip()
            if name_class == "scientific name":
                names[taxid] = name

    return ranks, parents, names


def rank_ancestor(taxid, target_rank, ranks, parents):
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


DOMAIN_TAXIDS = {
    "2":     "Bacteria",
    "2157":  "Archaea",
    "2759":  "Eukaryota",
    "10239": "Viruses",
}

_superkingdom_cache: dict = {}

def get_superkingdom(taxid, ranks, parents, tax_names):
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


def parse_taxid(field):
    m = TAXID_RE.search(field)
    return m.group(1) if m else "0"


TAX_LEVELS = ["species", "family", "genus"]


def _read_kraken_file(path, max_retries=3, retry_delay=10):
    """Read a contig .out file into (raw_bases, raw_counts, total_bases, total_contigs)."""
    last_exc = None
    for attempt in range(1, max_retries + 1):
        raw_bases   = defaultdict(int)
        raw_counts  = defaultdict(int)
        total_bases   = 0
        total_contigs = 0
        try:
            with open(path) as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 4:
                        continue
                    taxon_field = parts[2]
                    try:
                        contig_len = int(parts[3])
                    except ValueError:
                        continue
                    taxid = parse_taxid(taxon_field)
                    raw_bases[taxid]  += contig_len
                    raw_counts[taxid] += 1
                    total_bases   += contig_len
                    total_contigs += 1
            return raw_bases, raw_counts, total_bases, total_contigs
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
    """Parse one contig .out file; return rows at species, family and genus level."""
    try:
        raw_bases, raw_counts, total_bases, total_contigs = _read_kraken_file(kraken_out_path)
    except OSError as exc:
        print(f"  SKIPPING {kraken_out_path}: {exc}", file=sys.stderr)
        return []

    if total_bases == 0:
        return []

    rows = []
    for level in TAX_LEVELS:
        level_bases  = defaultdict(int)
        level_counts = defaultdict(int)

        for taxid, bases in raw_bases.items():
            if taxid == "0":
                continue
            ancestor = rank_ancestor(taxid, level, ranks, parents)
            if ancestor:
                level_bases[ancestor]  += bases
                level_counts[ancestor] += raw_counts[taxid]

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
                "bases":            bases,
                "contig_count":     level_counts[anc_taxid],
                "total_bases":      total_bases,
                "total_contigs":    total_contigs,
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
            # Strip _m1000_kraken2 suffix to get clean sample name
            sample = re.sub(r"_m\d+_kraken2$", "", os.path.basename(path).replace(".out", ""))
            print(f"Processing: {sample}")
            rows = process_sample(path, ranks, parents, tax_names)
            for r in rows:
                r["sample"] = sample
                r["group"]  = group
            all_rows.extend(rows)
            total_samples += 1
            if rows:
                tb = rows[0]["total_bases"]
                tc = rows[0]["total_contigs"]
                cb = rows[0]["classified_bases"]
                print(f"  total contigs: {tc:,}  total bases: {tb:,}  "
                      f"classified: {cb:,} ({cb/tb*100:.1f}%)")

    header = (
        "group\tsample\ttax_level\ttaxid\ttaxon\t"
        "bases\tcontig_count\ttotal_bases\ttotal_contigs\tclassified_bases\t"
        "bpm\trel_abund\tsuperkingdom"
    )
    with open(OUTPUT_FILE, "w") as fh:
        fh.write(header + "\n")
        for r in sorted(all_rows, key=lambda r: (r["group"], r["sample"], r["tax_level"], -r["rel_abund"])):
            fh.write(
                f"{r['group']}\t{r['sample']}\t{r['tax_level']}\t{r['taxid']}\t{r['taxon']}\t"
                f"{r['bases']}\t{r['contig_count']}\t{r['total_bases']}\t{r['total_contigs']}\t"
                f"{r['classified_bases']}\t{r['bpm']}\t{r['rel_abund']}\t{r['superkingdom']}\n"
            )

    print(f"\nOutput written to: {OUTPUT_FILE}")
    print(f"Total samples: {total_samples}, total rows: {len(all_rows)}")

    from collections import Counter
    sk_counts = Counter(r["superkingdom"] for r in all_rows if r["tax_level"] == "genus")
    print("\nSuperkingdom breakdown (genus-level rows):")
    for sk, n in sorted(sk_counts.items(), key=lambda x: -x[1]):
        print(f"  {sk}: {n:,}")


if __name__ == "__main__":
    main()

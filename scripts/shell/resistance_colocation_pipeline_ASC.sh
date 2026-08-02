#!/usr/bin/env bash
# =============================================================================
# Resistance Gene Co-location Analysis Pipeline
# =============================================================================
# Analyses all contigs in a folder for co-located resistance genes using:
#   - AMRFinderPlus
#   - RGI (CARD)
#   - Prokka
#   - IntegronFinder
#   - ISEScan
#   - PlasmidFinder
#   - MOB-suite
#   - Summary report generation
#
# Usage:
#   chmod +x resistance_colocation_pipeline.sh
#   ./resistance_colocation_pipeline.sh /path/to/contig/folder
#
# Default contig folder: ~/Desktop/contigs
# =============================================================================

set -euo pipefail

# ── Colour helpers ────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BLUE='\033[0;34m'; CYAN='\033[0;36m'; NC='\033[0m'

log()  { echo -e "${GREEN}[$(date '+%H:%M:%S')] ✔ $*${NC}"; }
info() { echo -e "${CYAN}[$(date '+%H:%M:%S')] ➜ $*${NC}"; }
warn() { echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠ $*${NC}"; }
err()  { echo -e "${RED}[$(date '+%H:%M:%S')] ✘ $*${NC}" >&2; }
section() { echo -e "\n${BLUE}════════════════════════════════════════${NC}"; \
            echo -e "${BLUE}  $*${NC}"; \
            echo -e "${BLUE}════════════════════════════════════════${NC}\n"; }

# Source conda.sh then activate the named env.
# Uses `conda info --base` to find the real install location (handles /opt/anaconda3,
# /opt/miniconda3, $HOME/anaconda3, etc. without hardcoding paths).
_conda_activate() {
    local env="$1"
    local base
    base="$(conda info --base 2>/dev/null)" || {
        # conda not on PATH yet — try common locations as a last resort
        local p
        for p in \
            "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" \
            "/opt/miniconda3" "/opt/anaconda3" "/opt/conda"; do
            [[ -f "$p/etc/profile.d/conda.sh" ]] && { base="$p"; break; }
        done
    }
    [[ -z "$base" ]] && { echo "ERROR: cannot locate conda base" >&2; return 1; }
    # shellcheck source=/dev/null
    source "$base/etc/profile.d/conda.sh"
    conda activate "$env"
}

# ── Input / output paths ──────────────────────────────────────────────────────
CONTIG_DIR="${1:-/Users/alicebradbury/Desktop/ASC_BFBB_Contigs_resistance_co-localisation}"
OUTPUT_DIR="${2:-/Users/alicebradbury/Desktop/ASC_BFBB_resistance_analysis}"
THREADS="${3:-$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)}"
CONDA_ENV="resistance_analysis"

# ── Validate input ────────────────────────────────────────────────────────────
if [[ ! -d "$CONTIG_DIR" ]]; then
    err "Contig directory not found: $CONTIG_DIR"
    echo "Usage: $0 [contig_dir] [output_dir] [threads]"
    echo "Example: $0 ~/Desktop/contigs ~/Desktop/results 8"
    exit 1
fi

CONTIG_COUNT=$(find "$CONTIG_DIR" -maxdepth 1 -name "*.fasta" 2>/dev/null | wc -l | tr -d ' ')
if [[ "$CONTIG_COUNT" -eq 0 ]]; then
    err "No .fasta files found in: $CONTIG_DIR"
    exit 1
fi

# ── Detect OS ─────────────────────────────────────────────────────────────────
OS="$(uname -s)"
case "$OS" in
    Linux*)  PLATFORM="linux" ;;
    Darwin*) PLATFORM="mac" ;;
    *)       err "Unsupported OS: $OS"; exit 1 ;;
esac

# =============================================================================
# STEP 0 — Install Conda / Mamba if absent
# =============================================================================
install_conda() {
    section "STEP 0 — Installing Conda"

    if command -v conda &>/dev/null; then
        log "Conda already installed: $(conda --version)"
        return
    fi

    info "Conda not found — installing Miniconda..."
    local installer
    if [[ "$PLATFORM" == "mac" ]]; then
        # Detect Apple Silicon vs Intel
        if [[ "$(uname -m)" == "arm64" ]]; then
            installer="Miniconda3-latest-MacOSX-arm64.sh"
        else
            installer="Miniconda3-latest-MacOSX-x86_64.sh"
        fi
    else
        installer="Miniconda3-latest-Linux-x86_64.sh"
    fi

    curl -fsSL "https://repo.anaconda.com/miniconda/$installer" -o /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
    rm /tmp/miniconda.sh

    # Initialise for current shell session
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda init bash 2>/dev/null || true
    [[ -f "$HOME/.zshrc" ]] && (conda init zsh 2>/dev/null || true)

    # Install mamba for faster dependency resolution
    conda install -y -n base -c conda-forge mamba
    log "Conda + Mamba installed successfully"
}

# =============================================================================
# STEP 1 — Create conda environment with all tools
# =============================================================================
create_environment() {
    section "STEP 1 — Creating Conda Environment: $CONDA_ENV"

    # Source conda if not yet active
    if ! command -v conda &>/dev/null; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    fi

    if conda env list | grep -q "^${CONDA_ENV} "; then
        warn "Environment '$CONDA_ENV' already exists — skipping creation"
        warn "To rebuild: conda env remove -n $CONDA_ENV"
    else
        info "Creating environment (this may take 10–20 minutes)..."

        # Use mamba if available, fall back to conda
        PKG_MGR="conda"
        command -v mamba &>/dev/null && PKG_MGR="mamba"

        # Core packages — all conda-solvable on osx-arm64.
        # rgi, prokka, and isescan are handled separately below.
        $PKG_MGR create -y -n "$CONDA_ENV" \
            -c conda-forge -c bioconda -c defaults \
            python=3.10 \
            ncbi-amrfinderplus \
            integron_finder \
            plasmidfinder \
            mob_suite \
            samtools \
            biopython \
            pandas \
            matplotlib \
            seaborn \
            blast \
            hmmer \
            prodigal

        log "Base environment created"

        # ── ISEScan: requires TensorFlow 1.x which has no ARM64 support.
        #   Skipped on macOS ARM64; on Linux/Intel it can be pip-installed.
        if [[ "$PLATFORM" == "linux" ]] || [[ "$(uname -m)" != "arm64" ]]; then
            info "Installing ISEScan via pip..."
            conda run -n "$CONDA_ENV" pip install isescan \
            && log "ISEScan installed" \
            || warn "ISEScan pip install failed — ISEScan step will be skipped"
        else
            warn "ISEScan skipped — TensorFlow 1.x dependency has no ARM64 support"
        fi

        # ── RGI: not on PyPI; install from CARD GitHub repo.
        #   NOTE: do NOT use CONDA_SUBDIR=osx-64 — mixing x86_64 packages
        #   into an arm64 env corrupts shared libraries (libexpat etc).
        #   Pin numpy/pandas after install — RGI upgrades them beyond what
        #   mob-suite supports (numpy<1.23.5, pandas<=1.5.3).
        info "Installing RGI from GitHub..."
        conda run -n "$CONDA_ENV" pip install \
            "git+https://github.com/arpcard/rgi.git" \
        && conda run -n "$CONDA_ENV" pip install \
            "numpy<1.23.5" "pandas==1.5.3" "dask[dataframe]<2024.1.0" \
        && log "RGI installed" \
        || warn "RGI install failed — CARD/RGI steps will be skipped"

        # ── Prokka / Bakta annotation ─────────────────────────────────────────
        # Prokka's osx-arm64 conda package has missing deps. Installing it with
        # CONDA_SUBDIR=osx-64 contaminates shared libs (breaks libexpat/Python).
        # Use bakta instead — it is ARM64-native and actively maintained.
        info "Installing Bakta (ARM64-native annotation tool)..."
        $PKG_MGR install -y -n "$CONDA_ENV" -c conda-forge -c bioconda bakta \
        && log "Bakta installed" \
        || {
            warn "Bakta install failed — trying prokka (Linux/Intel only)..."
            $PKG_MGR install -y -n "$CONDA_ENV" -c conda-forge -c bioconda prokka \
            && log "Prokka installed" \
            || warn "Annotation tool install failed — annotation step will be skipped"
        }

        # ── IntegronFinder: conda 2.0rc6 uses Bio.Seq.IUPAC removed in
        #   BioPython 1.78. Upgrade to the PyPI release which is compatible.
        info "Upgrading IntegronFinder to BioPython-compatible version..."
        conda run -n "$CONDA_ENV" pip install --upgrade "integron_finder" \
        && log "IntegronFinder upgraded" \
        || warn "IntegronFinder upgrade failed — integron detection may not work"

        log "Environment setup complete"
    fi
}

# =============================================================================
# STEP 2 — Download / update databases
# =============================================================================
setup_databases() {
    section "STEP 2 — Setting Up Databases"

    _conda_activate "$CONDA_ENV"

    local DB_DIR="$OUTPUT_DIR/databases"
    mkdir -p "$DB_DIR"

    # AMRFinderPlus database
    info "Updating AMRFinderPlus database..."
    amrfinder --update 2>/dev/null || {
        warn "AMRFinderPlus auto-update failed; trying manual download..."
        amrfinder_update --database "$DB_DIR/amrfinder_db" || \
            warn "AMRFinderPlus DB update failed — check internet connection"
    }

    # CARD database for RGI
    local CARD_DIR="$DB_DIR/card"
    if command -v rgi &>/dev/null; then
        info "Downloading CARD database for RGI..."
        mkdir -p "$CARD_DIR"
        if [[ ! -f "$CARD_DIR/card.json" ]]; then
            curl -fsSL "https://card.mcmaster.ca/latest/data" -o "$CARD_DIR/card-data.tar.bz2"
            tar -xjf "$CARD_DIR/card-data.tar.bz2" -C "$CARD_DIR"
            rgi load --card_json "$CARD_DIR/card.json" --local \
            && log "CARD database loaded" \
            || warn "rgi load failed — RGI analysis may be incomplete"
        else
            warn "CARD database already present — skipping download"
        fi
    else
        warn "rgi not installed — skipping CARD database download"
    fi

    # PlasmidFinder database
    info "Downloading PlasmidFinder database..."
    local PF_DIR="$DB_DIR/plasmidfinder_db"
    if [[ ! -d "$PF_DIR" ]]; then
        download-db.sh "$PF_DIR" 2>/dev/null || \
        python -m plasmidfinder.database --update "$PF_DIR" 2>/dev/null || \
            warn "PlasmidFinder DB download failed — run manually if needed"
    else
        warn "PlasmidFinder DB already present — skipping"
    fi

    # MOB-suite database
    info "Setting up MOB-suite database..."
    mob_init 2>/dev/null || warn "MOB-suite init failed — may already be set up"

    # Bakta database (light ~2 GB; skip if already present)
    local BAKTA_DB="$DB_DIR/bakta_db"
    if command -v bakta_db &>/dev/null; then
        if [[ ! -d "$BAKTA_DB/db-light" ]]; then
            info "Downloading Bakta database (light, ~2 GB)..."
            bakta_db download --type light --output "$BAKTA_DB" \
            && log "Bakta database downloaded" \
            || warn "Bakta DB download failed — annotation step will be skipped"
        else
            warn "Bakta database already present — skipping"
        fi
    else
        warn "bakta_db not found — skipping Bakta database download"
    fi

    log "Database setup complete"
}

# =============================================================================
# STEP 3 — Run analysis on each contig
# =============================================================================
run_analysis() {
    section "STEP 3 — Running Analysis on $CONTIG_COUNT Contigs"

    _conda_activate "$CONDA_ENV"

    local CONTIG_NUM=0
    local FAILED=()

    for CONTIG in "$CONTIG_DIR"/*.fasta; do
        [[ -f "$CONTIG" ]] || continue
        CONTIG_NUM=$((CONTIG_NUM + 1))
        SAMPLE=$(basename "$CONTIG" .fasta)

        echo ""
        info "[$CONTIG_NUM/$CONTIG_COUNT] Processing: $SAMPLE"

        local SAMPLE_DIR="$OUTPUT_DIR/results/$SAMPLE"
        mkdir -p "$SAMPLE_DIR"

        # ── 3a. AMRFinderPlus ─────────────────────────────────────────────────
        info "  Running AMRFinderPlus..."
        amrfinder \
            --nucleotide "$CONTIG" \
            --output "$SAMPLE_DIR/amrfinder.tsv" \
            --threads "$THREADS" \
            --name "$SAMPLE" \
            2>"$SAMPLE_DIR/amrfinder.log" \
        && log "  AMRFinderPlus: done" \
        || { warn "  AMRFinderPlus: failed (see $SAMPLE_DIR/amrfinder.log)"; FAILED+=("$SAMPLE:amrfinder"); }

        # ── 3b. RGI (CARD) ────────────────────────────────────────────────────
        info "  Running RGI (CARD)..."
        rgi main \
            --input_sequence "$CONTIG" \
            --output_file "$SAMPLE_DIR/rgi_output" \
            --input_type contig \
            --threads "$THREADS" \
            --clean \
            --local \
            2>"$SAMPLE_DIR/rgi.log" \
        && log "  RGI: done" \
        || { warn "  RGI: failed (see $SAMPLE_DIR/rgi.log)"; FAILED+=("$SAMPLE:rgi"); }

        # ── 3c. Annotation (Prokka preferred; bakta fallback on ARM64) ───────
        info "  Running annotation..."
        if command -v prokka &>/dev/null; then
            prokka \
                --outdir "$SAMPLE_DIR/prokka" \
                --prefix "$SAMPLE" \
                --cpus "$THREADS" \
                --force \
                --quiet \
                "$CONTIG" \
                2>"$SAMPLE_DIR/prokka.log" \
            && log "  Prokka: done" \
            || { warn "  Prokka: failed (see $SAMPLE_DIR/prokka.log)"; FAILED+=("$SAMPLE:prokka"); }
        elif command -v bakta &>/dev/null; then
            bakta \
                --output "$SAMPLE_DIR/prokka" \
                --prefix "$SAMPLE" \
                --threads "$THREADS" \
                --force \
                --db "$OUTPUT_DIR/databases/bakta_db/db-light" \
                "$CONTIG" \
                2>"$SAMPLE_DIR/prokka.log" \
            && log "  Bakta (prokka fallback): done" \
            || { warn "  Bakta: failed (see $SAMPLE_DIR/prokka.log)"; FAILED+=("$SAMPLE:annotation"); }
        else
            warn "  Neither prokka nor bakta available — skipping annotation"
            FAILED+=("$SAMPLE:annotation")
        fi

        # ── 3d. IntegronFinder ────────────────────────────────────────────────
        info "  Running IntegronFinder..."
        integron_finder \
            "$CONTIG" \
            --outdir "$SAMPLE_DIR/integrons" \
            --cpu "$THREADS" \
            --linear \
            2>"$SAMPLE_DIR/integron_finder.log" \
        && log "  IntegronFinder: done" \
        || { warn "  IntegronFinder: failed (see $SAMPLE_DIR/integron_finder.log)"; FAILED+=("$SAMPLE:integron_finder"); }

        # ── 3e. ISEScan (insertion sequences) ────────────────────────────────
        if command -v isescan.py &>/dev/null; then
            info "  Running ISEScan..."
            isescan.py \
                --seqfile "$CONTIG" \
                --output "$SAMPLE_DIR/isescan" \
                --nthread "$THREADS" \
                2>"$SAMPLE_DIR/isescan.log" \
            && log "  ISEScan: done" \
            || { warn "  ISEScan: failed (see $SAMPLE_DIR/isescan.log)"; FAILED+=("$SAMPLE:isescan"); }
        else
            warn "  ISEScan: not available (TensorFlow 1.x required; no ARM64 support) — skipping"
        fi

        # ── 3f. PlasmidFinder ────────────────────────────────────────────────
        info "  Running PlasmidFinder..."
        # Detect DB: prefer custom path, fall back to conda's default share location
        PF_DB="$OUTPUT_DIR/databases/plasmidfinder_db"
        if [[ ! -d "$PF_DB" ]]; then
            PF_DB=$(find "$(conda info --base)/envs/$CONDA_ENV/share" \
                -maxdepth 3 -type d -name "database" -path "*plasmidfinder*" \
                2>/dev/null | head -1)
        fi
        if [[ -n "$PF_DB" ]]; then
            mkdir -p "$SAMPLE_DIR/plasmidfinder"
            plasmidfinder.py \
                -i "$CONTIG" \
                -o "$SAMPLE_DIR/plasmidfinder" \
                -p "$PF_DB" \
                2>"$SAMPLE_DIR/plasmidfinder.log" \
            && log "  PlasmidFinder: done" \
            || { warn "  PlasmidFinder: failed (see $SAMPLE_DIR/plasmidfinder.log)"; FAILED+=("$SAMPLE:plasmidfinder"); }
        else
            warn "  PlasmidFinder: database not found — skipping"
            FAILED+=("$SAMPLE:plasmidfinder")
        fi

        # ── 3g. MOB-suite (plasmid typing) ───────────────────────────────────
        # mob-suite needs pandas<=1.5.3; integron_finder needs pandas>=2.
        # If current pandas is too new, use the isolated mob_suite_env instead.
        info "  Running MOB-suite..."
        local MOB_CMD="mob_recon"
        local pandas_compat
        pandas_compat=$(python -c \
            "import pandas as pd; v=pd.__version__.split('.')[:2]; \
             print('yes' if int(v[0]) < 2 else 'no')" 2>/dev/null)
        if [[ "$pandas_compat" != "yes" ]]; then
            MOB_CMD="conda run -n mob_suite_env mob_recon"
        fi
        $MOB_CMD \
            --infile "$CONTIG" \
            --outdir "$SAMPLE_DIR/mob_suite" \
            --num_threads "$THREADS" \
            2>"$SAMPLE_DIR/mob_suite.log" \
        && log "  MOB-suite: done" \
        || { warn "  MOB-suite: failed (see $SAMPLE_DIR/mob_suite.log)"; FAILED+=("$SAMPLE:mob_suite"); }

        log "  [$CONTIG_NUM/$CONTIG_COUNT] $SAMPLE — complete"
    done

    if [[ ${#FAILED[@]} -gt 0 ]]; then
        warn "The following tool runs failed:"
        for entry in "${FAILED[@]}"; do
            warn "  $entry"
        done
    else
        log "All tools completed successfully"
    fi
}

# =============================================================================
# STEP 4 — Generate co-location summary report
# =============================================================================
generate_report() {
    section "STEP 4 — Generating Co-location Summary Report"

    local REPORT_SCRIPT="$OUTPUT_DIR/generate_report.py"
    local RESULTS_DIR="$OUTPUT_DIR/results"

    cat > "$REPORT_SCRIPT" << 'PYEOF'
#!/usr/bin/env python3
"""
Resistance Gene Co-location Summary Report
Parses AMRFinderPlus and RGI outputs to identify co-located resistance genes.
"""

import os, sys, json, glob, csv
from pathlib import Path
from collections import defaultdict
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

RESULTS_DIR = Path(sys.argv[1])
OUTPUT_DIR  = RESULTS_DIR.parent
REPORT_DIR  = OUTPUT_DIR / "report"
REPORT_DIR.mkdir(exist_ok=True)

# ── Colour palette ────────────────────────────────────────────────────────────
COLOURS = {
    'AMRFinderPlus': '#2196F3',
    'RGI':           '#FF5722',
    'both':          '#4CAF50',
    'integron':      '#9C27B0',
    'IS_element':    '#FF9800',
    'plasmid':       '#795548',
}

def parse_amrfinder(tsv_path):
    """Parse AMRFinderPlus output TSV."""
    genes = []
    try:
        df = pd.read_csv(tsv_path, sep='\t')
        for _, row in df.iterrows():
            genes.append({
                'gene':    row.get('Gene symbol', row.get('Element symbol', 'unknown')),
                'contig':  row.get('Contig id', row.get('Sequence name', 'unknown')),
                'start':   int(row.get('Start', 0)),
                'end':     int(row.get('Stop', 0)),
                'strand':  row.get('Strand', '.'),
                'class':   row.get('Class', 'unknown'),
                'subclass':row.get('Subclass', 'unknown'),
                'source':  'AMRFinderPlus',
                'identity':float(row.get('% Identity to reference sequence', 0)),
            })
    except Exception as e:
        print(f"  Warning: could not parse {tsv_path}: {e}")
    return genes

def parse_rgi(txt_path):
    """Parse RGI output TSV."""
    genes = []
    txt_path = Path(txt_path)
    if not txt_path.exists():
        return genes
    try:
        df = pd.read_csv(txt_path, sep='\t')
        for _, row in df.iterrows():
            orf = str(row.get('ORF_ID', ''))
            parts = orf.rsplit('_', 2)
            contig = parts[0] if len(parts) >= 3 else orf
            start  = int(parts[1]) if len(parts) >= 3 else 0
            end    = int(parts[2]) if len(parts) >= 3 else 0
            genes.append({
                'gene':    row.get('Best_Hit_ARO', 'unknown'),
                'contig':  contig,
                'start':   start,
                'end':     end,
                'strand':  row.get('Orientation', '.'),
                'class':   row.get('AMR Gene Family', 'unknown'),
                'subclass':row.get('Drug Class', 'unknown'),
                'source':  'RGI',
                'identity':float(row.get('Best_Identities', 0)),
            })
    except Exception as e:
        print(f"  Warning: could not parse RGI output: {e}")
    return genes

def parse_integrons(sample_dir):
    """Parse IntegronFinder summary."""
    integrons = []
    for f in list(Path(sample_dir).glob("integrons/**/*.integrons")) + \
             list(Path(sample_dir).glob("integrons/**/*.summary")):
        try:
            df = pd.read_csv(f, sep='\t', comment='#')
            for _, row in df.iterrows():
                integrons.append({
                    'contig': row.get('ID_replicon', 'unknown'),
                    'start':  int(row.get('pos_beg', 0)),
                    'end':    int(row.get('pos_end', 0)),
                    'type':   row.get('type_integron', 'unknown'),
                    'element':row.get('annotation', ''),
                })
        except Exception:
            pass
    return integrons

def parse_is_elements(sample_dir):
    """Parse ISEScan predictions."""
    is_elements = []
    for f in Path(sample_dir).glob("isescan/*.tsv"):
        try:
            df = pd.read_csv(f, sep='\t')
            for _, row in df.iterrows():
                is_elements.append({
                    'contig': str(row.get('seqid', 'unknown')),
                    'start':  int(row.get('isBegin', 0)),
                    'end':    int(row.get('isEnd', 0)),
                    'family': row.get('family', 'unknown'),
                    'type':   row.get('type', 'unknown'),
                })
        except Exception:
            pass
    return is_elements

def find_colocation(genes, proximity_bp=50000):
    """
    Find pairs of resistance genes on the same contig within proximity_bp.
    Returns a list of co-located pairs with distance.
    """
    pairs = []
    by_contig = defaultdict(list)
    for g in genes:
        by_contig[g['contig']].append(g)

    for contig, contig_genes in by_contig.items():
        if len(contig_genes) < 2:
            continue
        contig_genes.sort(key=lambda x: x['start'])
        for i in range(len(contig_genes)):
            for j in range(i + 1, len(contig_genes)):
                g1, g2 = contig_genes[i], contig_genes[j]
                if g1['gene'] == g2['gene']:
                    continue
                dist = g2['start'] - g1['end']
                if dist <= proximity_bp:
                    pairs.append({
                        'contig':   contig,
                        'gene_1':   g1['gene'],
                        'class_1':  g1['class'],
                        'start_1':  g1['start'],
                        'end_1':    g1['end'],
                        'gene_2':   g2['gene'],
                        'class_2':  g2['class'],
                        'start_2':  g2['start'],
                        'end_2':    g2['end'],
                        'distance_bp': max(0, dist),
                        'source':   f"{g1['source']} + {g2['source']}",
                    })
    return pairs

def _draw_gene_track(ax, genes, integrons, is_elements, contig, length, show_labels=True):
    """Draw genes as filled rectangles on ax. Returns legend handles."""
    legend_handles = []
    seen_sources = set()

    # Backbone
    ax.barh(0, length, height=0.12, color='#4a4a6a', zorder=2, left=0)

    # IS elements (wide background band)
    for is_el in is_elements:
        if is_el['contig'] != contig:
            continue
        ax.barh(0, is_el['end'] - is_el['start'], left=is_el['start'],
                height=0.55, color=COLOURS['IS_element'], alpha=0.35, zorder=3)
        if 'IS_element' not in seen_sources:
            legend_handles.append(mpatches.Patch(
                color=COLOURS['IS_element'], alpha=0.6, label='IS Element'))
            seen_sources.add('IS_element')

    # Integrons (wide background band, different y)
    for integ in integrons:
        if integ['contig'] != contig:
            continue
        w = max(integ['end'] - integ['start'], 1)
        ax.barh(0, w, left=integ['start'],
                height=0.75, color=COLOURS['integron'], alpha=0.2, zorder=3)
        if 'integron' not in seen_sources:
            legend_handles.append(mpatches.Patch(
                color=COLOURS['integron'], alpha=0.5, label='Integron'))
            seen_sources.add('integron')

    # Resistance genes — filled rectangles with gene-height arrows
    row_y = {'AMRFinderPlus': 0.32, 'RGI': -0.32}
    min_vis_w = length * 0.008   # minimum 0.8 % of contig for visibility
    for g in genes:
        if g['contig'] != contig:
            continue
        y   = row_y.get(g['source'], 0)
        col = COLOURS.get(g['source'], '#aaaaaa')
        w   = max(g['end'] - g['start'], min_vis_w)
        # Filled rectangle
        ax.barh(y, w, left=g['start'], height=0.22,
                color=col, alpha=0.9, zorder=5)
        # Direction tick
        if g['strand'] in ('+', '1', 1, '1'):
            ax.annotate('', xy=(g['start'] + w, y),
                        xytext=(g['start'] + w * 0.7, y),
                        arrowprops=dict(arrowstyle='->', color='white',
                                        lw=1.2, mutation_scale=10),
                        zorder=6)
        else:
            ax.annotate('', xy=(g['start'], y),
                        xytext=(g['start'] + w * 0.3, y),
                        arrowprops=dict(arrowstyle='->', color='white',
                                        lw=1.2, mutation_scale=10),
                        zorder=6)
        if show_labels:
            ax.text(g['start'] + w / 2, y + 0.15,
                    g['gene'], ha='center', va='bottom',
                    fontsize=8, color='white', fontweight='bold',
                    rotation=30, zorder=7)
        if g['source'] not in seen_sources:
            legend_handles.append(mpatches.Patch(
                color=col, label=g['source']))
            seen_sources.add(g['source'])

    return legend_handles


def draw_contig_map(sample, contig, genes, integrons, is_elements, outpath, contig_len=None):
    """Draw a two-panel map: full contig overview + zoomed gene cluster view."""
    contig_genes  = [g for g in genes      if g['contig'] == contig]
    contig_integ  = [i for i in integrons  if i['contig'] == contig]
    contig_is     = [s for s in is_elements if s['contig'] == contig]

    all_positions = ([g['end'] for g in contig_genes] +
                     [i['end'] for i in contig_integ] +
                     [s['end'] for s in contig_is])
    length = contig_len or (max(all_positions) + 5000 if all_positions else 100000)

    # ── decide zoom window ────────────────────────────────────────────────────
    gene_starts = [g['start'] for g in contig_genes]
    gene_ends   = [g['end']   for g in contig_genes]
    if gene_starts:
        z_centre = (min(gene_starts) + max(gene_ends)) / 2
        z_span   = max(max(gene_ends) - min(gene_starts), 5000) * 3
        z_lo, z_hi = max(0, z_centre - z_span / 2), min(length, z_centre + z_span / 2)
    else:
        z_lo, z_hi = 0, length

    has_zoom = (z_hi - z_lo) < length * 0.5

    nrows = 2 if has_zoom else 1
    fig, axes = plt.subplots(nrows, 1,
                             figsize=(14, 2.8 * nrows),
                             gridspec_kw={'hspace': 0.55})
    if nrows == 1:
        axes = [axes]
    fig.patch.set_facecolor('#1a1a2e')

    # ── Panel 1: full contig ──────────────────────────────────────────────────
    ax0 = axes[0]
    ax0.set_facecolor('#1a1a2e')
    handles = _draw_gene_track(ax0, contig_genes, contig_integ, contig_is,
                               contig, length, show_labels=False)

    # Mark zoom region on overview
    if has_zoom:
        ax0.axvspan(z_lo, z_hi, ymin=0.05, ymax=0.95,
                    color='white', alpha=0.07, zorder=1, label='Zoom region')
        ax0.annotate('', xy=(z_lo, 0.62), xytext=(z_hi, 0.62),
                     xycoords=('data', 'axes fraction'),
                     textcoords=('data', 'axes fraction'),
                     arrowprops=dict(arrowstyle='<->', color='#aaaaaa', lw=0.8))

    ax0.set_xlim(-length * 0.01, length * 1.01)
    ax0.set_ylim(-0.65, 0.75)
    ax0.set_xlabel('Position (bp)', color='white', fontsize=8)
    ax0.set_title(f'{sample} — {contig}  (full, {length:,} bp)',
                  color='white', fontsize=9)
    ax0.tick_params(colors='white', labelsize=7)
    for sp in ax0.spines.values():
        sp.set_edgecolor('#555577')
    if handles:
        ax0.legend(handles=handles, loc='upper left',
                   facecolor='#2a2a4e', labelcolor='white', fontsize=7,
                   framealpha=0.8)

    # ── Panel 2: zoomed cluster ───────────────────────────────────────────────
    if has_zoom:
        ax1 = axes[1]
        ax1.set_facecolor('#1a1a2e')
        _draw_gene_track(ax1, contig_genes, contig_integ, contig_is,
                         contig, length, show_labels=True)
        ax1.set_xlim(z_lo, z_hi)
        ax1.set_ylim(-0.65, 0.75)
        ax1.set_xlabel('Position (bp)', color='white', fontsize=8)
        ax1.set_title(f'Zoomed: {z_lo:,.0f} – {z_hi:,.0f} bp',
                      color='white', fontsize=9)
        ax1.tick_params(colors='white', labelsize=7)
        for sp in ax1.spines.values():
            sp.set_edgecolor('#555577')

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()

# ── Main processing ───────────────────────────────────────────────────────────
print("\nGenerating co-location report...")

all_samples = []
all_pairs   = []

for sample_dir in sorted(RESULTS_DIR.iterdir()):
    if not sample_dir.is_dir():
        continue
    sample = sample_dir.name
    print(f"  Processing sample: {sample}")

    # Parse tool outputs
    genes     = []
    amr_file  = sample_dir / "amrfinder.tsv"
    if amr_file.exists():
        genes += parse_amrfinder(amr_file)

    genes += parse_rgi(sample_dir / "rgi_output.txt")

    integrons   = parse_integrons(sample_dir)
    is_elements = parse_is_elements(sample_dir)

    # Find co-located pairs
    pairs = find_colocation(genes)

    # Draw contig maps
    map_dir = REPORT_DIR / "maps" / sample
    map_dir.mkdir(parents=True, exist_ok=True)
    contigs_with_genes = set(g['contig'] for g in genes)
    for contig in contigs_with_genes:
        c_genes  = [g for g in genes      if g['contig'] == contig]
        c_integ  = [i for i in integrons  if i['contig'] == contig]
        c_is     = [s for s in is_elements if s['contig'] == contig]
        outpath  = map_dir / f"{contig.replace('/', '_')}.png"
        draw_contig_map(sample, contig, c_genes, c_integ, c_is, outpath)

    sample_data = {
        'sample':         sample,
        'n_genes':        len(genes),
        'n_pairs':        len(pairs),
        'n_integrons':    len(integrons),
        'n_is_elements':  len(is_elements),
        'genes':          genes,
        'pairs':          pairs,
    }
    all_samples.append(sample_data)
    for p in pairs:
        p['sample'] = sample
        all_pairs.append(p)

# ── TSV summary tables ────────────────────────────────────────────────────────
genes_rows = []
for s in all_samples:
    for g in s['genes']:
        g['sample'] = s['sample']
        genes_rows.append(g)

if genes_rows:
    pd.DataFrame(genes_rows).to_csv(REPORT_DIR / "all_resistance_genes.tsv",
                                    sep='\t', index=False)
    print(f"  Wrote all_resistance_genes.tsv ({len(genes_rows)} entries)")

if all_pairs:
    pd.DataFrame(all_pairs).to_csv(REPORT_DIR / "colocation_pairs.tsv",
                                   sep='\t', index=False)
    print(f"  Wrote colocation_pairs.tsv ({len(all_pairs)} pairs)")

# ── Summary bar chart ─────────────────────────────────────────────────────────
if all_samples:
    names  = [s['sample'] for s in all_samples]
    n_gene = [s['n_genes'] for s in all_samples]
    n_pair = [s['n_pairs'] for s in all_samples]

    x = range(len(names))
    fig, axes = plt.subplots(1, 2, figsize=(max(10, len(names) * 1.5), 5))
    fig.patch.set_facecolor('#1a1a2e')

    for ax, vals, colour, title in [
        (axes[0], n_gene, '#2196F3', 'Resistance Genes per Sample'),
        (axes[1], n_pair, '#4CAF50', 'Co-located Pairs per Sample'),
    ]:
        ax.set_facecolor('#1a1a2e')
        bars = ax.bar(x, vals, color=colour, edgecolor='#4a4a6a', linewidth=0.5)
        ax.set_xticks(list(x))
        ax.set_xticklabels(names, rotation=45, ha='right', color='white', fontsize=8)
        ax.set_title(title, color='white')
        ax.set_ylabel('Count', color='white')
        ax.tick_params(colors='white')
        for spine in ax.spines.values():
            spine.set_edgecolor('#555577')
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                    str(val), ha='center', va='bottom', color='white', fontsize=8)

    plt.tight_layout()
    plt.savefig(REPORT_DIR / "summary_chart.png", dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print("  Wrote summary_chart.png")

# ── HTML report ───────────────────────────────────────────────────────────────
def sample_html(s):
    pairs_html = ""
    if s['pairs']:
        rows = "".join(
            f"""<tr>
                <td>{p['contig']}</td>
                <td class='gene'>{p['gene_1']}</td><td>{p['class_1']}</td>
                <td>{p['start_1']:,}–{p['end_1']:,}</td>
                <td class='gene'>{p['gene_2']}</td><td>{p['class_2']}</td>
                <td>{p['start_2']:,}–{p['end_2']:,}</td>
                <td class='dist'>{p['distance_bp']:,} bp</td>
                <td>{p['source']}</td>
            </tr>"""
            for p in s['pairs']
        )
        pairs_html = f"""
        <h3>Co-located Pairs ({len(s['pairs'])})</h3>
        <table>
          <thead>
            <tr>
              <th>Contig</th><th>Gene 1</th><th>Class 1</th><th>Position 1</th>
              <th>Gene 2</th><th>Class 2</th><th>Position 2</th>
              <th>Distance</th><th>Source</th>
            </tr>
          </thead>
          <tbody>{rows}</tbody>
        </table>"""
    else:
        pairs_html = "<p class='none'>No co-located pairs detected within 50 kb</p>"

    # Map images
    map_dir = REPORT_DIR / "maps" / s['sample']
    maps_html = ""
    if map_dir.exists():
        for img in sorted(map_dir.glob("*.png")):
            rel = img.relative_to(REPORT_DIR)
            maps_html += f'<img src="{rel}" alt="{img.stem}" class="map-img">\n'

    return f"""
    <section class='sample'>
      <h2>{s['sample']}</h2>
      <div class='stats'>
        <span>🧬 Resistance genes: <b>{s['n_genes']}</b></span>
        <span>🔗 Co-located pairs: <b>{s['n_pairs']}</b></span>
        <span>🔄 Integrons: <b>{s['n_integrons']}</b></span>
        <span>📌 IS elements: <b>{s['n_is_elements']}</b></span>
      </div>
      {maps_html}
      {pairs_html}
    </section>"""

total_genes = sum(s['n_genes'] for s in all_samples)
total_pairs = sum(s['n_pairs'] for s in all_samples)

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Resistance Gene Co-location Report</title>
<style>
  :root {{
    --bg: #1a1a2e; --surface: #16213e; --card: #0f3460;
    --accent: #2196F3; --green: #4CAF50; --orange: #FF9800;
    --text: #e0e0e0; --muted: #9e9e9e;
  }}
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: 'Segoe UI', system-ui, sans-serif;
         background: var(--bg); color: var(--text); padding: 2rem; }}
  h1 {{ color: var(--accent); font-size: 1.8rem; margin-bottom: 0.3rem; }}
  .subtitle {{ color: var(--muted); margin-bottom: 2rem; }}
  .overview {{ display: flex; gap: 1rem; flex-wrap: wrap; margin-bottom: 2rem; }}
  .kpi {{ background: var(--card); border-radius: 10px; padding: 1rem 1.5rem;
          min-width: 140px; text-align: center; }}
  .kpi .val {{ font-size: 2rem; font-weight: bold; color: var(--accent); }}
  .kpi .lbl {{ color: var(--muted); font-size: 0.85rem; }}
  .summary-img {{ width: 100%; max-width: 900px; border-radius: 10px;
                  margin-bottom: 2rem; }}
  .sample {{ background: var(--surface); border-radius: 12px;
             padding: 1.5rem; margin-bottom: 2rem;
             border-left: 4px solid var(--accent); }}
  h2 {{ color: var(--accent); font-size: 1.2rem; margin-bottom: 0.75rem; }}
  h3 {{ color: var(--green); margin: 1rem 0 0.5rem; }}
  .stats {{ display: flex; gap: 1.5rem; flex-wrap: wrap;
            color: var(--muted); font-size: 0.9rem; margin-bottom: 1rem; }}
  .stats b {{ color: var(--text); }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.82rem;
           margin-top: 0.5rem; }}
  th {{ background: var(--card); color: var(--muted); text-align: left;
        padding: 0.5rem 0.75rem; border-bottom: 1px solid #333; }}
  td {{ padding: 0.4rem 0.75rem; border-bottom: 1px solid #222; }}
  tr:hover td {{ background: rgba(255,255,255,0.04); }}
  .gene {{ color: #81d4fa; font-family: monospace; }}
  .dist {{ color: var(--orange); }}
  .none {{ color: var(--muted); font-style: italic; padding: 0.5rem 0; }}
  .map-img {{ width: 100%; max-width: 900px; border-radius: 8px;
              margin: 0.5rem 0; display: block; }}
  footer {{ color: var(--muted); font-size: 0.8rem; margin-top: 2rem;
            text-align: center; }}
</style>
</head>
<body>
<h1>🧬 Resistance Gene Co-location Report</h1>
<p class="subtitle">Generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')} |
  Tools: AMRFinderPlus · RGI/CARD · Bakta · IntegronFinder ·
  PlasmidFinder · MOB-suite</p>

<div class="overview">
  <div class="kpi"><div class="val">{len(all_samples)}</div><div class="lbl">Samples</div></div>
  <div class="kpi"><div class="val">{total_genes}</div><div class="lbl">Resistance Genes</div></div>
  <div class="kpi"><div class="val">{total_pairs}</div><div class="lbl">Co-located Pairs</div></div>
</div>

<img src="summary_chart.png" class="summary-img" alt="Summary chart">

{''.join(sample_html(s) for s in all_samples)}

<footer>
  Pipeline: AMRFinderPlus · RGI (CARD) · Bakta · IntegronFinder ·
  PlasmidFinder · MOB-suite<br>
  Co-location threshold: ≤ 50,000 bp on same contig
</footer>
</body>
</html>"""

report_path = REPORT_DIR / "report.html"
report_path.write_text(html)
print(f"\n✔ Report written to: {report_path}")
print(f"  Open with: open '{report_path}'  (Mac) or xdg-open '{report_path}'  (Linux)")
PYEOF

    _conda_activate "$CONDA_ENV"
    python "$REPORT_SCRIPT" "$OUTPUT_DIR/results"
}

# =============================================================================
# MAIN
# =============================================================================
main() {
    echo ""
    echo -e "${BLUE}╔══════════════════════════════════════════════════╗${NC}"
    echo -e "${BLUE}║   Resistance Gene Co-location Analysis Pipeline  ║${NC}"
    echo -e "${BLUE}╚══════════════════════════════════════════════════╝${NC}"
    echo ""
    info "Contig directory : $CONTIG_DIR"
    info "Output directory : $OUTPUT_DIR"
    info "Threads          : $THREADS"
    info "Platform         : $PLATFORM"
    info "Contigs found    : $CONTIG_COUNT"
    echo ""

    mkdir -p "$OUTPUT_DIR"

    install_conda
    create_environment
    setup_databases
    run_analysis
    generate_report

    section "PIPELINE COMPLETE"
    log "Results : $OUTPUT_DIR/results/"
    log "Report  : $OUTPUT_DIR/report/report.html"
    echo ""
    echo -e "${GREEN}Open the report with:${NC}"
    if [[ "$PLATFORM" == "mac" ]]; then
        echo "  open '$OUTPUT_DIR/report/report.html'"
    else
        echo "  xdg-open '$OUTPUT_DIR/report/report.html'"
    fi
    echo ""
}

main "$@"
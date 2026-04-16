"""
alignment_visualization.py  (Windows-compatible version)
=========================================================
pysam is NOT available on Windows.  This module replaces every pysam call
with a lightweight, pure-Python in-memory alignment store so the rest of
the code works identically on Windows, macOS, and Linux.
"""

import random
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines
from theoretical_stats import compute_alpha   # local module – no change needed


# ---------------------------------------------------------------------------
# Lightweight in-memory alignment store (pysam replacement)
# ---------------------------------------------------------------------------

class AlignedRead:
    """Minimal stand-in for pysam.AlignedSegment."""
    __slots__ = ("query_name", "query_sequence",
                 "reference_start", "reference_end", "mapping_quality")

    def __init__(self, name, seq, start, mapq):
        self.query_name      = name
        self.query_sequence  = seq
        self.reference_start = start
        self.reference_end   = start + len(seq)
        self.mapping_quality = mapq


class InMemoryBAM:
    """Collects AlignedRead objects the same way pysam would write them to a file."""

    def __init__(self):
        self._reads = []

    def write(self, read: AlignedRead):
        self._reads.append(read)

    def fetch(self):
        return iter(self._reads)

    def __len__(self):
        return len(self._reads)


# Module-level cache: maps bam_filename -> InMemoryBAM
# This lets create_alignment_bam() "save" a BAM and other functions "open" it.
_BAM_STORE: dict[str, InMemoryBAM] = {}


# ---------------------------------------------------------------------------
# Core alignment helpers
# ---------------------------------------------------------------------------

def find_with_mismatches(genome: str, read: str, max_mismatch: int = 2) -> list[int]:
    """Return all start positions in genome where read aligns with <= max_mismatch mismatches."""
    positions = []
    L = len(read)
    for i in range(len(genome) - L + 1):
        segment = genome[i: i + L]
        mismatches = sum(1 for a, b in zip(segment, read) if a != b)
        if mismatches <= max_mismatch:
            positions.append(i)
    return positions


def create_alignment_bam(genome: str, reads: list[str],
                         bam_filename: str, max_mismatch: int = 2):
    """
    Realistic alignment pipeline (pysam-free, Windows-compatible).

    Stores results in the module-level _BAM_STORE dict under *bam_filename*.
    Calling code can pass any string as bam_filename; no actual file is written.
    """
    bam = InMemoryBAM()

    for i, read_seq in enumerate(reads):
        positions = find_with_mismatches(genome, read_seq, max_mismatch)

        if not positions:
            continue   # unmapped read

        # Simulate mapping quality
        if len(positions) == 1:
            mapq = 60
        elif len(positions) <= 5:
            mapq = 20
        else:
            mapq = 0   # highly repetitive

        # Drop very ambiguous reads (realistic behaviour)
        if mapq == 0 and random.random() < 0.9:
            continue
        if mapq == 20 and random.random() < 0.8:
            continue

        # Choose one mapping position (like real aligners do)
        pos = random.choice(positions[:5])

        read = AlignedRead(
            name=f"read_{i}",
            seq=read_seq,
            start=pos,
            mapq=mapq,
        )
        bam.write(read)

    _BAM_STORE[bam_filename] = bam


def _open_bam(bam_filename: str) -> InMemoryBAM:
    """Retrieve an in-memory BAM (raises KeyError if not created yet)."""
    if bam_filename not in _BAM_STORE:
        raise KeyError(
            f"BAM '{bam_filename}' not found in memory. "
            "Call create_alignment_bam() first."
        )
    return _BAM_STORE[bam_filename]


def _remove_bam(bam_filename: str):
    """Discard an in-memory BAM (mirrors os.remove() on real files)."""
    _BAM_STORE.pop(bam_filename, None)


# ---------------------------------------------------------------------------
# Contig & coverage extraction
# ---------------------------------------------------------------------------

def get_contigs(reads_data: list[dict]) -> list[tuple[int, int]]:
    """
    Identifies contigs by merging all overlapping reads.

    Args:
        reads_data: list of dicts with keys 'start' and 'end'.

    Returns:
        List of (start, end) tuples representing merged contig boundaries.
    """
    if not reads_data:
        return []

    sorted_reads = sorted(reads_data, key=lambda x: x['start'])
    contigs = []
    current_start = sorted_reads[0]['start']
    current_end   = sorted_reads[0]['end']

    for r in sorted_reads[1:]:
        if r['start'] <= current_end:
            current_end = max(current_end, r['end'])
        else:
            contigs.append((current_start, current_end))
            current_start = r['start']
            current_end   = r['end']

    contigs.append((current_start, current_end))
    return contigs


def get_coverage(genome_length: int, reads_data: list[dict]) -> list[int]:
    """
    Computes read depth at every position in the genome.

    The average coverage is theoretically alpha = N*L/G.
    """
    coverage = [0] * genome_length
    for r in reads_data:
        for pos in range(r['start'], r['end']):
            if 0 <= pos < genome_length:
                coverage[pos] += 1
    return coverage


# ---------------------------------------------------------------------------
# Effective coverage & bootstrap statistics
# ---------------------------------------------------------------------------

def effective_alpha(reads_data: list[dict], genome_len: int,
                    L_mean: float) -> tuple[float, float, float]:
    """
    Computes nominal alpha vs effective alpha (MAPQ-60 reads only).

    A region may have moderate nominal coverage, but if a non-trivial fraction
    comes from multi-mapped reads the effective assembly information is much
    smaller.  alpha_eff only counts uniquely-mapped reads (MAPQ == 60).

    Args:
        reads_data:  list of dicts with at least keys 'start', 'end', 'mapq'.
        genome_len:  genome length G.
        L_mean:      mean read length L used for the alpha formula.

    Returns:
        (alpha_nominal, alpha_effective, mapq_loss)
        mapq_loss = fraction of coverage lost to multi-mapping  (0-1).
    """
    if genome_len == 0 or not reads_data:
        return 0.0, 0.0, 0.0

    unique_reads    = [r for r in reads_data if r.get('mapq', 0) == 60]
    alpha_nominal   = len(reads_data)   * L_mean / genome_len
    alpha_effective = len(unique_reads) * L_mean / genome_len
    mapq_loss       = 1.0 - (alpha_effective / alpha_nominal) if alpha_nominal > 0 else 0.0
    return alpha_nominal, alpha_effective, mapq_loss


def bootstrap_ci(trial_counts: list[int],
                 n_boot: int = 500) -> tuple[float, float]:
    """
    Bootstrap 95 % confidence interval for the mean of *trial_counts*.

    Args:
        trial_counts: list of per-trial contig counts at one alpha value.
        n_boot:       number of bootstrap resamples (default 500).

    Returns:
        (lower_2.5%, upper_97.5%) percentile tuple.
    """
    counts = np.array(trial_counts, dtype=float)
    boots  = [
        np.mean(np.random.choice(counts, size=len(counts), replace=True))
        for _ in range(n_boot)
    ]
    return float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def plot_effective_alpha_comparison(results: list[dict], output_png: str):
    """
    Side-by-side bar chart: nominal alpha vs effective alpha for each condition.

    Args:
        results: list of dicts, each with keys:
                   'label'            – condition name (str)
                   'alpha_nominal'    – float
                   'alpha_effective'  – float
                   'mapq_loss'        – float (fraction, 0-1)
        output_png: filename to save the figure.
    """
    labels    = [r['label']           for r in results]
    nominals  = [r['alpha_nominal']   for r in results]
    effs      = [r['alpha_effective'] for r in results]
    losses    = [r['mapq_loss'] * 100 for r in results]   # convert to %

    x     = np.arange(len(labels))
    width = 0.35

    fig, (ax_bar, ax_loss) = plt.subplots(1, 2, figsize=(14, 5))

    # Left panel – nominal vs effective alpha bars
    bars_nom = ax_bar.bar(x - width / 2, nominals, width,
                          label='Nominal α', color='steelblue', alpha=0.85)
    bars_eff = ax_bar.bar(x + width / 2, effs,     width,
                          label='Effective α (MAPQ-60 only)',
                          color='darkorange', alpha=0.85)

    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels, fontsize=9)
    ax_bar.set_ylabel(r"Coverage depth $\alpha$")
    ax_bar.set_title("Nominal vs Effective α\n(only uniquely-mapped reads contribute)",
                     fontweight='bold')
    ax_bar.legend()
    ax_bar.grid(axis='y', linestyle='--', alpha=0.4)

    # Annotate bars with values
    for bar in bars_nom:
        ax_bar.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.02,
                    f"{bar.get_height():.2f}",
                    ha='center', va='bottom', fontsize=8)
    for bar in bars_eff:
        ax_bar.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.02,
                    f"{bar.get_height():.2f}",
                    ha='center', va='bottom', fontsize=8)

    # Right panel – MAPQ loss %
    colors = ['steelblue' if 'Non-Repeat' in l else 'tomato' for l in labels]
    ax_loss.bar(x, losses, color=colors, alpha=0.85)
    ax_loss.set_xticks(x)
    ax_loss.set_xticklabels(labels, fontsize=9)
    ax_loss.set_ylabel("Coverage lost to multi-mapping (%)")
    ax_loss.set_title("MAPQ Loss: fraction of α lost\nto ambiguous/multi-mapped reads",
                      fontweight='bold')
    ax_loss.grid(axis='y', linestyle='--', alpha=0.4)
    for i, v in enumerate(losses):
        ax_loss.text(i, v + 0.3, f"{v:.1f}%", ha='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(output_png)
    print(f"Effective alpha comparison plot saved to {output_png}")


# ---------------------------------------------------------------------------
# Visualisation helpers
# ---------------------------------------------------------------------------

def plot_contigs(ax, contigs, y=1, height=0.6, color='red',
                 show_labels=True, show_gaps=False):
    """
    Plots contigs as solid horizontal bars on a matplotlib axis.

    Args:
        ax:           matplotlib Axes object.
        contigs:      list of (start, end) tuples.
        y:            vertical position of the contig bars.
        height:       thickness of the contig bars (used for label offset).
        color:        bar colour.
        show_labels:  whether to label each contig (C0, C1, …).
        show_gaps:    whether to draw a dashed line across each gap.
    """
    if not contigs:
        return

    for i, (start, end) in enumerate(contigs):
        ax.plot([start, end], [y, y], color=color,
                linewidth=8, solid_capstyle='butt')
        if show_labels:
            ax.text((start + end) / 2, y + height,
                    f"C{i}", ha='center', fontsize=8, color=color)

    if show_gaps and len(contigs) > 1:
        for i in range(len(contigs) - 1):
            gap_start = contigs[i][1]
            gap_end   = contigs[i + 1][0]
            ax.plot([gap_start, gap_end], [y, y],
                    color='gray', linestyle='dashed', linewidth=2)


def plot_contigs_stacked(ax, contigs, base_y=1, spacing=0.5, color='red'):
    for i, (start, end) in enumerate(contigs):
        y = base_y + i * spacing
        ax.plot([start, end], [y, y], color=color,
                linewidth=6, solid_capstyle='butt')


def plot_lander_waterman_sweep(genome_len: int, read_len: float,
                                sweep_data_nr: list, sweep_data_r: list,
                                output_png: str,
                                ci_nr: list = None, ci_r: list = None):
    """
    Plots the Lander-Waterman sweep comparing theory to empirical results,
    with optional bootstrap 95% CI bands.

    Args:
        genome_len:    Genome size G.
        read_len:      Mean read length L.
        sweep_data_nr: List of (alpha, mean_contigs) for no-repeat genome.
        sweep_data_r:  List of (alpha, mean_contigs) for repeat genome.
        output_png:    Output filename.
        ci_nr:         Optional list of (lo, hi) CI tuples aligned with sweep_data_nr.
        ci_r:          Optional list of (lo, hi) CI tuples aligned with sweep_data_r.
    """
    alphas_nr  = [d[0] for d in sweep_data_nr]
    contigs_nr = [d[1] for d in sweep_data_nr]
    alphas_r   = [d[0] for d in sweep_data_r]
    contigs_r  = [d[1] for d in sweep_data_r]

    # Theoretical curve: E[contigs] = (G/L) * alpha * e^(-alpha)
    alpha_range = np.linspace(0, max(max(alphas_nr), max(alphas_r), 6), 100)
    theory_y    = (genome_len / read_len) * alpha_range * np.exp(-alpha_range)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    fig.suptitle(
        r"Lander-Waterman Model: $E[\# \text{contigs}] = N e^{-\alpha}$",
        fontsize=16
    )

    # Panel 1 – No Repeats
    ax1.plot(alpha_range, theory_y, '--', color='gray',
             label=r"Theory: $(G/L) \alpha e^{-\alpha}$")
    if ci_nr is not None:
        lows_nr  = [c[0] for c in ci_nr]
        highs_nr = [c[1] for c in ci_nr]
        ax1.fill_between(alphas_nr, lows_nr, highs_nr,
                         alpha=0.2, color='C0', label="95% bootstrap CI")
    ax1.stem(alphas_nr, contigs_nr,
             linefmt='lightblue', markerfmt='C0o', basefmt=" ",
             label="Empirical mean (No repeats)")
    ax1.scatter(alphas_nr, contigs_nr, color='C0', s=50, zorder=3)
    ax1.set_title("Lander-Waterman sweep — No repeats", fontweight='bold')
    ax1.set_xlabel(r"Read density $\alpha = NL/G$")
    ax1.set_ylabel("Number of contigs")
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.3)

    typical_alpha  = 3.25
    typical_theory = (genome_len / read_len) * typical_alpha * np.exp(-typical_alpha)
    ax1.annotate(
        f"α = {typical_alpha}\n(typical run)",
        xy=(typical_alpha, typical_theory),
        xytext=(typical_alpha + 0.5, typical_theory + 2),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
        fontsize=9, color='gray'
    )

    # Panel 2 – With Repeats
    ax2.plot(alpha_range, theory_y, '--', color='gray',
             label=r"Theory: $(G/L) \alpha e^{-\alpha}$")
    if ci_r is not None:
        lows_r  = [c[0] for c in ci_r]
        highs_r = [c[1] for c in ci_r]
        ax2.fill_between(alphas_r, lows_r, highs_r,
                         alpha=0.2, color='C3', label="95% bootstrap CI")
    ax2.stem(alphas_r, contigs_r,
             linefmt='mistyrose', markerfmt='C3o', basefmt=" ",
             label="Empirical mean (With repeats)")
    ax2.scatter(alphas_r, contigs_r, color='C3', s=50, zorder=3)
    ax2.set_title("Lander-Waterman sweep — With repeats", fontweight='bold')
    ax2.set_xlabel(r"Read density $\alpha = NL/G$")
    ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.3)

    ax2.annotate(
        f"α = {typical_alpha}\n(typical run)",
        xy=(typical_alpha, typical_theory),
        xytext=(typical_alpha + 0.5, typical_theory + 2),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
        fontsize=9, color='gray'
    )

    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='red')
    ax2.text(
        0.95, 0.95,
        "Dots above the curve:\nPoisson assumption violated\nby repeat structure",
        transform=ax2.transAxes, fontsize=9,
        verticalalignment='top', horizontalalignment='right',
        bbox=props, color='red'
    )

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_png)
    print(f"Lander-Waterman sweep plot saved to {output_png}")


# ---------------------------------------------------------------------------
# Main alignment visualisation pipeline
# ---------------------------------------------------------------------------

def visualize_alignment(genome: str, reads: list[str],
                        output_png: str, repeat_target=None):
    """
    Full pipeline: align reads -> extract contigs/coverage -> two-panel plot.

    Args:
        genome:        Reference genome sequence (G).
        reads:         List of sampled reads (N).
        output_png:    Output filename for the alignment plot.
        repeat_target: Optional – not used in current visualisation.
    """
    bam_key = "_temp_alignment_"
    create_alignment_bam(genome, reads, bam_key)

    bam = _open_bam(bam_key)
    reads_data = []
    for read in bam.fetch():
        reads_data.append({
            'start': read.reference_start,
            'end':   read.reference_end,
            'seq':   read.query_sequence,
            'mapq':  read.mapping_quality,
        })
    _remove_bam(bam_key)

    contigs       = get_contigs(reads_data)
    genome_length = len(genome)
    coverage      = get_coverage(genome_length, reads_data)

    fig, (ax_main, ax_cov) = plt.subplots(
        2, 1, figsize=(12, 10),
        gridspec_kw={'height_ratios': [3, 1]},
        sharex=True
    )

    # Reference bar
    ax_main.axhline(y=0, color='black', linewidth=4, label='Reference Genome')
    ax_main.text(-2, 0, "Reference",
                 verticalalignment='center', fontweight='bold')

    # Contig bar
    plot_contigs(ax_main, contigs, y=1, show_labels=True, show_gaps=True)
    ax_main.text(-2, 1, "Contigs",
                 color='red', verticalalignment='center',
                 fontweight='bold', fontsize=10)

    # Individual reads coloured by MAPQ
    for i, r in enumerate(reads_data):
        y_pos = -(i + 1)
        if r['mapq'] >= 60:
            color = 'steelblue'
        elif r['mapq'] >= 20:
            color = 'orange'
        else:
            color = 'lightgray'
        ax_main.plot([r['start'], r['end']], [y_pos, y_pos],
                     color=color, linewidth=2, alpha=0.8)

    ax_main.set_xlim(-5, genome_length + 5)
    ax_main.set_ylim(-len(reads_data) - 5, 3)
    ax_main.set_title('Shotgun Simulation: Realistic Alignment and Assembly')

    unique_line = mlines.Line2D([], [], color='steelblue',  linewidth=2, label='Unique (MAPQ 60)')
    multi_line  = mlines.Line2D([], [], color='orange',     linewidth=2, label='Multimapper (MAPQ 20)')
    lowq_line   = mlines.Line2D([], [], color='lightgray',  linewidth=2, label='Ambiguous (MAPQ 0)')
    contig_line = mlines.Line2D([], [], color='red',        linewidth=8, label='Contig (Assembled)')
    ref_line    = mlines.Line2D([], [], color='black',      linewidth=4, label='Reference Genome')

    ax_main.legend(
        handles=[ref_line, contig_line, unique_line, multi_line, lowq_line],
        loc='upper right', fontsize='small'
    )
    ax_main.set_yticks([])

    # Coverage track
    x_positions = range(genome_length)
    ax_cov.fill_between(x_positions, coverage, color='skyblue', alpha=0.3)
    ax_cov.plot(x_positions, coverage, color='steelblue', linewidth=1.5)
    ax_cov.set_xlabel('Position in Genome (bp)')
    ax_cov.set_ylabel('Read Depth')
    ax_cov.set_title('Read Coverage Track')
    ax_cov.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    plt.savefig(output_png)
    print(f"Realistic alignment and contig plot saved to {output_png}")


# ---------------------------------------------------------------------------
# Standalone test block
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from shortgun_simulation import generate_genome, generate_reads
    g = generate_genome(100)
    r = generate_reads(g, 10, 10, 30)
    visualize_alignment(g, r, "test_contig_plot.png")

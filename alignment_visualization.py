import random
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
import matplotlib.lines as mlines


# ─────────────────────────────────────────────────────────────────────────────
# Core alignment  (pure Python — no pysam, no BAM files)
# ─────────────────────────────────────────────────────────────────────────────

def find_with_mismatches(genome, read, max_mismatch=2):
    """Return all start positions where `read` aligns to `genome` with at most
    `max_mismatch` mismatches."""
    positions = []
    L = len(read)
    for i in range(len(genome) - L + 1):
        segment = genome[i:i + L]
        mismatches = sum(1 for a, b in zip(segment, read) if a != b)
        if mismatches <= max_mismatch:
            positions.append(i)
    return positions


def align_reads(genome, reads, max_mismatch=2):
    """
    Pure-Python replacement for pysam BAM alignment.
    Aligns reads to the genome with mismatch tolerance and returns a list of
    dicts with 'start', 'end', 'seq', and 'mapq' — same schema used downstream.

    Args:
        genome (str): Reference genome sequence.
        reads (list): List of reads (plain strings or (start, seq) tuples).
        max_mismatch (int): Maximum allowed mismatches for a valid alignment.

    Returns:
        list[dict]: Aligned reads with keys 'start', 'end', 'seq', 'mapq'.
    """
    reads_data = []

    for i, read_item in enumerate(reads):
        read_seq = read_item[1] if isinstance(read_item, tuple) else read_item

        positions = find_with_mismatches(genome, read_seq, max_mismatch)

        if not positions:
            continue  # unmapped read

        # Simulate mapping quality
        if len(positions) == 1:
            mapq = 60
        elif len(positions) <= 5:
            mapq = 20
        else:
            mapq = 0  # highly repetitive

        # Drop very ambiguous reads (realistic behaviour)
        if mapq == 0 and random.random() < 0.9:
            continue
        if mapq == 20 and random.random() < 0.8:
            continue

        # Choose one mapping position (like real aligners do)
        pos = random.choice(positions[:5])

        reads_data.append({
            'start': pos,
            'end':   pos + len(read_seq),
            'seq':   read_seq,
            'mapq':  mapq,
        })

    return reads_data


# ─────────────────────────────────────────────────────────────────────────────
# Contig & coverage helpers
# ─────────────────────────────────────────────────────────────────────────────

def get_contigs(reads_data):
    """
    Identifies contigs by merging all overlapping reads.
    In an ideal Poisson model the expected number of contigs is N * e^{-alpha}.

    Args:
        reads_data (list[dict]): Dicts with 'start' and 'end' positions.

    Returns:
        list[tuple]: (start, end) tuples for each merged contig.
    """
    if not reads_data:
        return []

    sorted_reads = sorted(reads_data, key=lambda x: x['start'])

    current_start = sorted_reads[0]['start']
    current_end   = sorted_reads[0]['end']
    contigs = []

    for r in sorted_reads[1:]:
        if r['start'] <= current_end:
            current_end = max(current_end, r['end'])
        else:
            contigs.append((current_start, current_end))
            current_start = r['start']
            current_end   = r['end']

    contigs.append((current_start, current_end))
    return contigs


def get_coverage(genome_length, reads_data):
    """
    Computes per-position read depth.
    The average coverage is theoretically alpha = NL/G.

    Args:
        genome_length (int): Length of the reference genome.
        reads_data (list[dict]): Aligned reads with 'start' and 'end'.

    Returns:
        list[int]: Coverage depth at each position.
    """
    coverage = [0] * genome_length
    for r in reads_data:
        for pos in range(r['start'], r['end']):
            if 0 <= pos < genome_length:
                coverage[pos] += 1
    return coverage


# ─────────────────────────────────────────────────────────────────────────────
# Visualisation helpers
# ─────────────────────────────────────────────────────────────────────────────

def plot_contigs(ax, contigs, y=1, height=0.6, color='red',
                 show_labels=True, show_gaps=False):
    """
    Plots contigs as thick horizontal bars on a matplotlib axis.

    Args:
        ax: Matplotlib axis.
        contigs (list): List of (start, end) tuples.
        y (float): Vertical position.
        height (float): Label offset above bar.
        color (str): Bar colour.
        show_labels (bool): Label each contig C0, C1, ...
        show_gaps (bool): Draw dashed lines across gaps between contigs.
    """
    if not contigs:
        return

    for i, (start, end) in enumerate(contigs):
        ax.plot([start, end], [y, y], color=color, linewidth=8,
                solid_capstyle='butt')
        if show_labels:
            ax.text((start + end) / 2, y + height, f"C{i}",
                    ha='center', fontsize=8, color=color)

    if show_gaps and len(contigs) > 1:
        for i in range(len(contigs) - 1):
            gap_start = contigs[i][1]
            gap_end   = contigs[i + 1][0]
            ax.plot([gap_start, gap_end], [y, y],
                    color='gray', linestyle='dashed', linewidth=2)


def plot_contigs_stacked(ax, contigs, base_y=1, spacing=0.5, color='red'):
    for i, (start, end) in enumerate(contigs):
        y = base_y + i * spacing
        ax.plot([start, end], [y, y], color=color, linewidth=6,
                solid_capstyle='butt')


# ─────────────────────────────────────────────────────────────────────────────
# Full alignment visualisation
# ─────────────────────────────────────────────────────────────────────────────

def visualize_alignment(genome, reads, output_png, repeat_target=None):
    """
    Aligns reads to the genome (in memory) and produces a two-panel plot:
      - Top:    individual reads coloured by mapping quality + merged contigs.
      - Bottom: read coverage (depth) track.

    Args:
        genome (str): Reference genome sequence.
        reads (list): Sampled reads (strings or (start, seq) tuples).
        output_png (str): Output filename.
        repeat_target (str, optional): Unused; kept for API compatibility.
    """
    reads_data    = align_reads(genome, reads)
    contigs       = get_contigs(reads_data)
    genome_length = len(genome)
    coverage      = get_coverage(genome_length, reads_data)

    fig, (ax_main, ax_cov) = plt.subplots(
        2, 1, figsize=(12, 10),
        gridspec_kw={'height_ratios': [3, 1]}, sharex=True
    )

    # 1. Reference line
    ax_main.axhline(y=0, color='black', linewidth=4)
    ax_main.text(-2, 0, "Reference", verticalalignment='center', fontweight='bold')

    # 2. Contigs
    plot_contigs(ax_main, contigs, y=1, show_labels=True, show_gaps=True)
    ax_main.text(-2, 1, "Contigs", color='red',
                 verticalalignment='center', fontweight='bold', fontsize=10)

    # 3. Individual reads coloured by MAPQ
    for i, r in enumerate(reads_data):
        if r['mapq'] >= 60:
            color = 'steelblue'
        elif r['mapq'] >= 20:
            color = 'orange'
        else:
            color = 'lightgray'
        ax_main.plot([r['start'], r['end']], [-(i + 1), -(i + 1)],
                     color=color, linewidth=2, alpha=0.8)

    ax_main.set_xlim(-5, genome_length + 5)
    ax_main.set_ylim(-len(reads_data) - 5, 3)
    ax_main.set_title('Shotgun Simulation: Realistic Alignment and Assembly')
    ax_main.set_yticks([])

    legend_handles = [
        mlines.Line2D([], [], color='black',     linewidth=4, label='Reference Genome'),
        mlines.Line2D([], [], color='red',       linewidth=8, label='Contig (Assembled)'),
        mlines.Line2D([], [], color='steelblue', linewidth=2, label='Unique (MAPQ 60)'),
        mlines.Line2D([], [], color='orange',    linewidth=2, label='Multimapper (MAPQ 20)'),
        mlines.Line2D([], [], color='lightgray', linewidth=2, label='Ambiguous (MAPQ 0)'),
    ]
    ax_main.legend(handles=legend_handles, loc='upper right', fontsize='small')

    # 4. Coverage track
    x_pos = range(genome_length)
    ax_cov.fill_between(x_pos, coverage, color='skyblue', alpha=0.3)
    ax_cov.plot(x_pos, coverage, color='steelblue', linewidth=1.5)
    ax_cov.set_xlabel('Position in Genome (bp)')
    ax_cov.set_ylabel('Read Depth')
    ax_cov.set_title('Read Coverage Track')
    ax_cov.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    print(f"Realistic alignment and contig plot saved to {output_png}")


# ─────────────────────────────────────────────────────────────────────────────
# Lander-Waterman sweep plot  (with 95 % bootstrap CI)
# ─────────────────────────────────────────────────────────────────────────────

def _bootstrap_ci(trials, n_boot=1000, ci=95):
    """
    Compute a bootstrap confidence interval for the mean of trial values.

    Args:
        trials (list[float]): Raw per-trial contig counts.
        n_boot (int): Number of bootstrap resamples.
        ci (float): Confidence level (95 = 95 % CI).

    Returns:
        tuple: (mean, lower_bound, upper_bound)
    """
    if not trials:
        return 0.0, 0.0, 0.0
    n = len(trials)
    boot_means = sorted(
        sum(trials[int(random.random() * n)] for _ in range(n)) / n
        for _ in range(n_boot)
    )
    lo_idx = int((1 - ci / 100) / 2 * n_boot)
    hi_idx = min(int((1 + ci / 100) / 2 * n_boot), n_boot - 1)
    return sum(trials) / n, boot_means[lo_idx], boot_means[hi_idx]


def plot_lander_waterman_sweep(genome_len, read_len, sweep_data_nr, sweep_data_r, output_png):
    """
    Plots the Lander-Waterman sweep with 95 % bootstrap CI shading.

    Each entry in sweep_data_nr / sweep_data_r must be one of:
      - (alpha, [trial_count_1, trial_count_2, ...])   <- preferred, CI computed
      - (alpha, mean_float)                             <- legacy, no CI band

    Args:
        genome_len (int):   Genome size G.
        read_len (float):   Mean read length L.
        sweep_data_nr (list): Per-alpha data for the non-repeat genome.
        sweep_data_r  (list): Per-alpha data for the repeat genome.
        output_png (str):   Output filename.
    """
    def unpack(sweep_data):
        alphas, means, lows, highs = [], [], [], []
        for alpha, payload in sweep_data:
            alphas.append(alpha)
            if isinstance(payload, list):
                m, lo, hi = _bootstrap_ci(payload)
            else:
                m, lo, hi = payload, payload, payload   # legacy scalar
            means.append(m)
            lows.append(lo)
            highs.append(hi)
        return alphas, means, lows, highs

    alphas_nr, means_nr, lows_nr, highs_nr = unpack(sweep_data_nr)
    alphas_r,  means_r,  lows_r,  highs_r  = unpack(sweep_data_r)

    alpha_max   = max(max(alphas_nr), max(alphas_r), 6)
    alpha_range = np.linspace(0, alpha_max, 200)
    theory_y    = (genome_len / read_len) * alpha_range * np.exp(-alpha_range)

    typical_alpha  = 3.25
    typical_theory = (genome_len / read_len) * typical_alpha * np.exp(-typical_alpha)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    fig.suptitle(
        r"Lander-Waterman Model: $E[\#\mathrm{contigs}] = Ne^{-\alpha}$",
        fontsize=16
    )

    # ── Panel 1: No Repeats ──────────────────────────────────────────────────
    ax1.plot(alpha_range, theory_y, '--', color='gray', linewidth=1.5,
             label=r"Theory: $(G/L)\,\alpha e^{-\alpha}$")
    ax1.fill_between(alphas_nr, lows_nr, highs_nr,
                     color='steelblue', alpha=0.25, label='95% bootstrap CI')
    for a, m in zip(alphas_nr, means_nr):
        ax1.plot([a, a], [0, m], color='lightblue', linewidth=1, zorder=1)
    ax1.scatter(alphas_nr, means_nr, color='C0', s=50, zorder=3,
                label='Empirical mean (No repeats)')
    ax1.set_title("Lander-Waterman sweep — No repeats", fontweight='bold')
    ax1.set_xlabel(r"Read density $\alpha = NL/G$")
    ax1.set_ylabel("Number of contigs")
    ax1.legend(fontsize=8)
    ax1.grid(True, linestyle='--', alpha=0.3)
    ax1.annotate(
        f"α = {typical_alpha}\n(typical run)",
        xy=(typical_alpha, typical_theory),
        xytext=(typical_alpha + 0.5, typical_theory + 2),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
        fontsize=9, color='gray'
    )

    # ── Panel 2: With Repeats ────────────────────────────────────────────────
    ax2.plot(alpha_range, theory_y, '--', color='gray', linewidth=1.5,
             label=r"Theory: $(G/L)\,\alpha e^{-\alpha}$")
    ax2.fill_between(alphas_r, lows_r, highs_r,
                     color='salmon', alpha=0.30, label='95% bootstrap CI')
    for a, m in zip(alphas_r, means_r):
        ax2.plot([a, a], [0, m], color='mistyrose', linewidth=1, zorder=1)
    ax2.scatter(alphas_r, means_r, color='C3', s=50, zorder=3,
                label='Empirical mean (With repeats)')
    ax2.set_title("Lander-Waterman sweep — With repeats", fontweight='bold')
    ax2.set_xlabel(r"Read density $\alpha = NL/G$")
    ax2.legend(fontsize=8)
    ax2.grid(True, linestyle='--', alpha=0.3)
    ax2.annotate(
        f"α = {typical_alpha}\n(typical run)",
        xy=(typical_alpha, typical_theory),
        xytext=(typical_alpha + 0.5, typical_theory + 2),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
        fontsize=9, color='gray'
    )
    props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='red')
    ax2.text(0.97, 0.97,
             "Dots above the curve:\nPoisson assumption violated\nby repeat structure",
             transform=ax2.transAxes, fontsize=9,
             verticalalignment='top', horizontalalignment='right',
             bbox=props, color='red')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_png, dpi=150)
    print(f"Lander-Waterman sweep plot saved to {output_png}")


# ─────────────────────────────────────────────────────────────────────────────
# Quick self-test
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    from shortgun_simulation import generate_genome, generate_reads
    g = generate_genome(100)
    r = generate_reads(g, 10, 10, 30)
    visualize_alignment(g, r, "test_contig_plot.png")

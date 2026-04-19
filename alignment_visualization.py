"""
alignment_visualization.py
===========================
Pure-Python shotgun alignment and visualization pipeline.
No pysam / BAM dependency — works on Windows, macOS, and Linux.

Alignment results are kept as plain Python dicts:
    {'start': int, 'end': int, 'seq': str, 'mapq': int}
"""

import random
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np


# ---------------------------------------------------------------------------
# Alignment helpers (replaces pysam BAM pipeline)
# ---------------------------------------------------------------------------

def find_with_mismatches(genome, read, max_mismatch=2):
    """Return all start positions where *read* aligns to *genome* with at
    most *max_mismatch* substitutions."""
    positions = []
    L = len(read)
    for i in range(len(genome) - L + 1):
        mismatches = sum(1 for a, b in zip(genome[i:i + L], read) if a != b)
        if mismatches <= max_mismatch:
            positions.append(i)
    return positions


def align_reads(genome, reads, max_mismatch=2):
    """
    Align a list of reads to *genome* and return a list of alignment dicts.

    Replaces the previous pysam / BAM pipeline entirely.
    Mapping-quality logic mirrors the original BAM approach:
      - 1 hit  → MAPQ 60  (unique)
      - 2-5    → MAPQ 20  (low-confidence multi-mapper, 80 % dropped)
      - >5     → MAPQ 0   (highly repetitive,            90 % dropped)

    Args:
        genome (str):  Reference sequence.
        reads  (list): Strings or (start, seq) tuples produced by generate_reads*.
        max_mismatch (int): Allowed substitutions per alignment.

    Returns:
        list[dict]: Alignment records with keys start, end, seq, mapq.
    """
    aligned = []
    for i, read_item in enumerate(reads):
        read_seq = read_item[1] if isinstance(read_item, tuple) else read_item

        positions = find_with_mismatches(genome, read_seq, max_mismatch)
        if not positions:
            continue                          # unmapped

        n_hits = len(positions)
        if n_hits == 1:
            mapq = 60
        elif n_hits <= 5:
            mapq = 20
        else:
            mapq = 0

        # Simulate realistic aligner dropping of ambiguous reads
        if mapq == 0 and random.random() < 0.9:
            continue
        if mapq == 20 and random.random() < 0.8:
            continue

        pos = random.choice(positions[:5])
        aligned.append({
            'start': pos,
            'end':   pos + len(read_seq),
            'seq':   read_seq,
            'mapq':  mapq,
        })

    return aligned


# ---------------------------------------------------------------------------
# Contig / coverage utilities
# ---------------------------------------------------------------------------

def get_contigs(reads_data):
    """
    Merge overlapping reads into contigs via a sweep-line algorithm.

    Args:
        reads_data (list[dict]): Alignment records with 'start' and 'end'.

    Returns:
        list[tuple]: (start, end) pairs for each merged contig.
    """
    if not reads_data:
        return []

    sorted_reads = sorted(reads_data, key=lambda x: x['start'])
    contigs = []
    cur_start = sorted_reads[0]['start']
    cur_end   = sorted_reads[0]['end']

    for r in sorted_reads[1:]:
        if r['start'] <= cur_end:
            cur_end = max(cur_end, r['end'])
        else:
            contigs.append((cur_start, cur_end))
            cur_start = r['start']
            cur_end   = r['end']

    contigs.append((cur_start, cur_end))
    return contigs


def get_coverage(genome_length, reads_data):
    """
    Compute per-base read depth.

    Returns:
        list[int]: Coverage value at each position 0 … genome_length-1.
    """
    coverage = [0] * genome_length
    for r in reads_data:
        for pos in range(r['start'], min(r['end'], genome_length)):
            coverage[pos] += 1
    return coverage


# ---------------------------------------------------------------------------
# Contig plotting helpers
# ---------------------------------------------------------------------------

def plot_contigs(ax, contigs, y=1, height=0.6, color='red',
                 show_labels=True, show_gaps=False):
    """Draw contigs as thick horizontal bars on *ax*."""
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
            ax.plot([contigs[i][1], contigs[i + 1][0]], [y, y],
                    color='gray', linestyle='dashed', linewidth=2)


def plot_contigs_stacked(ax, contigs, base_y=1, spacing=0.5, color='red'):
    for i, (start, end) in enumerate(contigs):
        y = base_y + i * spacing
        ax.plot([start, end], [y, y], color=color, linewidth=6,
                solid_capstyle='butt')


# ---------------------------------------------------------------------------
# Main alignment visualisation (replaces the old BAM-based pipeline)
# ---------------------------------------------------------------------------

def visualize_alignment(genome, reads, output_png, repeat_target=None):
    """
    Align *reads* to *genome*, discover contigs, and save a two-panel plot.

    Panel 1 – Individual reads (coloured by MAPQ) + merged contigs.
    Panel 2 – Per-base read-depth track.

    Args:
        genome      (str):  Reference genome sequence.
        reads       (list): Reads as strings or (start, seq) tuples.
        output_png  (str):  Output image path.
        repeat_target (str, optional): Unused; kept for API compatibility.
    """
    reads_data = align_reads(genome, reads)

    contigs       = get_contigs(reads_data)
    genome_length = len(genome)
    coverage      = get_coverage(genome_length, reads_data)

    fig, (ax_main, ax_cov) = plt.subplots(
        2, 1, figsize=(12, 10),
        gridspec_kw={'height_ratios': [3, 1]},
        sharex=True,
    )

    # ── Reference bar ────────────────────────────────────────────────────────
    ax_main.axhline(y=0, color='black', linewidth=4)
    ax_main.text(-2, 0, "Reference", verticalalignment='center', fontweight='bold')

    # ── Contigs ──────────────────────────────────────────────────────────────
    plot_contigs(ax_main, contigs, y=1, show_labels=True, show_gaps=True)
    ax_main.text(-2, 1, "Contigs", color='red', verticalalignment='center',
                 fontweight='bold', fontsize=10)

    # ── Individual reads coloured by MAPQ ────────────────────────────────────
    for i, r in enumerate(reads_data):
        color = ('steelblue' if r['mapq'] >= 60
                 else 'orange' if r['mapq'] >= 20
                 else 'lightgray')
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

    # ── Coverage track ───────────────────────────────────────────────────────
    x = range(genome_length)
    ax_cov.fill_between(x, coverage, color='skyblue', alpha=0.3)
    ax_cov.plot(x, coverage, color='steelblue', linewidth=1.5)
    ax_cov.set_xlabel('Position in Genome (bp)')
    ax_cov.set_ylabel('Read Depth')
    ax_cov.set_title('Read Coverage Track')
    ax_cov.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout()
    plt.savefig(output_png)
    plt.close(fig)
    print(f"Alignment plot saved to {output_png}")


# ---------------------------------------------------------------------------
# Lander-Waterman sweep plot
# ---------------------------------------------------------------------------

def plot_lander_waterman_sweep(genome_len, read_len, sweep_data_nr, sweep_data_r, output_png):
    """
    Plot the Lander-Waterman sweep with 95 % CI shading on both panels.

    Args:
        genome_len    (int):   Genome size G.
        read_len      (float): Mean read length L.
        sweep_data_nr (list):  [(alpha, mean_contigs, ci_95), …] — no-repeat genome.
        sweep_data_r  (list):  [(alpha, mean_contigs, ci_95), …] — repeat genome.
        output_png    (str):   Output filename.
    """
    alphas_nr  = np.array([d[0] for d in sweep_data_nr])
    contigs_nr = np.array([d[1] for d in sweep_data_nr])
    ci_nr      = np.array([d[2] if len(d) > 2 else 0 for d in sweep_data_nr])

    alphas_r  = np.array([d[0] for d in sweep_data_r])
    contigs_r = np.array([d[1] for d in sweep_data_r])
    ci_r      = np.array([d[2] if len(d) > 2 else 0 for d in sweep_data_r])

    alpha_max   = max(alphas_nr.max(), alphas_r.max(), 6.0)
    alpha_range = np.linspace(0, alpha_max, 300)
    theory_y    = (genome_len / read_len) * alpha_range * np.exp(-alpha_range)

    typical_alpha  = 3.25
    typical_theory = (genome_len / read_len) * typical_alpha * np.exp(-typical_alpha)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
    fig.suptitle(
        r"Lander-Waterman Model: $E[\#\mathrm{contigs}] = N e^{-\alpha}$",
        fontsize=16,
    )

    # ── Panel 1: No Repeats ──────────────────────────────────────────────────
    ax1.plot(alpha_range, theory_y, '--', color='gray', linewidth=1.5,
             label=r"Theory: $(G/L)\,\alpha\,e^{-\alpha}$", zorder=1)

    ax1.fill_between(
        alphas_nr,
        np.maximum(0, contigs_nr - ci_nr),
        contigs_nr + ci_nr,
        color='steelblue', alpha=0.25,
        label="95% bootstrap CI",
    )
    ax1.scatter(alphas_nr, contigs_nr, color='C0', s=50, zorder=4,
                label="Empirical mean (No repeats)")

    ax1.set_title("Lander-Waterman sweep \u2014 No repeats", fontweight='bold')
    ax1.set_xlabel(r"Read density $\alpha = NL/G$")
    ax1.set_ylabel("Number of contigs")
    ax1.legend(fontsize=9)
    ax1.grid(True, linestyle='--', alpha=0.3)
    ax1.annotate(
        f"\u03b1 = {typical_alpha}\n(typical run)",
        xy=(typical_alpha, typical_theory),
        xytext=(typical_alpha + 0.5, typical_theory + 2),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
        fontsize=9, color='gray',
    )

    # ── Panel 2: With Repeats ────────────────────────────────────────────────
    ax2.plot(alpha_range, theory_y, '--', color='gray', linewidth=1.5,
             label=r"Theory: $(G/L)\,\alpha\,e^{-\alpha}$", zorder=1)

    ax2.fill_between(
        alphas_r,
        np.maximum(0, contigs_r - ci_r),
        contigs_r + ci_r,
        color='salmon', alpha=0.35,
        label="95% bootstrap CI",
    )
    ax2.scatter(alphas_r, contigs_r, color='C3', s=50, zorder=4,
                label="Empirical mean (With repeats)")

    ax2.set_title("Lander-Waterman sweep \u2014 With repeats", fontweight='bold')
    ax2.set_xlabel(r"Read density $\alpha = NL/G$")
    ax2.legend(fontsize=9)
    ax2.grid(True, linestyle='--', alpha=0.3)
    ax2.annotate(
        f"\u03b1 = {typical_alpha}\n(typical run)",
        xy=(typical_alpha, typical_theory),
        xytext=(typical_alpha + 0.5, typical_theory + 2),
        arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.2"),
        fontsize=9, color='gray',
    )

    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='red')
    ax2.text(
        0.95, 0.95,
        "Dots above the curve:\nPoisson assumption violated\nby repeat structure",
        transform=ax2.transAxes, fontsize=9,
        verticalalignment='top', horizontalalignment='right',
        bbox=props, color='red',
    )

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(output_png)
    plt.close(fig)
    print(f"Lander-Waterman sweep plot saved to {output_png}")


# ---------------------------------------------------------------------------
# Quick self-test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    from shortgun_simulation import generate_genome, generate_reads
    g = generate_genome(200)
    r = generate_reads(g, 20, 20, 5)
    visualize_alignment(g, r, "test_contig_plot.png")
    print("Self-test complete.")

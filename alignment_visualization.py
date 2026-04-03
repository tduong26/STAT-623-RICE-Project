"""
alignment_visualization.py
---------------------------
pysam-free alignment, contig detection, coverage calculation, and plotting.

pysam has been removed (no Windows / Colab binary available).
All alignment bookkeeping uses plain Python dicts and lists — functionally
identical to the original BAM-based approach.
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection


# ---------------------------------------------------------------------------
# Core alignment & assembly primitives
# ---------------------------------------------------------------------------

def create_alignment(genome: str, reads: list) -> list:
    """Align reads to the genome using exact string matching.

    Args:
        genome: Reference genome sequence.
        reads:  List of shotgun read sequences.
    Returns:
        List of dicts with keys 'start', 'end', 'seq' for each mapped read.
    """
    aligned = []
    for read_seq in reads:
        pos = genome.find(read_seq)
        if pos != -1:
            aligned.append({
                "start": pos,
                "end":   pos + len(read_seq),
                "seq":   read_seq,
            })
    return aligned


def get_contigs(reads_data: list) -> list:
    """Identify contigs by merging all overlapping / adjacent reads.

    In an ideal Poisson model the expected number of contigs is N · e^{-α}.

    Args:
        reads_data: List of dicts with 'start' and 'end' keys.
    Returns:
        List of (start, end) tuples representing merged contig boundaries.
    """
    if not reads_data:
        return []

    sorted_reads  = sorted(reads_data, key=lambda x: x["start"])
    current_start = sorted_reads[0]["start"]
    current_end   = sorted_reads[0]["end"]
    contigs       = []

    for r in sorted_reads[1:]:
        if r["start"] <= current_end:
            current_end = max(current_end, r["end"])
        else:
            contigs.append((current_start, current_end))
            current_start = r["start"]
            current_end   = r["end"]

    contigs.append((current_start, current_end))
    return contigs


def get_coverage(genome_length: int, reads_data: list) -> list:
    """Compute per-position read depth across the genome.

    Average coverage is theoretically α = N·L/G.

    Args:
        genome_length: Total length of the reference genome.
        reads_data:    List of dicts with 'start' and 'end' keys.
    Returns:
        List of ints (depth at each position), length == genome_length.
    """
    coverage = [0] * genome_length
    for r in reads_data:
        for pos in range(r["start"], r["end"]):
            if 0 <= pos < genome_length:
                coverage[pos] += 1
    return coverage


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def visualize_alignment(genome: str, reads: list, output_png: str,
                         repeat_target: str = None) -> None:
    """Create and save a two-panel alignment + coverage plot.

    Panel 1 – Individual reads and merged contigs coloured by depth.
    Panel 2 – Read coverage (depth) track.

    Args:
        genome:        Reference genome sequence.
        reads:         List of sampled read sequences.
        output_png:    Output filename for the plot.
        repeat_target: Optional repeat sequence (accepted but not plotted).
    """
    reads_data = create_alignment(genome, reads)
    contigs    = get_contigs(reads_data)
    genome_len = len(genome)
    coverage   = get_coverage(genome_len, reads_data)

    fig, (ax_main, ax_cov) = plt.subplots(
        2, 1, figsize=(12, 10),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )

    # Reference backbone
    ax_main.axhline(y=0, color="black", linewidth=4)
    ax_main.text(-2, 0, "Reference", verticalalignment="center",
                 fontweight="bold")

    # Contigs coloured by depth using a LineCollection
    max_cov  = max(coverage) if coverage and max(coverage) > 0 else 1
    points   = np.array([np.arange(genome_len),
                          np.ones(genome_len)]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    for start, end in contigs:
        end = min(end, genome_len)
        if end <= start:
            continue
        lc = LineCollection(segments[start:end], cmap="Reds",
                             norm=plt.Normalize(0, max_cov))
        lc.set_array(np.array(coverage[start:end]))
        lc.set_linewidth(6)
        ax_main.add_collection(lc)

    ax_main.text(-2, 1, "Contig(s)", color="red",
                 verticalalignment="center", fontweight="bold")

    # Individual reads
    for i, r in enumerate(reads_data):
        ax_main.plot([r["start"], r["end"]], [-(i + 1), -(i + 1)],
                     color="skyblue", linewidth=2)

    ax_main.set_xlim(-5, genome_len + 5)
    ax_main.set_ylim(-len(reads_data) - 2, 3)
    ax_main.set_title("Shotgun Simulation: Alignment and Contig Discovery")
    ax_main.set_yticks([])

    # Coverage track
    x_pos = range(genome_len)
    ax_cov.fill_between(x_pos, coverage, color="orange", alpha=0.5)
    ax_cov.plot(x_pos, coverage, color="darkorange", linewidth=1.5)
    ax_cov.set_xlabel("Position in Genome (bp)")
    ax_cov.set_ylabel("Read Depth")
    ax_cov.set_title("Read Coverage Track")
    ax_cov.grid(True, linestyle="--", alpha=0.6)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    plt.close(fig)
    print(f"Alignment plot saved to '{output_png}'")


def plot_contigs(ax, contigs: list, y: float = 1, height: float = 0.6,
                 color: str = "red", show_labels: bool = True,
                 show_gaps: bool = False) -> None:
    """Draw contigs as solid horizontal bars on a matplotlib Axes.

    Args:
        ax:          Matplotlib Axes object.
        contigs:     List of (start, end) tuples.
        y:           Vertical position of the bars.
        height:      Vertical offset for contig labels above the bar.
        color:       Bar colour.
        show_labels: Label each contig C0, C1, …
        show_gaps:   Draw dashed lines between consecutive contigs.
    """
    if not contigs:
        return

    for i, (start, end) in enumerate(contigs):
        ax.plot([start, end], [y, y], color=color, linewidth=8,
                solid_capstyle="butt")
        if show_labels:
            ax.text((start + end) / 2, y + height, f"C{i}",
                    ha="center", fontsize=8, color=color)

    if show_gaps and len(contigs) > 1:
        for i in range(len(contigs) - 1):
            gap_start = contigs[i][1]
            gap_end   = contigs[i + 1][0]
            ax.plot([gap_start, gap_end], [y, y],
                    color="gray", linestyle="dashed", linewidth=2)

"""
shortgun_simulation.py
----------------------
Main entry point for the STAT 623 Shotgun Sequencing Simulation.

Run with:
    python shortgun_simulation.py

Two simulations are executed back-to-back:
  1. Random genome (no repeats)  — Lander-Waterman Poisson model holds
  2. Repeat-rich genome          — Poisson assumption is violated

Both produce alignment plots, contig maps, and an overlap graph.
All results are then combined into two summary figures:
  • summary_comparison.png      (6-panel overview)
  • lander_waterman_sweep.png   (standalone α-sweep, high resolution)
  • extended_summary.png        (5 recommended statistical plots)
"""

import matplotlib.pyplot as plt

from alignment_visualization import (
    create_alignment, get_contigs, get_coverage,
    visualize_alignment, plot_contigs,
)
from extended_plots import (
    generate_extended_summary_plots,
    plot_lander_waterman_sweep,
)
from genome_utils import (
    generate_genome, generate_genome_with_repeats, generate_reads,
)
from overlap import build_overlap_graph, draw_graph
from summary_plots import generate_summary_plots
from theoretical_stats import (
    compute_alpha, expected_contigs, expected_contig_length,
)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _collect_stats(genome: str, reads: list) -> dict:
    """Align reads and collect all metrics needed for plots and printing.

    Returns a dict with keys:
        genome, reads_data, contigs, coverage,
        N, G, L, alpha, exp_contigs, exp_contig_len,
        mapped_reads, contig_lengths, read_lengths
    """
    reads_data = create_alignment(genome, reads)
    contigs    = get_contigs(reads_data)
    coverage   = get_coverage(len(genome), reads_data)

    N     = len(reads)
    G     = len(genome)
    L     = max(len(r) for r in reads) if reads else 0
    alpha = compute_alpha(N, L, G)

    return {
        "genome":         genome,
        "reads_data":     reads_data,
        "contigs":        contigs,
        "coverage":       coverage,
        "N":              N,
        "G":              G,
        "L":              L,
        "alpha":          alpha,
        "exp_contigs":    expected_contigs(N, alpha),
        "exp_contig_len": expected_contig_length(G, N, alpha),
        "mapped_reads":   len(reads_data),
        "contig_lengths": [e - s for s, e in contigs],
        "read_lengths":   [len(r) for r in reads],
    }


def _print_stats(label: str, s: dict) -> None:
    """Pretty-print simulation statistics to stdout."""
    print(f"\n{'='*54}")
    print(f"  {label}")
    print(f"{'='*54}")
    print(f"  Genome length    : {s['G']} bp")
    print(f"  Number of reads  : {s['N']}")
    print(f"  Mapped reads     : {s['mapped_reads']}  "
          f"({100*s['mapped_reads']/s['N']:.1f}%)")
    print(f"  Number of contigs: {len(s['contigs'])}")
    if s["contig_lengths"]:
        print(f"  Max contig length: {max(s['contig_lengths'])} bp")
        print(f"  Avg contig length: "
              f"{sum(s['contig_lengths'])/len(s['contig_lengths']):.1f} bp")
    print(f"  Max coverage     : {max(s['coverage'])}")
    print(f"  Mean coverage    : "
          f"{sum(s['coverage'])/len(s['coverage']):.2f}x")
    print(f"  [Theory] α       = {s['alpha']:.4f}")
    print(f"  [Theory] Expected contigs    ≈ {s['exp_contigs']:.2f}")
    print(f"  [Theory] Expected contig len ≈ {s['exp_contig_len']:.2f} bp")


# ---------------------------------------------------------------------------
# Simulations
# ---------------------------------------------------------------------------

def shotgun_simulation_without_repeats() -> dict:
    """Simulate shotgun sequencing on a purely random genome.

    Generates a random genome of length G = 2000 bp, samples N = 1200 reads
    of length 50–100 bp, builds the overlap graph, produces alignment and
    contig plots, and returns a complete stats dict.

    Returns:
        Stats dict (see _collect_stats for keys).
    """
    genome = generate_genome(2000)
    reads  = generate_reads(genome, num_reads=1200, min_len=50, max_len=100)

    graph = build_overlap_graph(reads, min_length=5)
    draw_graph(graph, "overlap_structure_graph_without_repeats")

    visualize_alignment(genome, reads, "alignment_plot_without_repeats.png")

    s = _collect_stats(genome, reads)
    _print_stats("Genome WITHOUT repeats", s)

    fig, ax = plt.subplots(figsize=(12, 4))
    plot_contigs(ax, s["contigs"], show_labels=True, show_gaps=True)
    ax.set_xlim(0, len(genome))
    ax.set_title("Contig Assembly Results (Without Repeats)")
    ax.set_xlabel("Genome Position")
    ax.set_yticks([])
    plt.tight_layout()
    plt.savefig("contigs_only_without_repeats.png", dpi=150)
    plt.close(fig)
    print("Contig plot saved to 'contigs_only_without_repeats.png'")

    return s


def shotgun_simulation_with_repeats() -> dict:
    """Simulate shotgun sequencing on a repeat-rich genome.

    Uses the same parameters as the plain-genome simulation but generates
    a genome with LINE / SINE / LTR / tandem repeat elements. Repeat
    elements longer than the read length L violate the Poisson assumption,
    causing empirical contig counts to diverge from the N·e^{-α} prediction.

    Returns:
        Stats dict (see _collect_stats for keys).
    """
    genome = generate_genome_with_repeats(2000)
    reads  = generate_reads(genome, num_reads=1200, min_len=50, max_len=100)

    graph = build_overlap_graph(reads, min_length=5)
    draw_graph(graph, "overlap_structure_graph_with_repeats")

    visualize_alignment(genome, reads, "alignment_plot_with_repeats.png",
                        repeat_target="CAG" * 5)

    s = _collect_stats(genome, reads)
    _print_stats("Genome WITH repeats", s)

    fig, ax = plt.subplots(figsize=(12, 4))
    plot_contigs(ax, s["contigs"], show_labels=True, show_gaps=True)
    ax.set_xlim(0, len(genome))
    ax.set_title("Contig Assembly Results (With Repeats)")
    ax.set_xlabel("Genome Position")
    ax.set_yticks([])
    plt.tight_layout()
    plt.savefig("contigs_only_with_repeats.png", dpi=150)
    plt.close(fig)
    print("Contig plot saved to 'contigs_only_with_repeats.png'")

    return s


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    """Run both simulations, then produce all summary and extended plots."""
    s0 = shotgun_simulation_without_repeats()
    s1 = shotgun_simulation_with_repeats()

    generate_summary_plots(s0, s1, output_png="summary_comparison.png")

    plot_lander_waterman_sweep(
        G          = s0["G"],
        min_len    = min(s0["read_lengths"]),
        max_len    = max(s0["read_lengths"]),
        n_trials   = 3,
        output_png = "lander_waterman_sweep.png",
    )

    generate_extended_summary_plots(
        s0,
        s1,
        output_png   = "extended_summary.png",
        sweep_trials = 2,        # raise to 4–5 for smoother α-sweep dots
    )

    print("\n" + "="*54)
    print("  All output files")
    print("="*54)
    for fn in [
        "overlap_structure_graph_without_repeats.png",
        "overlap_structure_graph_with_repeats.png",
        "alignment_plot_without_repeats.png",
        "alignment_plot_with_repeats.png",
        "contigs_only_without_repeats.png",
        "contigs_only_with_repeats.png",
        "summary_comparison.png",
        "lander_waterman_sweep.png",
        "extended_summary.png",
    ]:
        print(f"  • {fn}")


if __name__ == "__main__":
    main()

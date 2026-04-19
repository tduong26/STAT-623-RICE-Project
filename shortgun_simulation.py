"""
shortgun_simulation.py
======================
Shotgun-sequencing simulation pipeline.
No pysam / BAM dependency — works on Windows, macOS, and Linux.
All alignment is performed in-memory via alignment_visualization.align_reads().
"""

import random
import os
import numpy as np
import matplotlib.pyplot as plt

from overlap import build_overlap_graph, draw_graph
from alignment_visualization import (
    align_reads,
    visualize_alignment,
    get_contigs,
    get_coverage,
    plot_contigs,
    plot_lander_waterman_sweep,
)
from theoretical_stats import compute_alpha, expected_contigs, expected_contig_length


# ---------------------------------------------------------------------------
# Genome generation
# ---------------------------------------------------------------------------

def generate_genome(length):
    """
    Generate a purely random DNA sequence (the 'ground-truth' genome).

    Args:
        length (int): Total length G of the reference genome.
    """
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def generate_genome_with_repeats(length, snp_rate=0.01, repeat_sizes=None, weights=None):
    """
    Generate a bio-realistic genome with repeat families (LINE, SINE, LTR, DNA).

    Args:
        length       (int):   Total genome length.
        snp_rate     (float): Per-base mutation rate applied to each repeat copy.
        repeat_sizes (dict):  Custom sizes per repeat family.
        weights      (list):  [(type, probability), …] for genome composition.
    """
    bases = ['A', 'C', 'G', 'T']

    def gen_random(l):
        return ''.join(random.choice(bases) for _ in range(l))

    def mutate(seq, rate):
        return ''.join(
            random.choice(bases) if random.random() < rate else b
            for b in seq
        )

    if repeat_sizes is None:
        repeat_sizes = {"LINE": 800, "SINE": 150, "LTR": 400, "DNA": 300}

    # One template per family → highly identical copies across the genome
    repeat_library = {
        name: [gen_random(size) for _ in range(1)]
        for name, size in repeat_sizes.items()
    }

    if weights is None:
        weights = [
            ("LINE",   0.20),
            ("SINE",   0.13),
            ("LTR",    0.08),
            ("DNA",    0.03),
            ("TANDEM", 0.05),
            ("UNIQUE", 0.51),
        ]

    def pick_type():
        r = random.random()
        cumulative = 0.0
        for t, w in weights:
            cumulative += w
            if r < cumulative:
                return t
        return "UNIQUE"

    genome = []
    current_len = 0

    while current_len < length:
        t = pick_type()
        if t in repeat_library:
            seq = mutate(random.choice(repeat_library[t]), snp_rate)
        elif t == "TANDEM":
            unit = gen_random(random.randint(1, 6))
            seq  = unit * random.randint(5, 20)
        else:
            seq = gen_random(100)

        genome.append(seq)
        current_len += len(seq)

    return ''.join(genome)[:length]


# ---------------------------------------------------------------------------
# Read generation
# ---------------------------------------------------------------------------

def generate_reads_theory(genome, num_reads, min_len, max_len):
    """
    Simulate shotgun sequencing: sample *num_reads* reads of uniform-random
    length in [min_len, max_len] from random positions in *genome*.

    Returns:
        list[str]: Read sequences.
    """
    reads = []
    genome_len = len(genome)
    for _ in range(num_reads):
        read_len = random.randint(min_len, max_len)
        if read_len >= genome_len:
            continue
        start = random.randint(0, genome_len - read_len)
        reads.append(genome[start:start + read_len])
    return reads


def generate_reads(genome, num_fragments, mean_len, std_len):
    """
    Simulate mechanical fragmentation with Gaussian fragment lengths.

    Returns:
        list[tuple]: (start_position, sequence) pairs.
    """
    fragments  = []
    genome_len = len(genome)
    for _ in range(num_fragments):
        frag_len = max(1, int(random.gauss(mean_len, std_len)))
        if frag_len >= genome_len:
            continue
        start = random.randint(0, genome_len - frag_len)
        fragments.append((start, genome[start:start + frag_len]))
    return fragments


# ---------------------------------------------------------------------------
# Result display helper
# ---------------------------------------------------------------------------

def display_simulation_results(label, genome, reads, output_contig_png):
    """
    Align reads, print theoretical vs. empirical stats, and save a contig plot.

    Note: the old *bam_filename* argument has been removed — alignment is now
    done entirely in memory via align_reads().

    Args:
        label            (str):  Human-readable scenario label.
        genome           (str):  Reference genome sequence.
        reads            (list): Read strings or (start, seq) tuples.
        output_contig_png (str): Path for the saved contig-only plot.
    """
    reads_data    = align_reads(genome, reads)
    contigs       = get_contigs(reads_data)
    dropped_count = len(reads) - len(reads_data)

    print(f"\n--- {label} ---")
    print(f"Number of reads sampled : {len(reads)}")
    print(f"Number of reads aligned : {len(reads_data)} "
          f"({dropped_count} dropped due to mapping ambiguity)")
    print(f"Number of contigs       : {len(contigs)}")
    if contigs:
        lengths = [e - s for s, e in contigs]
        print(f"Max contig length       : {max(lengths)} bp")
        print(f"Avg contig length       : {sum(lengths) / len(lengths):.2f} bp")

    # Theoretical expectations
    L = (max(len(r[1]) if isinstance(r, tuple) else len(r) for r in reads)
         if reads else 0)
    G     = len(genome)
    N     = len(reads)
    alpha = compute_alpha(N, L, G)
    print(f"[Theory] alpha = {alpha:.4f}")
    print(f"[Theory] Expected # contigs ≈ {expected_contigs(N, alpha):.2f}")

    # Contig-only plot
    # Use a taller figure and explicit subplots_adjust to avoid tight_layout warnings
    fig, ax = plt.subplots(figsize=(12, 5))
    plot_contigs(ax, contigs, show_labels=True, show_gaps=True)
    ax.set_xlim(0, G)
    ax.set_ylim(0.5, 2.5)          # fixed y range: enough room for labels above bars
    ax.set_title(f"Contig Assembly Results ({label})", pad=12)
    ax.set_xlabel("Genome Position (bp)")
    ax.set_yticks([])
    fig.subplots_adjust(left=0.04, right=0.98, top=0.88, bottom=0.12)
    plt.savefig(output_contig_png)
    plt.close(fig)
    print(f"Contig plot saved to {output_contig_png}")


# ---------------------------------------------------------------------------
# Scenario 1 – Human mitochondrial genome assembly
# ---------------------------------------------------------------------------

def run_alignment_comparison():
    """
    Scenario 1: Human Mitochondrial Genome Assembly (~16,569 bp).
    High unique density, small tandem repeats and D-loop variations.
    Standard Illumina read lengths (100–150 bp).
    """
    print("\n>>> Scenario 1: Human Mitochondrial Genome Assembly <<<")
    G      = 16569
    L_min  = 100
    L_max  = 150
    N      = 2500   # ~18× coverage

    # 1a. Non-repeat genome
    genome_nr = generate_genome(G)
    reads_nr  = generate_reads(genome_nr, N, L_min, L_max)
    visualize_alignment(genome_nr, reads_nr, "alignment_nr_scenario1.png")
    display_simulation_results(
        "S1: Non-Repeat Alignment (mtDNA size)",
        genome_nr, reads_nr,
        "contigs_s1_nr.png",
    )

    # 1b. Repeat genome (mtDNA: HVR + tandem repeats)
    repeat_sizes = {"HVR": 150}
    weights_r = [
        ("HVR",    0.03),   # ~3 % hypervariable region
        ("TANDEM", 0.05),   # ~5 % short tandem repeats
        ("UNIQUE", 0.92),   # ~92 % unique coding sequence
    ]
    genome_r = generate_genome_with_repeats(
        G, snp_rate=0.01,
        repeat_sizes=repeat_sizes, weights=weights_r,
    )
    reads_r = generate_reads(genome_r, N, L_min, L_max)
    visualize_alignment(genome_r, reads_r, "alignment_r_scenario1.png")
    display_simulation_results(
        "S1: Realistic mtDNA Assembly",
        genome_r, reads_r,
        "contigs_s1_r.png",
    )


# ---------------------------------------------------------------------------
# Scenario 2 – Human nuclear genome segment (fragmentation study)
# ---------------------------------------------------------------------------

def run_assembly_comparison():
    """
    Scenario 2: Contig Structure Comparison – 4,000 bp nuclear genome segment.
    Nuclear human repeats (~45-50 % repetitive): SINE/Alu, LINE, LTR, DNA.
    Standard coverage: ~30× (N=600, lengths 50–100 bp).
    """
    print("\n>>> Scenario 2: Contig Structure Comparison (Fragmentation) <<<")
    G     = 4000
    L_min = 50
    L_max = 100
    N     = 600

    # 2a. Non-repeat genome
    genome_nr = generate_genome(G)
    reads_nr  = generate_reads(genome_nr, N, L_min, L_max)
    visualize_alignment(genome_nr, reads_nr, "alignment_nr_scenario2.png")
    display_simulation_results(
        "S2: Non-Repeat Assembly",
        genome_nr, reads_nr,
        "contigs_s2_nr.png",
    )

    # 2b. Repeat genome (nuclear human proportions ~45 % repeat)
    repeat_sizes = {"LINE": 1000, "SINE": 300, "LTR": 500, "DNA": 300}
    weights_r = [
        ("LINE",   0.20),
        ("SINE",   0.13),
        ("LTR",    0.08),
        ("DNA",    0.03),
        ("UNIQUE", 0.55),
    ]
    genome_r = generate_genome_with_repeats(
        G, snp_rate=0.005,
        repeat_sizes=repeat_sizes, weights=weights_r,
    )
    reads_r = generate_reads(genome_r, N, L_min, L_max)
    visualize_alignment(genome_r, reads_r, "alignment_r_scenario2.png")
    display_simulation_results(
        "S2: Repeat Assembly (Fragmentation focused)",
        genome_r, reads_r,
        "contigs_s2_r.png",
    )


# ---------------------------------------------------------------------------
# Lander-Waterman alpha sweep
# ---------------------------------------------------------------------------

def run_lander_waterman_sweep(genome_r, output_png, label, num_trials=10):
    """
    Sweep read densities (alpha) and record mean ± 95 % CI contig counts
    for both a non-repeat baseline and the supplied repeat-rich genome.

    Args:
        genome_r   (str):  Repeat-containing genome sequence.
        output_png (str):  Output plot filename.
        label      (str):  Human-readable label for progress output.
        num_trials (int):  Trials per alpha value (more → tighter CI).
    """
    print(f"\n>>> Lander-Waterman Sweep: {label} ({num_trials} trials) <<<")
    G      = len(genome_r)
    L_min  = 50
    L_max  = 100
    L_mean = (L_min + L_max) / 2

    alpha_values = [
        0.2, 0.5, 0.8, 1.2, 1.5, 1.8,
        2.2, 2.6, 3.0, 3.25, 3.6, 4.0,
        4.4, 4.8, 5.2, 5.6, 6.0,
    ]

    genome_nr    = generate_genome(G)   # non-repeat baseline
    sweep_data_nr = []
    sweep_data_r  = []

    for alpha in alpha_values:
        N = int(alpha * G / L_mean)

        counts_nr = []
        counts_r  = []

        for _ in range(num_trials):
            # Non-repeat trial
            reads_nr   = generate_reads_theory(genome_nr, N, L_min, L_max)
            aligned_nr = align_reads(genome_nr, reads_nr)
            counts_nr.append(len(get_contigs(aligned_nr)))

            # Repeat trial
            reads_r_t  = generate_reads_theory(genome_r, N, L_min, L_max)
            aligned_r  = align_reads(genome_r, reads_r_t)
            counts_r.append(len(get_contigs(aligned_r)))

        mean_nr = float(np.mean(counts_nr))
        ci_nr   = 1.96 * float(np.std(counts_nr)) / np.sqrt(num_trials) if num_trials > 1 else 0.0

        mean_r  = float(np.mean(counts_r))
        ci_r    = 1.96 * float(np.std(counts_r))  / np.sqrt(num_trials) if num_trials > 1 else 0.0

        sweep_data_nr.append((alpha, mean_nr, ci_nr))
        sweep_data_r.append( (alpha, mean_r,  ci_r))

        print(f"  alpha={alpha:.2f}: NR={mean_nr:.2f}±{ci_nr:.2f}  "
              f"R={mean_r:.2f}±{ci_r:.2f}")

    plot_lander_waterman_sweep(G, L_mean, sweep_data_nr, sweep_data_r, output_png)


# ---------------------------------------------------------------------------
# Read-length sweep (constant coverage)
# ---------------------------------------------------------------------------

def run_read_length_sweep(genome_r, target_coverage, output_png, label, num_trials=3):
    """
    Sweep read lengths to show how longer reads bridge repeats.
    Coverage depth is held constant at *target_coverage*×.

    Args:
        genome_r        (str):  Repeat-containing genome.
        target_coverage (int):  Desired coverage depth (e.g. 30).
        output_png      (str):  Output plot filename.
        label           (str):  Human-readable label.
        num_trials      (int):  Trials per read-length point.
    """
    print(f"\n>>> Read-Length Sweep: {label} ({num_trials} trials) <<<")
    G = len(genome_r)

    length_ranges = [
        ( 50,  75),   # shorter than all repeats
        ( 75, 100),   # standard Illumina short reads
        (150, 200),   # approaching SINE size (300 bp)
        (250, 300),   # bridges SINEs and DNA transposons
        (350, 400),   # bridges LTRs (400 bp)
        (750, 800),   # bridges truncated LINEs (800 bp)
    ]

    results = []

    for L_min, L_max in length_ranges:
        L_mean = (L_min + L_max) / 2
        N      = int((target_coverage * G) / L_mean)   # constant coverage

        trial_counts = []
        for _ in range(num_trials):
            reads      = generate_reads_theory(genome_r, N, L_min, L_max)
            aligned    = align_reads(genome_r, reads)
            contigs    = get_contigs(aligned)
            trial_counts.append(len(contigs))

        avg = float(np.mean(trial_counts))
        std = float(np.std(trial_counts)) if num_trials > 1 else 0.0
        results.append((L_mean, avg, std))
        print(f"  L={L_mean:.0f} bp  N={N}: avg contigs={avg:.2f} ± {std:.2f}")

    l_means = [r[0] for r in results]
    contigs = [r[1] for r in results]
    stds    = [r[2] for r in results]

    plt.figure(figsize=(10, 6))
    plt.plot(l_means, contigs, marker='o', linewidth=2,
             color='darkorange', label='Average contigs')
    plt.fill_between(
        l_means,
        [max(0, c - s) for c, s in zip(contigs, stds)],
        [c + s         for c, s in zip(contigs, stds)],
        color='darkorange', alpha=0.2, label='\u00b11 Std Dev',
    )
    plt.axvline(x=300, color='gray', linestyle='--', label='SINE (Alu) size (300 bp)')
    plt.axvline(x=400, color='gray', linestyle=':',  label='LTR size (400 bp)')
    plt.axvline(x=800, color='gray', linestyle='-.', label='LINE size (800 bp)')

    plt.title(f"Impact of Read Length on Genome Assembly\n"
              f"(Constant {target_coverage}\u00d7 Coverage)")
    plt.xlabel("Mean Read Length (bp)")
    plt.ylabel("Number of Contigs (Fragmentation)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_png)
    plt.close()
    print(f"Read-length sweep plot saved to {output_png}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    """Run all scenarios, sweeps, alignment plots, and contig plots."""
    random.seed(42)
    np.random.seed(42)

    # =========================================================================
    # Scenario 1 - Human Mitochondrial Genome (16,569 bp)
    # =========================================================================
    G_s1 = 16569
    repeat_sizes_s1 = {"HVR": 150}
    weights_s1 = [("HVR", 0.03), ("TANDEM", 0.05), ("UNIQUE", 0.92)]

    genome_s1_nr = generate_genome(G_s1)
    genome_s1_r  = generate_genome_with_repeats(
        G_s1, snp_rate=0.01,
        repeat_sizes=repeat_sizes_s1, weights=weights_s1,
    )

    # S1 alignment visualisations
    print("\n>>> S1: Alignment visualisations <<<")
    reads_s1_nr = generate_reads(genome_s1_nr, 2500, 100, 150)
    visualize_alignment(genome_s1_nr, reads_s1_nr, "alignment_nr_scenario1.png")

    reads_s1_r = generate_reads(genome_s1_r, 2500, 100, 150)
    visualize_alignment(genome_s1_r, reads_s1_r, "alignment_r_scenario1.png")

    # S1 contig-only plots
    display_simulation_results(
        "S1: Non-Repeat (mtDNA size)", genome_s1_nr, reads_s1_nr,
        "contigs_s1_nr.png",
    )
    display_simulation_results(
        "S1: Realistic mtDNA Assembly", genome_s1_r, reads_s1_r,
        "contigs_s1_r.png",
    )

    # S1 Lander-Waterman sweep
    run_lander_waterman_sweep(
        genome_s1_r, "lander_waterman_sweep_s1.png",
        "Scenario 1 (Human mtDNA)", num_trials=5,
    )

    # =========================================================================
    # Scenario 2 - Human Nuclear Genome Segment (4,000 bp)
    # =========================================================================
    G_s2 = 4000
    repeat_sizes_s2 = {"LINE": 1000, "SINE": 300, "LTR": 500, "DNA": 300}
    weights_s2 = [
        ("LINE", 0.20), ("SINE", 0.13), ("LTR", 0.08),
        ("DNA",  0.03), ("UNIQUE", 0.55),
    ]

    genome_s2_nr = generate_genome(G_s2)
    genome_s2_r  = generate_genome_with_repeats(
        G_s2, snp_rate=0.005,
        repeat_sizes=repeat_sizes_s2, weights=weights_s2,
    )

    # S2 alignment visualisations
    print("\n>>> S2: Alignment visualisations <<<")
    reads_s2_nr = generate_reads(genome_s2_nr, 600, 50, 100)
    visualize_alignment(genome_s2_nr, reads_s2_nr, "alignment_nr_scenario2.png")

    reads_s2_r = generate_reads(genome_s2_r, 600, 50, 100)
    visualize_alignment(genome_s2_r, reads_s2_r, "alignment_r_scenario2.png")

    # S2 contig-only plots
    display_simulation_results(
        "S2: Non-Repeat Assembly", genome_s2_nr, reads_s2_nr,
        "contigs_s2_nr.png",
    )
    display_simulation_results(
        "S2: Repeat Assembly (Fragmentation)", genome_s2_r, reads_s2_r,
        "contigs_s2_r.png",
    )

    # S2 Lander-Waterman sweep
    run_lander_waterman_sweep(
        genome_s2_r, "lander_waterman_sweep_s2.png",
        "Scenario 2 (Human Nuclear)", num_trials=5,
    )

    # S2 Read-length sweep
    run_read_length_sweep(
        genome_s2_r, target_coverage=30,
        output_png="read_length_sweep_s2.png",
        label="Scenario 2 Length Impact", num_trials=5,
    )


if __name__ == "__main__":
    main()

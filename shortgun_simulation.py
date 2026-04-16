"""
shortgun_simulation.py  (Windows-compatible version)
=====================================================
All pysam calls have been replaced with the pure-Python in-memory alignment
store defined in alignment_visualization.py.  No BAM files are written to
disk; everything is kept in memory, which also makes the code faster on
Windows where file I/O is slower than on Unix systems.
"""

import random
import matplotlib.pyplot as plt

from overlap import build_overlap_graph, draw_graph
from alignment_visualization import (
    visualize_alignment,
    get_contigs,
    get_coverage,
    plot_contigs,
    create_alignment_bam,
    plot_lander_waterman_sweep,
    plot_effective_alpha_comparison,
    effective_alpha,
    bootstrap_ci,
    _open_bam,      # in-memory BAM accessor
    _remove_bam,    # in-memory BAM cleanup
)
from theoretical_stats import (
    compute_alpha,
    expected_contigs,
    expected_contig_length,
)


# ---------------------------------------------------------------------------
# Genome generation
# ---------------------------------------------------------------------------

def generate_genome(length: int) -> str:
    """
    Generates a purely random DNA sequence (the 'ground truth' genome).

    Args:
        length (int): Total length G of the reference genome.
    """
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def generate_genome_with_repeats(length: int, snp_rate: float = 0.01,
                                  repeat_sizes: dict = None,
                                  weights: list = None) -> str:
    """
    Generates a bio-realistic genome with repeat families (LINE, SINE, LTR, DNA).

    Args:
        length (int):           Total length of the genome.
        snp_rate (float):       Mutation rate applied to each repeat copy.
        repeat_sizes (dict):    Custom lengths for each repeat family.
        weights (list):         Probability weights for each repeat type.
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

    # One template per family – makes copies highly identical across the genome
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
        return "UNIQUE"   # fallback

    genome = []
    current_len = 0

    while current_len < length:
        t = pick_type()

        if t in repeat_library:
            family = random.choice(repeat_library[t])
            seq = mutate(family, snp_rate)
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

def generate_reads(genome: str, num_reads: int,
                   min_len: int, max_len: int) -> list[str]:
    """
    Simulates shotgun sequencing by generating random reads from a genome.
    """
    reads = []
    genome_len = len(genome)
    for _ in range(num_reads):
        read_len = random.randint(min_len, max_len)
        start    = random.randint(0, genome_len - read_len)
        reads.append(genome[start: start + read_len])
    return reads


# ---------------------------------------------------------------------------
# Shared result display helper
# ---------------------------------------------------------------------------

def display_simulation_results(label: str, genome: str, reads: list[str],
                                bam_key: str, output_contig_png: str):
    """
    Align reads, extract stats, print theoretical vs empirical, save contig plot.

    Also computes and prints effective alpha (MAPQ-60 only) alongside nominal alpha.
    *bam_key* is just an identifier string; no file is written to disk.
    """
    create_alignment_bam(genome, reads, bam_key)

    bam = _open_bam(bam_key)
    reads_data = [
        {
            'start': r.reference_start,
            'end':   r.reference_end,
            'mapq':  r.mapping_quality,
        }
        for r in bam.fetch()
    ]
    _remove_bam(bam_key)

    contigs  = get_contigs(reads_data)
    coverage = get_coverage(len(genome), reads_data)

    dropped_count = len(reads) - len(reads_data)
    L_mean = sum(len(r) for r in reads) / len(reads) if reads else 0

    print(f"\n--- {label} ---")
    print(f"Number of reads sampled: {len(reads)}")
    print(
        f"Number of reads aligned: {len(reads_data)} "
        f"({dropped_count} dropped due to mapping ambiguity)"
    )
    print(f"Number of contigs assembled: {len(contigs)}")
    if contigs:
        lengths = [e - s for s, e in contigs]
        print(f"Max contig length: {max(lengths)} bp")
        print(f"Average contig length: {sum(lengths) / len(lengths):.2f} bp")

    # Theoretical expectations
    L = max(len(r) for r in reads) if reads else 0
    G = len(genome)
    N = len(reads)
    alpha_t    = compute_alpha(N, L, G)
    exp_ctgs   = expected_contigs(N, alpha_t)
    print(f"[Theory] alpha (read density) = {alpha_t:.4f}")
    print(f"[Theory] Expected # contigs ~ {exp_ctgs:.2f}")

    # Effective alpha (MAPQ-60 only)
    a_nom, a_eff, loss = effective_alpha(reads_data, G, L_mean)
    print(f"[Coverage] Nominal  alpha = {a_nom:.4f}")
    print(f"[Coverage] Effective alpha (MAPQ-60) = {a_eff:.4f}")
    print(f"[Coverage] MAPQ loss = {loss * 100:.1f}% of coverage is from multi-mapped reads")

    # Save contig-only plot
    fig, ax = plt.subplots(figsize=(12, 4))
    plot_contigs(ax, contigs, show_labels=True, show_gaps=True)
    ax.set_xlim(0, len(genome))
    ax.set_title(f"Contig Assembly Results ({label})")
    ax.set_xlabel("Genome Position (bp)")
    ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(output_contig_png)
    plt.close(fig)
    print(f"Contig-only plot saved to {output_contig_png}")

    # Return reads_data so callers can collect effective_alpha results
    return reads_data


# ---------------------------------------------------------------------------
# Scenario runners
# ---------------------------------------------------------------------------

def run_alignment_comparison():
    """
    Scenario 1: Highlighting Mapping Quality (MAPQ Colors)
    - High coverage (N=600).
    - Realistic repeats, short elements (50-100 bp reads).
    Also produces an effective-alpha comparison plot.
    """
    print("\n>>> Scenario 1: Alignment Visualization (MAPQ Colors) <<<")
    G        = 4000
    L_min, L_max = 50, 100
    L_mean   = (L_min + L_max) / 2
    N        = 600

    eff_alpha_results = []

    # 1. Non-repeat genome
    genome_nr = generate_genome(G)
    reads_nr  = generate_reads(genome_nr, N, L_min, L_max)
    visualize_alignment(genome_nr, reads_nr, "alignment_nr_scenario1.png")
    rd_nr = display_simulation_results(
        "S1: Non-Repeat Alignment",
        genome_nr, reads_nr,
        "bam_s1_nr", "contigs_s1_nr.png"
    )
    a_nom, a_eff, loss = effective_alpha(rd_nr, G, L_mean)
    eff_alpha_results.append({
        'label': 'S1 Non-Repeat',
        'alpha_nominal': a_nom, 'alpha_effective': a_eff, 'mapq_loss': loss
    })

    # 2. Repeat genome (realistic repeats)
    repeat_sizes = {"LINE": 1000, "SINE": 300, "LTR": 500, "DNA": 300}
    genome_r  = generate_genome_with_repeats(G, snp_rate=0.005,
                                              repeat_sizes=repeat_sizes)
    reads_r   = generate_reads(genome_r, N, L_min, L_max)
    visualize_alignment(genome_r, reads_r, "alignment_r_scenario1.png")
    rd_r = display_simulation_results(
        "S1: Repeat Alignment (MAPQ focused)",
        genome_r, reads_r,
        "bam_s1_r", "contigs_s1_r.png"
    )
    a_nom, a_eff, loss = effective_alpha(rd_r, G, L_mean)
    eff_alpha_results.append({
        'label': 'S1 Repeat',
        'alpha_nominal': a_nom, 'alpha_effective': a_eff, 'mapq_loss': loss
    })

    return eff_alpha_results


def run_assembly_comparison():
    """
    Scenario 2: Highlighting Contig Fragmentation (Gaps)
    - Intermediate coverage (N=150).
    - Long repeat elements (400-800 bp).
    Also returns effective-alpha results for the combined comparison plot.
    """
    print("\n>>> Scenario 2: Contig Structure Comparison (Fragmentation) <<<")
    G        = 2000
    L_min, L_max = 50, 100
    L_mean   = (L_min + L_max) / 2
    N        = 150

    eff_alpha_results = []

    # 1. Non-repeat genome
    genome_nr = generate_genome(G)
    reads_nr  = generate_reads(genome_nr, N, L_min, L_max)
    visualize_alignment(genome_nr, reads_nr, "alignment_nr_scenario2.png")
    rd_nr = display_simulation_results(
        "S2: Non-Repeat Assembly",
        genome_nr, reads_nr,
        "bam_s2_nr", "contigs_s2_nr.png"
    )
    a_nom, a_eff, loss = effective_alpha(rd_nr, G, L_mean)
    eff_alpha_results.append({
        'label': 'S2 Non-Repeat',
        'alpha_nominal': a_nom, 'alpha_effective': a_eff, 'mapq_loss': loss
    })

    # 2. Repeat genome (large elements + very high frequency)
    weights_r = [
        ("LINE",   0.40),
        ("SINE",   0.30),
        ("LTR",    0.10),
        ("DNA",    0.10),
        ("TANDEM", 0.00),
        ("UNIQUE", 0.10),
    ]
    genome_r = generate_genome_with_repeats(G, snp_rate=0.005,
                                             repeat_sizes=None,
                                             weights=weights_r)
    reads_r  = generate_reads(genome_r, N, L_min, L_max)
    visualize_alignment(genome_r, reads_r, "alignment_r_scenario2.png")
    rd_r = display_simulation_results(
        "S2: Repeat Assembly (Fragmentation focused)",
        genome_r, reads_r,
        "bam_s2_r", "contigs_s2_r.png"
    )
    a_nom, a_eff, loss = effective_alpha(rd_r, G, L_mean)
    eff_alpha_results.append({
        'label': 'S2 Repeat',
        'alpha_nominal': a_nom, 'alpha_effective': a_eff, 'mapq_loss': loss
    })

    return eff_alpha_results


def run_lander_waterman_sweep(genome_r: str, output_png: str,
                               label: str, num_trials: int = 10):
    """
    Sweeps read densities (alpha) and plots theoretical vs empirical contig counts
    with bootstrap 95% CI bands.

    Averages over *num_trials* independent trials for both a non-repeat and the
    provided repeat-rich genome.
    """
    print(f"\n>>> Running Lander-Waterman Sweep: {label} ({num_trials} trials) <<<")
    G      = len(genome_r)
    L_min, L_max = 50, 100
    L_mean = (L_min + L_max) / 2

    alpha_values = [
        0.2, 0.5, 0.8, 1.2, 1.5, 1.8, 2.2, 2.6, 3.0,
        3.25, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0
    ]

    sweep_data_nr = []
    sweep_data_r  = []
    ci_nr_list    = []
    ci_r_list     = []
    genome_nr     = generate_genome(G)   # non-repeat baseline

    for alpha in alpha_values:
        N = int(alpha * G / L_mean)

        trial_counts_nr = []
        trial_counts_r  = []

        for trial in range(num_trials):
            # ---- Non-repeat trial ----
            reads_nr = generate_reads(genome_nr, N, L_min, L_max)
            bam_key_nr = f"sweep_nr_{alpha}_{trial}"
            create_alignment_bam(genome_nr, reads_nr, bam_key_nr)
            bam_nr = _open_bam(bam_key_nr)
            reads_data_nr = [
                {'start': r.reference_start, 'end': r.reference_end}
                for r in bam_nr.fetch()
            ]
            _remove_bam(bam_key_nr)
            trial_counts_nr.append(len(get_contigs(reads_data_nr)))

            # ---- Repeat trial ----
            reads_r = generate_reads(genome_r, N, L_min, L_max)
            bam_key_r = f"sweep_r_{alpha}_{trial}"
            create_alignment_bam(genome_r, reads_r, bam_key_r)
            bam_r = _open_bam(bam_key_r)
            reads_data_r = [
                {'start': r.reference_start, 'end': r.reference_end}
                for r in bam_r.fetch()
            ]
            _remove_bam(bam_key_r)
            trial_counts_r.append(len(get_contigs(reads_data_r)))

        mean_nr = sum(trial_counts_nr) / num_trials
        mean_r  = sum(trial_counts_r)  / num_trials
        sweep_data_nr.append((alpha, mean_nr))
        sweep_data_r.append((alpha, mean_r))

        # Bootstrap 95% CI (only meaningful with enough trials)
        if num_trials >= 3:
            ci_nr_list.append(bootstrap_ci(trial_counts_nr))
            ci_r_list.append(bootstrap_ci(trial_counts_r))

        print(f"  Alpha {alpha:.2f}: NR_mean={mean_nr:.2f}, R_mean={mean_r:.2f}")

    ci_nr = ci_nr_list if ci_nr_list else None
    ci_r  = ci_r_list  if ci_r_list  else None
    plot_lander_waterman_sweep(G, L_mean, sweep_data_nr, sweep_data_r,
                               output_png, ci_nr=ci_nr, ci_r=ci_r)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    """Executes all comparison scenarios and the Lander-Waterman sweeps."""
    eff_s1 = run_alignment_comparison()
    eff_s2 = run_assembly_comparison()

    # Combined effective-alpha comparison plot (all 4 conditions side-by-side)
    all_eff = eff_s1 + eff_s2
    plot_effective_alpha_comparison(all_eff, "effective_alpha_comparison.png")
    print("\nEffective alpha comparison plot saved to effective_alpha_comparison.png")

    # Sweep for Scenario 1 (realistic repeats, G=4000)
    G_s1         = 4000
    repeat_s1    = {"LINE": 1000, "SINE": 300, "LTR": 500, "DNA": 300}
    genome_s1    = generate_genome_with_repeats(G_s1, snp_rate=0.005,
                                                 repeat_sizes=repeat_s1)
    run_lander_waterman_sweep(
        genome_s1, "lander_waterman_sweep_s1.png",
        "Scenario 1 (Realistic Repeats)", num_trials=3
    )

    # Sweep for Scenario 2 (long slightly-mutated repeats, G=2000)
    G_s2     = 2000
    weights_s2 = [
        ("LINE",   0.40), ("SINE",   0.30), ("LTR",    0.10),
        ("DNA",    0.10), ("TANDEM", 0.00), ("UNIQUE", 0.10),
    ]
    genome_s2 = generate_genome_with_repeats(G_s2, snp_rate=0.005,
                                              repeat_sizes=None,
                                              weights=weights_s2)
    run_lander_waterman_sweep(
        genome_s2, "lander_waterman_sweep_s2.png",
        "Scenario 2 (Long Repeats)", num_trials=5
    )


if __name__ == "__main__":
    main()

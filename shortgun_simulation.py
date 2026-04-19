import random
import matplotlib.pyplot as plt
from overlap import build_overlap_graph, draw_graph
from alignment_visualization import (
    visualize_alignment, get_contigs, get_coverage,
    plot_contigs, align_reads, plot_lander_waterman_sweep
)
from theoretical_stats import compute_alpha, expected_contigs, expected_contig_length


# ─────────────────────────────────────────────────────────────────────────────
# Genome generators
# ─────────────────────────────────────────────────────────────────────────────

def generate_genome(length):
    """
    Generates a purely random DNA sequence (the 'ground truth' genome).

    Args:
        length (int): Total length G of the reference genome.
    """
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def generate_genome_with_repeats(length, snp_rate=0.01, repeat_sizes=None, weights=None):
    """
    Generates a bio-realistic genome with repeat families (LINE, SINE, LTR, DNA).

    Args:
        length (int): Total genome length.
        snp_rate (float): Mutation rate applied to each repeat copy.
        repeat_sizes (dict, optional): Custom lengths for each repeat family.
        weights (list, optional): Probability weights for each repeat type.
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

    # One template per family makes repeats highly identical across the genome
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
        cumulative = 0
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


# ─────────────────────────────────────────────────────────────────────────────
# Read generators
# ─────────────────────────────────────────────────────────────────────────────

def generate_reads_theory(genome, num_reads, min_len, max_len):
    """
    Simulates shotgun sequencing: uniform random reads, uniform start positions.
    """
    reads = []
    genome_len = len(genome)
    for _ in range(num_reads):
        read_len = random.randint(min_len, max_len)
        start    = random.randint(0, genome_len - read_len)
        reads.append(genome[start:start + read_len])
    return reads


def generate_reads(genome, num_fragments, mean_len, std_len):
    """
    Gaussian-length fragmentation model.
    Returns (start, sequence) tuples.
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


# ─────────────────────────────────────────────────────────────────────────────
# Display helpers
# ─────────────────────────────────────────────────────────────────────────────

def display_simulation_results(label, genome, reads, output_contig_png):
    """
    Aligns reads in memory (no BAM files), prints theoretical vs empirical
    statistics, and saves a contig-only plot.

    Args:
        label (str): Description label for this run.
        genome (str): Reference genome string.
        reads (list): Reads (strings or (start, seq) tuples).
        output_contig_png (str): Output filename for the contig plot.
    """
    reads_data = align_reads(genome, reads)
    contigs    = get_contigs(reads_data)
    coverage   = get_coverage(len(genome), reads_data)

    dropped_count = len(reads) - len(reads_data)

    print(f"\n--- {label} ---")
    print(f"Number of reads sampled:  {len(reads)}")
    print(f"Number of reads aligned:  {len(reads_data)} "
          f"({dropped_count} dropped due to mapping ambiguity)")
    print(f"Number of contigs assembled: {len(contigs)}")
    if contigs:
        lengths = [e - s for s, e in contigs]
        print(f"Max contig length:     {max(lengths)} bp")
        print(f"Average contig length: {sum(lengths)/len(lengths):.2f} bp")

    if reads:
        L = max(len(r[1]) if isinstance(r, tuple) else len(r) for r in reads)
    else:
        L = 0
    G     = len(genome)
    N     = len(reads)
    alpha = compute_alpha(N, L, G)
    print(f"[Theory] alpha (read density) = {alpha:.4f}")
    print(f"[Theory] Expected # contigs  ~ {expected_contigs(N, alpha):.2f}")

    fig, ax = plt.subplots(figsize=(12, 4))
    plot_contigs(ax, contigs, show_labels=True, show_gaps=True)
    ax.set_xlim(0, len(genome))
    ax.set_title(f"Contig Assembly Results ({label})")
    ax.set_xlabel("Genome Position (bp)")
    ax.set_yticks([])
    plt.tight_layout()
    plt.savefig(output_contig_png, dpi=150)
    print(f"Contig-only plot saved to {output_contig_png}")


# ─────────────────────────────────────────────────────────────────────────────
# Scenario runners
# ─────────────────────────────────────────────────────────────────────────────

def run_alignment_comparison():
    """
    Scenario 1: Human Mitochondrial Genome Assembly
      - Realistic length (~16,569 bp).
      - High unique density, small tandem repeats & D-loop variations.
      - Standard Illumina read lengths (100-150 bp).
    """
    print("\n>>> Scenario 1: Human Mitochondrial Genome Assembly <<<")
    G = 16569
    L_min, L_max = 100, 150
    N = 2500   # ~18x coverage

    # 1a. Non-repeat genome
    genome_nr = generate_genome(G)
    reads_nr  = generate_reads(genome_nr, N, L_min, L_max)
    visualize_alignment(genome_nr, reads_nr, "alignment_nr_scenario1.png")
    display_simulation_results("S1: Non-Repeat Alignment (mtDNA size)",
                               genome_nr, reads_nr, "contigs_s1_nr.png")

    # 1b. Repeat genome (mtDNA realistic repeats)
    repeat_sizes = {"HVR": 150}
    weights_r = [
        ("HVR",    0.03),
        ("TANDEM", 0.05),
        ("UNIQUE", 0.92),
    ]
    genome_r = generate_genome_with_repeats(
        G, snp_rate=0.01, repeat_sizes=repeat_sizes, weights=weights_r)
    reads_r  = generate_reads(genome_r, N, L_min, L_max)
    visualize_alignment(genome_r, reads_r, "alignment_r_scenario1.png")
    display_simulation_results("S1: Realistic mtDNA Assembly",
                               genome_r, reads_r, "contigs_s1_r.png")


def run_assembly_comparison():
    """
    Scenario 2: Contig Structure Comparison (Fragmentation)
      - Human genome segment (4 000 bp).
      - Nuclear human repeats (~45-50 % repetitive).
      - Standard coverage: ~30x.
    """
    print("\n>>> Scenario 2: Contig Structure Comparison (Fragmentation) <<<")
    G = 4000
    L_min, L_max = 50, 100
    N = 600

    # 2a. Non-repeat genome
    genome_nr = generate_genome(G)
    reads_nr  = generate_reads(genome_nr, N, L_min, L_max)
    visualize_alignment(genome_nr, reads_nr, "alignment_nr_scenario2.png")
    display_simulation_results("S2: Non-Repeat Assembly",
                               genome_nr, reads_nr, "contigs_s2_nr.png")

    # 2b. Repeat genome (Nuclear human genome proportions ~45 % repeat)
    repeat_sizes = {"LINE": 1000, "SINE": 300, "LTR": 500, "DNA": 300}
    weights_r = [
        ("LINE",   0.20),
        ("SINE",   0.13),
        ("LTR",    0.08),
        ("DNA",    0.03),
        ("UNIQUE", 0.55),
    ]
    genome_r = generate_genome_with_repeats(
        G, snp_rate=0.005, repeat_sizes=repeat_sizes, weights=weights_r)
    reads_r  = generate_reads(genome_r, N, L_min, L_max)
    visualize_alignment(genome_r, reads_r, "alignment_r_scenario2.png")
    display_simulation_results("S2: Repeat Assembly (Fragmentation focused)",
                               genome_r, reads_r, "contigs_s2_r.png")


def run_lander_waterman_sweep(genome_r, output_png, label, num_trials=10):
    """
    Sweeps read density (alpha) and collects per-trial contig counts for both
    a non-repeat baseline and the provided repeat-rich genome.

    Raw trial lists are passed to plot_lander_waterman_sweep so it can
    compute 95 % bootstrap CIs automatically.

    Args:
        genome_r (str): Repeat-rich genome string.
        output_png (str): Output filename for the sweep plot.
        label (str): Human-readable label for console output.
        num_trials (int): Number of independent trials per alpha value.
    """
    print(f"\n>>> Running Lander-Waterman Sweep: {label} ({num_trials} trials) <<<")
    G      = len(genome_r)
    L_min, L_max = 50, 100
    L_mean = (L_min + L_max) / 2
    alpha_values = [
        0.2, 0.5, 0.8, 1.2, 1.5, 1.8,
        2.2, 2.6, 3.0, 3.25, 3.6, 4.0,
        4.4, 4.8, 5.2, 5.6, 6.0,
    ]

    sweep_data_nr = []
    sweep_data_r  = []
    genome_nr     = generate_genome(G)   # non-repeat baseline

    for alpha in alpha_values:
        N = int(alpha * G / L_mean)

        trial_counts_nr = []
        trial_counts_r  = []

        for _ in range(num_trials):
            # Non-repeat trial — align in memory
            reads_nr     = generate_reads_theory(genome_nr, N, L_min, L_max)
            reads_data_nr = align_reads(genome_nr, reads_nr)
            trial_counts_nr.append(len(get_contigs(reads_data_nr)))

            # Repeat trial — align in memory
            reads_r      = generate_reads_theory(genome_r, N, L_min, L_max)
            reads_data_r  = align_reads(genome_r, reads_r)
            trial_counts_r.append(len(get_contigs(reads_data_r)))

        mean_nr = sum(trial_counts_nr) / num_trials
        mean_r  = sum(trial_counts_r)  / num_trials

        # Store raw lists — the plot function uses them to compute bootstrap CIs
        sweep_data_nr.append((alpha, trial_counts_nr))
        sweep_data_r.append((alpha, trial_counts_r))

        print(f"  alpha={alpha:.2f}  NR_mean={mean_nr:.2f}  R_mean={mean_r:.2f}")

    plot_lander_waterman_sweep(G, L_mean, sweep_data_nr, sweep_data_r, output_png)


def run_read_length_sweep(genome_r, target_coverage, output_png, label, num_trials=3):
    """
    Sweeps read lengths to show how longer reads bridge repeats and reduce
    fragmentation, while maintaining constant coverage depth.

    Args:
        genome_r (str): Repeat-rich genome string.
        target_coverage (int): Target mean coverage depth (e.g. 30).
        output_png (str): Output filename.
        label (str): Human-readable label for console output.
        num_trials (int): Number of independent trials per read-length bin.
    """
    print(f"\n>>> Running Read Length Sweep: {label} ({num_trials} trials) <<<")
    G = len(genome_r)

    length_ranges = [
        (50,  75),    # shorter than all repeats
        (75,  100),   # standard Illumina
        (150, 200),   # approaching SINE size (300)
        (250, 300),   # bridges SINEs (300) and DNA (200)
        (350, 400),   # bridges LTRs (400)
        (750, 800),   # bridges truncated LINEs (800)
    ]

    results = []

    for L_min, L_max in length_ranges:
        L_mean = (L_min + L_max) / 2
        N      = int((target_coverage * G) / L_mean)

        trial_counts = []
        for _ in range(num_trials):
            reads      = generate_reads_theory(genome_r, N, L_min, L_max)
            reads_data = align_reads(genome_r, reads)
            trial_counts.append(len(get_contigs(reads_data)))

        avg_contigs = sum(trial_counts) / num_trials
        results.append((L_mean, avg_contigs))
        print(f"  Read length {L_mean:.0f} bp  (N={N}): avg contigs = {avg_contigs:.2f}")

    l_means = [r[0] for r in results]
    contigs  = [r[1] for r in results]

    plt.figure(figsize=(10, 6))
    plt.plot(l_means, contigs, marker='o', linewidth=2, color='darkorange')
    plt.axvline(x=300, color='gray', linestyle='--', label='SINE (Alu) size (300 bp)')
    plt.axvline(x=400, color='gray', linestyle=':',  label='LTR size (400 bp)')
    plt.axvline(x=800, color='gray', linestyle='-.', label='LINE size (800 bp)')
    plt.title(f"Impact of Read Length on Genome Assembly\n(Constant {target_coverage}x Coverage)")
    plt.xlabel("Mean Read Length (bp)")
    plt.ylabel("Number of Contigs (Fragmentation)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_png, dpi=150)
    print(f"Read length sweep plot saved to {output_png}")


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    """
    Runs ALL scenarios and sweeps, producing these output PNGs:

    Alignment & contig plots (Scenario 1 — mtDNA):
        alignment_nr_scenario1.png
        alignment_r_scenario1.png
        contigs_s1_nr.png
        contigs_s1_r.png

    Alignment & contig plots (Scenario 2 — Nuclear):
        alignment_nr_scenario2.png
        alignment_r_scenario2.png
        contigs_s2_nr.png
        contigs_s2_r.png

    Lander-Waterman sweep plots (with 95% bootstrap CI):
        lander_waterman_sweep_s1.png
        lander_waterman_sweep_s2.png

    Read-length sweep plot:
        read_length_sweep_s2.png
    """

    # ── Scenario 1: Human Mitochondrial Genome ────────────────────────────────
    # Produces alignment + contig PNGs for mtDNA (16,569 bp)
    print("\n========== SCENARIO 1: Alignment & Contigs ==========")
    run_alignment_comparison()

    # Lander-Waterman sweep for Scenario 1
    print("\n========== SCENARIO 1: Lander-Waterman Sweep ==========")
    G_s1 = 16569
    repeat_sizes_s1 = {"HVR": 150}
    weights_s1 = [("HVR", 0.03), ("TANDEM", 0.05), ("UNIQUE", 0.92)]
    genome_s1 = generate_genome_with_repeats(
        G_s1, snp_rate=0.01, repeat_sizes=repeat_sizes_s1, weights=weights_s1)
    run_lander_waterman_sweep(
        genome_s1, "lander_waterman_sweep_s1.png",
        "Scenario 1 (Human mtDNA)", num_trials=10)

    # ── Scenario 2: Human Nuclear Genome Segment ──────────────────────────────
    # Produces alignment + contig PNGs for a 4,000 bp nuclear segment
    print("\n========== SCENARIO 2: Alignment & Contigs ==========")
    run_assembly_comparison()

    # Lander-Waterman sweep for Scenario 2
    print("\n========== SCENARIO 2: Lander-Waterman Sweep ==========")
    G_s2 = 4000
    repeat_sizes_s2 = {"LINE": 1000, "SINE": 300, "LTR": 500, "DNA": 300}
    weights_s2 = [
        ("LINE", 0.20), ("SINE", 0.13), ("LTR", 0.08),
        ("DNA",  0.03), ("UNIQUE", 0.55),
    ]
    genome_s2 = generate_genome_with_repeats(
        G_s2, snp_rate=0.005, repeat_sizes=repeat_sizes_s2, weights=weights_s2)
    run_lander_waterman_sweep(
        genome_s2, "lander_waterman_sweep_s2.png",
        "Scenario 2 (Human Nuclear)", num_trials=5)

    # Read-length sweep (Scenario 2 extension)
    print("\n========== SCENARIO 2: Read Length Sweep ==========")
    run_read_length_sweep(
        genome_s2, target_coverage=30,
        output_png="read_length_sweep_s2.png",
        label="Scenario 2 Length Impact", num_trials=5)


if __name__ == "__main__":
    main()

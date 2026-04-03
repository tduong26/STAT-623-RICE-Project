"""
genome_utils.py
---------------
Genome generation and shotgun read sampling utilities.

Two genome models are provided:
  - generate_genome()              : purely random (uniform ACGT)
  - generate_genome_with_repeats() : realistic repeat structure
      simulating LINE / SINE / LTR / DNA transposon families and
      tandem repeats, mimicking eukaryotic genome composition.
"""

import random


def generate_genome(length: int) -> str:
    """Generate a purely random DNA sequence (the 'ground truth' genome).

    Under this model read start positions are truly uniform U(0, G), so the
    Lander-Waterman Poisson assumptions hold exactly.

    Args:
        length: Total length G of the reference genome in bp.
    Returns:
        Random DNA string of the requested length.
    """
    bases = ["A", "C", "G", "T"]
    return "".join(random.choice(bases) for _ in range(length))


def generate_genome_with_repeats(length: int) -> str:
    """Generate a DNA sequence with realistic repeat structure.

    Simulates LINE, SINE, LTR, and DNA transposon families plus tandem
    repeats at proportions loosely inspired by mammalian genome composition.
    Repeats are created once as 'family' sequences, then inserted as mutated
    copies throughout the genome.

    When repeat elements are longer than the read length L, the uniform-start
    assumption breaks down: multiple reads map to identical positions, inflating
    the empirical contig count above the Lander-Waterman prediction.

    Args:
        length: Total genome length in bp.
    Returns:
        DNA string of exactly `length` bp.
    """
    bases = ["A", "C", "G", "T"]

    def gen_random(l: int) -> str:
        return "".join(random.choice(bases) for _ in range(l))

    def mutate(seq: str, snp_rate: float = 0.1) -> str:
        """Apply random SNPs at the given per-base rate."""
        return "".join(
            random.choice(bases) if random.random() < snp_rate else b
            for b in seq
        )

    # ── Repeat-element families ─────────────────────────────────────────────
    repeat_library = {
        "LINE": [gen_random(6000) for _ in range(5)],   # ~20% of genome
        "SINE": [gen_random(300)  for _ in range(10)],  # ~13%
        "LTR":  [gen_random(2000) for _ in range(5)],   # ~8%
        "DNA":  [gen_random(2500) for _ in range(3)],   # ~3%
    }

    weights = [
        ("LINE",   0.20),
        ("SINE",   0.13),
        ("LTR",    0.08),
        ("DNA",    0.03),
        ("TANDEM", 0.05),
        ("UNIQUE", 0.51),
    ]

    def pick_type() -> str:
        r, cumulative = random.random(), 0.0
        for t, w in weights:
            cumulative += w
            if r < cumulative:
                return t
        return "UNIQUE"

    # ── Assemble genome by concatenating random segments ────────────────────
    genome, current_len = [], 0
    while current_len < length:
        t = pick_type()
        if t in repeat_library:
            seq = mutate(random.choice(repeat_library[t]))
        elif t == "TANDEM":
            unit = gen_random(random.randint(1, 6))
            seq  = unit * random.randint(5, 20)
        else:
            seq = gen_random(100)
        genome.append(seq)
        current_len += len(seq)

    return "".join(genome)[:length]


def generate_reads(genome: str, num_reads: int,
                   min_len: int, max_len: int) -> list:
    """Simulate shotgun sequencing by sampling random substrings.

    Each read is a substring of the genome with a uniformly random start
    position and a uniformly random length in [min_len, max_len].

    Args:
        genome:    Full genome sequence to sample from.
        num_reads: Number of reads N to generate.
        min_len:   Minimum read length (bp).
        max_len:   Maximum read length (bp); also used as L in theory formulae.
    Returns:
        List of read sequences (strings).
    """
    reads, genome_len = [], len(genome)
    for _ in range(num_reads):
        read_len = random.randint(min_len, max_len)
        start    = random.randint(0, genome_len - read_len)
        reads.append(genome[start : start + read_len])
    return reads

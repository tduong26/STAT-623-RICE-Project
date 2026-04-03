"""
theoretical_stats.py
--------------------
Lander-Waterman theoretical predictions for shotgun sequencing assembly.

All three formulas assume:
  - Read start positions are uniformly distributed: U(0, G)
  - Read count in any interval follows a Poisson process
  - α = N·L / G  (read density / coverage parameter)
"""

import math


def compute_alpha(N: int, L: int, G: int) -> float:
    """Compute the read density parameter α = N * L / G.

    Args:
        N: Number of reads.
        L: Read length (maximum read length used as an approximation).
        G: Genome length.
    Returns:
        α as a float.
    """
    if G == 0:
        return 0.0
    return N * L / G


def expected_contigs(N: int, alpha: float) -> float:
    """Expected number of contigs ≈ N · e^{-α}.

    Each read has probability e^{-α} of being the leftmost read of a new
    contig (i.e. no other read ends just before it).

    Args:
        N:     Number of reads.
        alpha: Read density parameter α = N·L/G.
    Returns:
        Expected contig count as a float.
    """
    return N * math.exp(-alpha)


def expected_contig_length(G: int, N: int, alpha: float) -> float:
    """Expected contig length ≈ G · (1 − e^{-α}) / (N · e^{-α}).

    Derived as:
        E[length] = total_covered_bases / expected_number_of_contigs
                  = G·(1 − e^{-α})  /  N·e^{-α}

    Args:
        G:     Genome length in bp.
        N:     Number of reads.
        alpha: Read density parameter α = N·L/G.
    Returns:
        Expected contig length in bp.
    """
    if N == 0:
        return 0.0
    return G * (1 - math.exp(-alpha)) / (N * math.exp(-alpha))

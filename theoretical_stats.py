import math

def compute_alpha(N: int, L: int, G: int) -> float:
    """Compute the read density parameter α = N * L / G.

    Args:
        N: Number of reads.
        L: Read length (we use the maximum read length as an approximation).
        G: Genome length.
    Returns:
        α as a float.
    """
    if G == 0:
        return 0.0
    return N * L / G

def expected_contigs(N: int, alpha: float) -> float:
    """Expected number of contigs ≈ N * e^{-α}.

    Args:
        N: Number of reads.
        alpha: Read density parameter.
    Returns:
        Expected contig count.
    """
    return N * math.exp(-alpha)

def expected_contig_length(G: int, N: int, alpha: float) -> float:
    """Expected contig length ≈ G * (1 - e^{-α}) / (N * e^{-α}).

    Args:
        G: Genome length.
        N: Number of reads.
        alpha: Read density parameter.
    Returns:
        Expected contig length.
    """
    if N == 0:
        return 0.0
    return G * (1 - math.exp(-alpha)) / (N * math.exp(-alpha))

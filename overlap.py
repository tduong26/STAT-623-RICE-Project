"""
overlap.py
----------
Overlap graph construction and visualisation for shotgun read assembly.

An overlap graph has one node per read; a directed edge (a → b) exists when
the suffix of read a matches a prefix of read b by at least `min_length` bp.
"""

import matplotlib.pyplot as plt
import networkx as nx


def overlap(a: str, b: str, min_length: int = 3) -> int:
    """Check if the end of string 'a' overlaps with the beginning of string 'b'.

    Args:
        a:          First DNA sequence.
        b:          Second DNA sequence.
        min_length: Minimum required overlap length.
    Returns:
        Length of the overlap (> 0), or 0 if none found.
    """
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1


def build_overlap_graph(reads: list, min_length: int = 3) -> dict:
    """Construct an adjacency-list overlap graph of the given reads.

    Args:
        reads:      List of DNA read sequences.
        min_length: Minimum overlap length required to add a directed edge.
    Returns:
        dict mapping each read → list of (neighbour_read, overlap_length).
    """
    graph = {}
    for a in reads:
        graph[a] = []
        for b in reads:
            if a != b:
                olen = overlap(a, b, min_length)
                if olen > 0:
                    graph[a].append((b, olen))
    return graph


def draw_graph(graph: dict, figure_name: str) -> None:
    """Visualise the overlap graph with NetworkX / Matplotlib.

    Skips rendering when the graph exceeds 250 nodes to avoid freezing.

    Args:
        graph:       Overlap graph adjacency list (from build_overlap_graph).
        figure_name: Used as the plot title and saved PNG filename (no extension).
    """
    G = nx.DiGraph()
    for src, neighbours in graph.items():
        for dest, weight in neighbours:
            G.add_edge(src, dest, weight=weight)

    if len(G.nodes) > 250:
        print(
            f"Graph too large ({len(G.nodes)} nodes) – "
            f"skipping overlap graph for '{figure_name}'."
        )
        return

    pos    = nx.spring_layout(G, seed=42)
    labels = nx.get_edge_attributes(G, "weight")

    fig, ax = plt.subplots(figsize=(10, 7))
    nx.draw(G, pos, ax=ax, with_labels=False, node_size=500)
    nx.draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=labels)
    ax.set_title(figure_name.replace("_", " ").title())
    plt.savefig(f"{figure_name}.png", dpi=150)
    plt.close(fig)
    print(f"Overlap graph saved to '{figure_name}.png'")

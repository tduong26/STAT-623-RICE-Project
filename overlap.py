import networkx as nx
import matplotlib.pyplot as plt

def overlap(a, b, min_length=3):
    """
    Checks if the end of string 'a' overlaps with the beginning of string 'b'.

    Args:
        a (str): The first DNA sequence.
        b (str): The second DNA sequence.
        min_length (int): The minimum required overlap length.

    Returns:
        int: The length of the overlap if found, otherwise 0.
    """
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def build_overlap_graph(reads, min_length=3):
    """
    Constructs an adjacency list representing the overlap graph of the given reads.

    Args:
        reads (list): A list of DNA read sequences.
        min_length (int): The minimum overlap length to consider as an edge.

    Returns:
        dict: A dictionary where keys are reads and values are lists of
              (neighbor_read, overlap_length).
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

def draw_graph(graph, figure_name):
    """
    Visualizes the overlap graph using NetworkX and Matplotlib.

    Args:
        graph (dict): The overlap graph adjacency list.
        figure_name (str): The name to use for the title and saved file.
    """
    G = nx.DiGraph()

    for src in graph:
        for dest, weight in graph[src]:
            G.add_edge(src, dest, weight=weight)

    pos = nx.spring_layout(G)
    labels = nx.get_edge_attributes(G, 'weight')

    # Disable plot rendering if node count is excessive (would freeze terminal)
    if len(G.nodes) > 250:
        print(
            f"Graph too large ({len(G.nodes)} nodes) to reasonably draw plot, "
            f"skipping overlap graph rendering for {figure_name} (would otherwise freeze)."
        )
        return

    fig, ax = plt.subplots(figsize=(10, 7))

    nx.draw(G, pos, ax=ax, with_labels=False, node_size=500)
    nx.draw_networkx_edge_labels(G, pos, ax=ax, edge_labels=labels)

    ax.set_title(figure_name.replace("_", " ").title())

    plt.savefig(f"{figure_name}.png")
    plt.close(fig)

# STAT 623 — Shotgun Sequencing Simulation

Simulates the Lander-Waterman model of shotgun genome sequencing, comparing
assembly outcomes on a purely random genome vs a repeat-rich genome.

---

## Project structure

```
STAT623_project/
│
├── shortgun_simulation.py      ← main entry point  (run this)
│
├── theoretical_stats.py        ← α, E[# contigs], E[contig length]
├── genome_utils.py             ← genome + read generation
├── overlap.py                  ← overlap graph construction & drawing
├── alignment_visualization.py  ← pysam-free alignment, contigs, coverage
├── summary_plots.py            ← 6-panel comparison figure
├── extended_plots.py           ← 5 recommended statistical plots + α sweep
│
├── requirements.txt
└── README.md
```

---

## Installation

```bash
pip install -r requirements.txt
```

> **pysam has been removed.** The original code used pysam (Linux-only).
> All alignment bookkeeping now uses plain Python dicts — functionally
> identical, but fully compatible with Windows and Google Colab.

---

## Usage

```bash
python shortgun_simulation.py
```

This runs both simulations and saves **9 PNG files**:

| File | Contents |
|---|---|
| `alignment_plot_without_repeats.png` | Read alignment + coverage track (random genome) |
| `alignment_plot_with_repeats.png` | Read alignment + coverage track (repeat genome) |
| `overlap_structure_graph_without_repeats.png` | Overlap graph (random genome) |
| `overlap_structure_graph_with_repeats.png` | Overlap graph (repeat genome) |
| `contigs_only_without_repeats.png` | Contig map only (random genome) |
| `contigs_only_with_repeats.png` | Contig map only (repeat genome) |
| `summary_comparison.png` | 6-panel side-by-side summary |
| `lander_waterman_sweep.png` | High-resolution α-sweep (standalone) |
| `extended_summary.png` | 5 recommended statistical plots |

---

## Simulation parameters

All parameters are in `shortgun_simulation.py` inside the two simulation
functions. Key values:

| Parameter | Default | Description |
|---|---|---|
| `G` | 2000 bp | Genome length |
| `N` | 1200 | Number of reads |
| `min_len` | 50 bp | Minimum read length |
| `max_len` | 100 bp | Maximum read length (= L in theory) |
| `min_length` | 5 bp | Minimum overlap to create a graph edge |

---

## Theoretical background

The Lander-Waterman model assumes read start positions are
uniform over U(0, G), so read counts in any interval follow a Poisson
distribution with rate α = NL/G.

| Formula | Meaning |
|---|---|
| α = NL/G | Read density (coverage parameter) |
| E[# contigs] = N · e^{−α} | Expected number of assembled contigs |
| E[contig length] = G·(1−e^{−α}) / (N·e^{−α}) | Expected length per contig |
| P(covered) = 1 − e^{−α} | Fraction of genome with ≥ 1 read |

**Effect of repeats:** When repeat elements are longer than the read length L,
the uniform-start assumption breaks. Reads from multiple genomic positions map
identically, causing empirical contig counts to exceed N·e^{−α} and per-position
coverage to deviate from Poisson(α).

---

## Extended plots (`extended_plots.py`)

`generate_extended_summary_plots(s0, s1)` produces a 2×3 figure:

| Panel | Plot | Theory tested |
|---|---|---|
| ① | Lander-Waterman α sweep | E[# contigs] = N·e^{−α} |
| ② | Coverage depth vs Poisson(α) PMF | Poisson depth assumption |
| ③ | Genome covered fraction vs α | 1 − e^{−α} completeness |
| ④ | Gap-length distribution vs Exponential | Inter-gap distribution |
| ⑤ | N50 / N90 assembly quality | Assembly fragmentation |
| ⑥ | Contig-length CDF | Cumulative assembly completeness |

# STAT 623 — Shotgun Sequencing Simulation

Simulates the Lander-Waterman model of shotgun genome sequencing, comparing
assembly outcomes on a purely random genome vs a repeat-rich genome.

> **Course:** STAT 623 — Rice University  
> **Environment:** Jupyter Notebook / Google Colab  

---

## Project structure

```
STAT623_project/
│
├── shortgun_simulation.ipynb      ← ▶ START HERE — runs the full simulation
│
├── theoretical_stats.ipynb        ← Lander-Waterman formulas (α, E[contigs], E[length])
├── genome_utils.ipynb             ← Genome generation + shotgun read sampling
├── overlap.ipynb                  ← Overlap graph construction & drawing
├── alignment_visualization.ipynb  ← Alignment, contig detection, coverage & plots
├── summary_plots.ipynb            ← 6-panel side-by-side comparison figure
├── extended_plots.ipynb           ← 5 recommended statistical plots + α sweep
│
├── requirements.txt
└── README.md
```

---

## Notebook overview

| Notebook | Cells | What's inside |
|---|---|---|
| `shortgun_simulation.ipynb` | 13 | Main entry point — imports all modules, runs both simulations, saves all 9 output figures |
| `alignment_visualization.ipynb` | 12 | `create_alignment` · `get_contigs` · `get_coverage` · `visualize_alignment` · `plot_contigs` |
| `extended_plots.ipynb` | 12 | `compute_nx` · `extract_gaps` · `_empirical_contig_count` · `plot_lander_waterman_sweep` · `generate_extended_summary_plots` |
| `genome_utils.ipynb` | 8 | `generate_genome` · `generate_genome_with_repeats` · `generate_reads` |
| `overlap.ipynb` | 8 | `overlap` · `build_overlap_graph` · `draw_graph` |
| `summary_plots.ipynb` | 4 | `generate_summary_plots` |
| `theoretical_stats.ipynb` | 8 | `compute_alpha` · `expected_contigs` · `expected_contig_length` |

---

## Quick start

### Google Colab
1. Upload all `.ipynb` files to a single Colab folder (or clone this repo)
2. Open `shortgun_simulation.ipynb`
3. Uncomment and run the install cell at the top:
   ```python
   !pip install networkx matplotlib numpy
   ```
4. Run all cells — all 9 output PNGs will be saved to the working directory

### VS Code / JupyterLab (local)
```bash
# Install dependencies
pip install -r requirements.txt

# Open the main notebook
jupyter notebook shortgun_simulation.ipynb
```

> **Note:** `pysam` has been removed. The original code used pysam (Linux-only).
> All alignment bookkeeping now uses plain Python dicts — fully compatible
> with Windows, macOS, and Google Colab.

---

## Output files

Running `shortgun_simulation.ipynb` end-to-end produces **9 PNG files**:

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

All parameters are in `shortgun_simulation.ipynb` inside the two simulation
cells. Key values:

| Parameter | Default | Description |
|---|---|---|
| `G` | 2000 bp | Genome length |
| `N` | 1200 | Number of reads |
| `min_len` | 50 bp | Minimum read length |
| `max_len` | 100 bp | Maximum read length (= L in theory) |
| `min_length` | 5 bp | Minimum overlap to create a graph edge |

---

## Theoretical background

The Lander-Waterman model assumes read start positions are uniform over
U(0, G), so read counts in any interval follow a Poisson distribution
with rate α = NL/G.

| Formula | Meaning |
|---|---|
| α = NL/G | Read density (coverage parameter) |
| E[# contigs] = N · e^{−α} | Expected number of assembled contigs |
| E[contig length] = G·(1−e^{−α}) / (N·e^{−α}) | Expected length per contig |
| P(covered) = 1 − e^{−α} | Fraction of genome with ≥ 1 read |

**Effect of repeats:** When repeat elements are longer than the read length L,
the uniform-start assumption breaks. Reads from multiple genomic positions map
identically, causing empirical contig counts to exceed N·e^{−α} and
per-position coverage to deviate from Poisson(α).

---

## Extended plots (`extended_plots.ipynb`)

`generate_extended_summary_plots(s0, s1)` produces a 2×3 figure:

| Panel | Plot | Theory tested |
|---|---|---|
| ① | Lander-Waterman α sweep | E[# contigs] = N·e^{−α} |
| ② | Coverage depth vs Poisson(α) PMF | Poisson depth assumption |
| ③ | Genome covered fraction vs α | 1 − e^{−α} completeness |
| ④ | Gap-length distribution vs Exponential | Inter-gap distribution |
| ⑤ | N50 / N90 assembly quality | Assembly fragmentation |
| ⑥ | Contig-length CDF | Cumulative assembly completeness |

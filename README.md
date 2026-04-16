# STAT 623 Project — Shotgun Sequencing Simulation

This project simulates the shotgun sequencing process, builds an overlap graph
for genome assembly, and quantifies the effects of repetitive elements on contig
assembly, read mapping quality, and overlap graph complexity.

> **Platform note:** All four files are Windows-compatible. `pysam` has been
> replaced with a pure-Python in-memory alignment store — no BAM files are
> written to disk.

---

## Installation

Install the required packages with pip (no conda required):

```bash
pip install matplotlib numpy networkx
```

| Package | Purpose |
|---------|---------|
| `matplotlib` | All plots (alignments, contigs, coverage, sweeps) |
| `numpy` | Theoretical curves, bootstrap resampling |
| `networkx` | Overlap graph construction and visualisation |

`random`, `os`, and `math` are standard library — no install needed.

---

## Usage

Run all scenarios from a single entry point:

```bash
python shortgun_simulation.py
```

Run the alignment visualisation standalone test:

```bash
python alignment_visualization.py
```

---

## Project Files

| File | Role |
|------|------|
| `shortgun_simulation.py` | **Main entry point.** Genome generation, read sampling, all scenario runners, and the Lander-Waterman sweep |
| `alignment_visualization.py` | Alignment pipeline, contig/coverage extraction, effective-α metrics, bootstrap CI, and all plot functions |
| `theoretical_stats.py` | Lander-Waterman helper functions (α, expected contigs, expected contig length) |
| `overlap.py` | Overlap graph construction and NetworkX visualisation |

---

## File Reference

### `theoretical_stats.py`

Stateless helper functions implementing the Lander-Waterman model.

| Function | Description |
|----------|-------------|
| `compute_alpha(N, L, G)` | Read density α = N × L / G |
| `expected_contigs(N, alpha)` | Expected contig count ≈ N e^{−α} |
| `expected_contig_length(G, N, alpha)` | Expected contig length from the Lander-Waterman formula |

---

### `overlap.py`

Builds and draws the overlap graph from a list of DNA reads.

| Function | Description |
|----------|-------------|
| `overlap(a, b, min_length)` | Returns the length of the suffix-prefix overlap between reads `a` and `b`, or 0 if none |
| `build_overlap_graph(reads, min_length)` | Constructs a directed adjacency-list graph over all read pairs; edges exist where overlap ≥ `min_length` |
| `draw_graph(graph, figure_name)` | Renders the graph with NetworkX/Matplotlib and saves it as a PNG. Skips rendering automatically if the graph exceeds 250 nodes to avoid freezing |

---

### `alignment_visualization.py`

Contains the full alignment pipeline (pysam-free), all contig/coverage
extraction logic, the new statistical metrics, and every plot function.

#### In-memory alignment store (pysam replacement)

| Class / Function | Description |
|-----------------|-------------|
| `AlignedRead` | Lightweight stand-in for `pysam.AlignedSegment` (stores name, sequence, start, end, MAPQ) |
| `InMemoryBAM` | Collects `AlignedRead` objects in memory instead of writing a `.bam` file |
| `create_alignment_bam(genome, reads, key)` | Aligns reads against the genome with up to 2 mismatches; simulates MAPQ scoring and probabilistic read dropping; stores result under `key` in the in-memory store |
| `_open_bam(key)` / `_remove_bam(key)` | Retrieve / discard an in-memory BAM by key |

#### Core extraction

| Function | Description |
|----------|-------------|
| `find_with_mismatches(genome, read, max_mismatch)` | Returns all alignment positions where `read` matches `genome` with ≤ `max_mismatch` substitutions |
| `get_contigs(reads_data)` | Merges overlapping aligned reads into contigs; returns a list of `(start, end)` tuples |
| `get_coverage(genome_length, reads_data)` | Computes per-position read depth across the genome |

#### Statistical metrics (new additions)

| Function | Description |
|----------|-------------|
| `effective_alpha(reads_data, genome_len, L_mean)` | Computes **nominal α** (all aligned reads) and **effective α** (MAPQ-60 uniquely-mapped reads only), plus the **MAPQ loss fraction** — the share of coverage that is unusable due to multi-mapping |
| `bootstrap_ci(trial_counts, n_boot=500)` | Bootstrap 95% confidence interval for the mean of a list of trial counts; returns `(2.5th percentile, 97.5th percentile)` |

#### Plot functions

| Function | Output file | Description |
|----------|-------------|-------------|
| `visualize_alignment(genome, reads, output_png)` | e.g. `alignment_nr_scenario1.png` | Two-panel plot: top panel shows individual reads coloured by MAPQ (steelblue = unique, orange = multi-mapper, gray = ambiguous) with contig bars; bottom panel shows the read depth track |
| `plot_contigs(ax, contigs, ...)` | — | Draws contig blocks on a matplotlib axis; optionally labels each contig and draws dashed gap lines |
| `plot_lander_waterman_sweep(...)` | e.g. `lander_waterman_sweep_s1.png` | Two-panel figure: **left** = no-repeat genome (blue stems vs gray dashed theory curve); **right** = repeat genome (red stems vs theory curve). Bootstrap 95% CI bands shown when trial data is provided |
| `plot_effective_alpha_comparison(results, output_png)` | `effective_alpha_comparison.png` | Side-by-side bar chart: nominal α vs effective α per condition (left panel), and MAPQ loss % per condition (right panel) |

---

### `shortgun_simulation.py`

Main entry point. Imports from all three other files.

#### Genome generators

| Function | Description |
|----------|-------------|
| `generate_genome(length)` | Generates a purely random DNA sequence with no structure |
| `generate_genome_with_repeats(length, snp_rate, repeat_sizes, weights)` | Generates a bio-realistic genome with four repeat families (LINE, SINE, LTR, DNA transposons), tandem repeats, and unique sequence. Each family uses a single source template; copies are mutated at `snp_rate` |
| `generate_genome_with_repeat_fraction(length, repeat_fraction, repeat_unit_len, snp_rate)` | Generates a genome where **exactly** `repeat_fraction` of the sequence is copies of one repeat template and the rest is unique. Used by the Part (c) density sweep for precise experimental control |

#### Read generator

| Function | Description |
|----------|-------------|
| `generate_reads(genome, num_reads, min_len, max_len)` | Simulates shotgun sequencing: samples `num_reads` random substrings of length uniformly drawn from `[min_len, max_len]` |

#### Scenario runners

| Function | Output files | Description |
|----------|-------------|-------------|
| `display_simulation_results(label, genome, reads, bam_key, output_png)` | e.g. `contigs_s1_nr.png` | Aligns reads, prints empirical stats vs Lander-Waterman theoretical predictions, prints nominal and effective α, saves a contig-only plot. Returns `reads_data` for downstream use |
| `run_alignment_comparison()` | `alignment_nr_scenario1.png`, `alignment_r_scenario1.png`, `contigs_s1_nr.png`, `contigs_s1_r.png` | **Scenario 1** — high coverage (N=600, G=4000 bp). Compares a non-repeat genome against a realistic repeat genome (LINE 1000 bp, SINE 300 bp, LTR 500 bp). Highlights MAPQ colour differences |
| `run_assembly_comparison()` | `alignment_nr_scenario2.png`, `alignment_r_scenario2.png`, `contigs_s2_nr.png`, `contigs_s2_r.png` | **Scenario 2** — intermediate coverage (N=150, G=2000 bp). Uses heavy repeat weights (LINE 40%, SINE 30%) to demonstrate contig fragmentation |
| `run_lander_waterman_sweep(genome_r, output_png, label, num_trials)` | `lander_waterman_sweep_s1.png`, `lander_waterman_sweep_s2.png` | For a given repeat genome, sweeps α from 0.2 to 6.0, runs `num_trials` trials per point, and plots empirical contig counts vs the Lander-Waterman theoretical curve for both repeat and non-repeat genomes. Bootstrap 95% CI bands included |
| `run_repeat_density_sweep(output_png, num_trials)` | `repeat_density_sweep.png` | **Part (c)** — fixes G=2000, N=150, L∈[50,100] and sweeps repeat fraction from 0% to 90% in 10 steps. Measures four quantities at each level: contig count, overlap graph edge count, effective α, and MAPQ loss. Produces a 4-panel figure with CI bands |

---

## Output Files

Running `python shortgun_simulation.py` produces the following PNG files in
the working directory:

| File | Produced by | Description |
|------|-------------|-------------|
| `alignment_nr_scenario1.png` | `run_alignment_comparison` | Full alignment plot, non-repeat genome, Scenario 1 |
| `alignment_r_scenario1.png` | `run_alignment_comparison` | Full alignment plot, repeat genome, Scenario 1 |
| `contigs_s1_nr.png` | `run_alignment_comparison` | Contig-only plot, non-repeat, Scenario 1 |
| `contigs_s1_r.png` | `run_alignment_comparison` | Contig-only plot, repeat, Scenario 1 |
| `alignment_nr_scenario2.png` | `run_assembly_comparison` | Full alignment plot, non-repeat genome, Scenario 2 |
| `alignment_r_scenario2.png` | `run_assembly_comparison` | Full alignment plot, repeat genome, Scenario 2 |
| `contigs_s2_nr.png` | `run_assembly_comparison` | Contig-only plot, non-repeat, Scenario 2 |
| `contigs_s2_r.png` | `run_assembly_comparison` | Contig-only plot, repeat, Scenario 2 |
| `effective_alpha_comparison.png` | `main` | Nominal vs effective α bar chart across all 4 conditions |
| `lander_waterman_sweep_s1.png` | `run_lander_waterman_sweep` | Theory vs empirical contig counts, Scenario 1 genome, with CI bands |
| `lander_waterman_sweep_s2.png` | `run_lander_waterman_sweep` | Theory vs empirical contig counts, Scenario 2 genome, with CI bands |
| `repeat_density_sweep.png` | `run_repeat_density_sweep` | Part (c) 4-panel sweep: contigs, graph edges, effective α, MAPQ loss vs repeat fraction |

---

## Theoretical Background

### Lander-Waterman Model

| Symbol | Meaning |
|--------|---------|
| G | Genome length (bp) |
| N | Number of reads |
| L | Read length (mean or max) |
| α | Read density: α = N × L / G |

Key predictions:

```
E[# contigs]       = N · e^{−α}   ≈ (G/L) · α · e^{−α}
E[contig length]   = G · (1 − e^{−α}) / (N · e^{−α})
```

### Effective Coverage (new metric)

The **nominal α** counts all aligned reads. The **effective α** counts only
uniquely-mapped reads (MAPQ = 60):

```
α_eff = N_unique × L_mean / G
MAPQ loss = 1 − (α_eff / α_nominal)
```

A genome with high repeat content can have a nominal α that looks adequate
while the effective α is dramatically lower — the paper's core verbal claim,
now quantified.

### Bootstrap 95% Confidence Intervals

For each alpha value in a sweep, `num_trials` independent trials are run.
The CI is computed by resampling those trial counts 500 times with replacement
and taking the 2.5th and 97.5th percentiles of the bootstrapped means. The
wider CI bands on the repeat-genome panels directly show that repeat structure
increases **variance** in assembly outcomes, not just the mean.

---

## Simulation Design

1. **Genome generation** — `generate_genome_with_repeats` builds a genome by
   interleaving blocks drawn from four repeat families and unique random
   sequence. Each family uses a single source template so copies are highly
   identical; `snp_rate` controls how much each copy drifts from the template.

2. **Alignment pipeline** — `find_with_mismatches` allows up to 2
   substitutions per read position. Reads with one unique match get MAPQ 60;
   reads with 2–5 matches get MAPQ 20 and are dropped with 80% probability;
   reads with 6+ matches get MAPQ 0 and are dropped with 90% probability.
   This mirrors real short-read aligner behaviour.

3. **Overlap graph** — built with a minimum overlap of 5 bp (3 bp default);
   rendered via NetworkX. Rendering is skipped automatically when the graph
   exceeds 250 nodes.

4. **Part (c) repeat sweep** — `generate_genome_with_repeat_fraction`
   places repeat blocks and unique blocks in alternating fashion so repeat
   content is spread uniformly across the genome, giving a controlled
   experimental variable.

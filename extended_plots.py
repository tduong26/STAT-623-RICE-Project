"""
extended_plots.py
-----------------
Five recommended statistical plots that go beyond the basic summary.

① Lander-Waterman α sweep       — theory curve vs empirical scatter (both genomes)
② Coverage vs Poisson(α) PMF   — tests the Poisson depth assumption
③ Coverage fraction vs α        — completeness curve 1 − e^{-α}
④ Gap-length distribution        — histogram + fitted Exponential PDF
⑤ N50 / N90 assembly quality    — standard bioinformatics metric
⑥ Contig-length CDF              — step-CDF with N50/N90 intersection markers

Helper functions
----------------
compute_nx(contig_lengths, x)    — N50 when x=0.5, N90 when x=0.9
extract_gaps(contigs)            — gap widths between consecutive contigs
_empirical_contig_count(...)     — internal helper for the α sweep
plot_lander_waterman_sweep(...)  — standalone high-resolution α-sweep figure
"""

import math

import matplotlib.pyplot as plt
import numpy as np

from alignment_visualization import create_alignment, get_contigs
from genome_utils import generate_genome, generate_genome_with_repeats, generate_reads
from theoretical_stats import compute_alpha


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def compute_nx(contig_lengths: list, x: float) -> int:
    """Compute the Nx assembly metric (N50 when x=0.5, N90 when x=0.9).

    Nx is the length L* such that contigs of length ≥ L* account for at
    least a fraction x of the total assembled sequence.

    Args:
        contig_lengths: List of contig lengths in bp.
        x:              Fraction threshold, e.g. 0.5 for N50.
    Returns:
        Nx value in bp, or 0 if the list is empty.
    """
    if not contig_lengths:
        return 0
    sorted_l = sorted(contig_lengths, reverse=True)
    target, cumsum = sum(sorted_l) * x, 0
    for l in sorted_l:
        cumsum += l
        if cumsum >= target:
            return l
    return sorted_l[-1]


def extract_gaps(contigs: list) -> list:
    """Return gap lengths (bp) between consecutive contigs.

    Args:
        contigs: List of (start, end) tuples in ascending order.
    Returns:
        List of positive gap sizes in bp.
    """
    return [
        contigs[i + 1][0] - contigs[i][1]
        for i in range(len(contigs) - 1)
        if contigs[i + 1][0] > contigs[i][1]
    ]


def _empirical_contig_count(genome: str, num_reads: int,
                             min_len: int, max_len: int) -> tuple:
    """Run one sequencing trial and return (alpha, empirical_contig_count).

    Internal helper used by the α-sweep functions.

    Args:
        genome:    Reference genome to sample reads from.
        num_reads: Number of reads for this trial.
        min_len:   Minimum read length.
        max_len:   Maximum read length (also used as L for alpha).
    Returns:
        (alpha, n_contigs) tuple.
    """
    reads      = generate_reads(genome, num_reads, min_len, max_len)
    reads_data = create_alignment(genome, reads)
    contigs    = get_contigs(reads_data)
    alpha      = compute_alpha(num_reads, max_len, len(genome))
    return alpha, len(contigs)


# ---------------------------------------------------------------------------
# Standalone Lander-Waterman sweep (high-resolution, two-panel figure)
# ---------------------------------------------------------------------------

def plot_lander_waterman_sweep(
        G: int = 2000,
        min_len: int = 50,
        max_len: int = 100,
        n_trials: int = 3,
        output_png: str = "lander_waterman_sweep.png") -> None:
    """Plot the Lander-Waterman contig-count curve against empirical data.

    Theory: E[# contigs] = N · e^{-α} = (G/L) · α · e^{-α}

    Sweeps N so that α spans 0.2 → 6.0. Each point is averaged over
    n_trials independent trials to reduce sampling noise. Both genome types
    are plotted so the deviation caused by repeats is clearly visible.

    Args:
        G:          Genome length in bp.
        min_len:    Minimum read length.
        max_len:    Maximum read length (= L in the α formula).
        n_trials:   Independent trials per N value (increase for smoother dots).
        output_png: Output filename.
    """
    L             = max_len
    alpha_targets = np.linspace(0.2, 6.0, 18)
    N_values      = np.maximum(np.round(alpha_targets * G / L).astype(int), 1)
    alpha_curve   = np.linspace(0, 6.5, 300)

    print("Generating genomes for α sweep …")
    genome_plain  = generate_genome(G)
    genome_repeat = generate_genome_with_repeats(G)

    emp_pa, emp_pc = [], []
    emp_ra, emp_rc = [], []

    for N in N_values:
        cp_acc = cr_acc = 0.0
        for _ in range(n_trials):
            _, cp = _empirical_contig_count(genome_plain,  N, min_len, max_len)
            _, cr = _empirical_contig_count(genome_repeat, N, min_len, max_len)
            cp_acc += cp
            cr_acc += cr
        alpha_n = N * L / G
        emp_pa.append(alpha_n);  emp_pc.append(cp_acc / n_trials)
        emp_ra.append(alpha_n);  emp_rc.append(cr_acc / n_trials)
        print(f"  N={N:5d}  α={alpha_n:.2f}  "
              f"plain={cp_acc/n_trials:.1f}  repeat={cr_acc/n_trials:.1f}")

    C0, C1 = "#2196F3", "#F44336"
    fig, axes = plt.subplots(1, 2, figsize=(14, 5),
                             gridspec_kw={"wspace": 0.35})

    for ax, (emp_a, emp_c), color, label in [
        (axes[0], (emp_pa, emp_pc), C0, "No repeats"),
        (axes[1], (emp_ra, emp_rc), C1, "With repeats"),
    ]:
        theory_y = (G / L) * alpha_curve * np.exp(-alpha_curve)
        ax.plot(alpha_curve, theory_y, color="#888", lw=2, ls="--",
                label=r"Theory: $(G/L)\,\alpha\,e^{-\alpha}$", zorder=2)
        ax.scatter(emp_a, emp_c, color=color, s=55, zorder=3,
                   edgecolors="white", lw=0.6, label=f"Empirical ({label})")

        theory_at_emp = (G / L) * np.array(emp_a) * np.exp(-np.array(emp_a))
        for xa, yc, yt in zip(emp_a, emp_c, theory_at_emp):
            ax.plot([xa, xa], [yt, yc], color=color, alpha=0.25, lw=1, zorder=1)

        mid_alpha = N_values[len(N_values) // 2] * L / G
        ax.annotate(f"α = {mid_alpha:.2f}\n(typical run)",
                    xy=(mid_alpha, (G/L)*mid_alpha*math.exp(-mid_alpha)),
                    xytext=(mid_alpha + 0.6, (G/L)*mid_alpha*math.exp(-mid_alpha) + (G/L)*0.07),
                    fontsize=8, color="#555",
                    arrowprops=dict(arrowstyle="->", color="#999", lw=0.8,
                                   connectionstyle="arc3,rad=0.2"))

        ax.set_xlabel(r"Read density  $\alpha = NL/G$", fontsize=10)
        ax.set_ylabel("Number of contigs", fontsize=10)
        ax.set_title(f"Lander-Waterman sweep — {label}", fontweight="bold")
        ax.legend(fontsize=8)
        ax.set_xlim(0, 6.8);  ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(True, ls="--", alpha=0.4)

    axes[1].text(0.97, 0.95,
                 "Dots above the curve:\nPoisson assumption violated\nby repeat structure",
                 transform=axes[1].transAxes, fontsize=8, color=C1,
                 ha="right", va="top",
                 bbox=dict(boxstyle="round,pad=0.4", fc="white",
                           ec=C1, alpha=0.85, lw=0.8))
    fig.suptitle(r"Lander-Waterman Model: $E[\#\,\mathrm{contigs}] = N e^{-\alpha}$",
                 fontsize=13, fontweight="bold")
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nLander-Waterman sweep saved to '{output_png}'")


# ---------------------------------------------------------------------------
# All-in-one extended summary (2 × 3 figure)
# ---------------------------------------------------------------------------

def generate_extended_summary_plots(
        s0: dict,
        s1: dict,
        output_png: str = "extended_summary.png",
        sweep_trials: int = 2) -> None:
    """Produce a 2 × 3 figure with all five recommended statistical plots.

    Args:
        s0:           Stats dict from shotgun_simulation_without_repeats().
        s1:           Stats dict from shotgun_simulation_with_repeats().
        output_png:   Output filename.
        sweep_trials: Trials per N value in the α-sweep panel
                      (2 = fast; raise to 4–5 for smoother empirical dots).
    """
    C0   = "#2196F3"
    C1   = "#F44336"
    GRAY = "#888888"

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle("Extended Statistical Analysis — Shotgun Sequencing Simulation",
                 fontsize=14, fontweight="bold", y=1.01)

    G  = s0["G"]
    L  = max(s0["read_lengths"])
    a0 = s0["alpha"]
    a1 = s1["alpha"]

    # ── [0,0]  ① Lander-Waterman α sweep ───────────────────────────────────
    ax            = axes[0, 0]
    alpha_targets = np.linspace(0.3, 5.5, 14)
    N_sweep       = np.maximum(np.round(alpha_targets * G / L).astype(int), 1)
    alpha_curve   = np.linspace(0, 6.2, 400)
    theory_y      = (G / L) * alpha_curve * np.exp(-alpha_curve)

    print("\n[Extended plots] Running α sweep …")
    pt_ap, pt_cp, pt_ar, pt_cr = [], [], [], []

    for N in N_sweep:
        cp_acc = cr_acc = 0
        for _ in range(sweep_trials):
            _, cp = _empirical_contig_count(s0["genome"], N, 50, L)
            _, cr = _empirical_contig_count(s1["genome"], N, 50, L)
            cp_acc += cp;  cr_acc += cr
        alpha_n = N * L / G
        pt_ap.append(alpha_n);  pt_cp.append(cp_acc / sweep_trials)
        pt_ar.append(alpha_n);  pt_cr.append(cr_acc / sweep_trials)

    ax.plot(alpha_curve, theory_y, color=GRAY, lw=2, ls="--",
            label=r"Theory: $(G/L)\,\alpha\,e^{-\alpha}$", zorder=2)
    ax.scatter(pt_ap, pt_cp, color=C0, s=50, zorder=4,
               edgecolors="white", lw=0.7, label="No repeats")
    ax.scatter(pt_ar, pt_cr, color=C1, s=50, zorder=4, marker="^",
               edgecolors="white", lw=0.7, label="With repeats")

    for a_op, col in [(a0, C0), (a1, C1)]:
        ax.axvline(a_op, color=col, lw=0.9, ls=":", alpha=0.65)

    theory_at_emp = (G / L) * np.array(pt_ar) * np.exp(-np.array(pt_ar))
    for xa, yc, yt in zip(pt_ar, pt_cr, theory_at_emp):
        ax.plot([xa, xa], [yt, yc], color=C1, alpha=0.20, lw=1.2, zorder=1)

    ax.text(0.97, 0.97,
            "Dots above curve:\nPoisson assumption\nviolated by repeats",
            transform=ax.transAxes, fontsize=7.5, color=C1,
            ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.35", fc="white",
                      ec=C1, alpha=0.85, lw=0.8))
    ax.set_xlabel(r"Read density  $\alpha = NL/G$", fontsize=9)
    ax.set_ylabel("Number of contigs", fontsize=9)
    ax.set_title("① Lander-Waterman α sweep", fontweight="bold")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(0, 6.5);  ax.set_ylim(bottom=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, ls="--", alpha=0.35)
    print("  α sweep complete.")

    # ── [0,1]  ② Coverage vs Poisson(α) PMF ────────────────────────────────
    ax        = axes[0, 1]
    bar_width = 0.38

    for offset, s, col, lbl in [
        (-bar_width / 2, s0, C0, "No repeats"),
        ( bar_width / 2, s1, C1, "With repeats"),
    ]:
        cov       = np.array(s["coverage"])
        alpha_val = s["alpha"]
        p99       = int(np.percentile(cov, 99)) + 1
        counts    = np.bincount(cov, minlength=p99 + 1)[:p99 + 1]
        k_vals    = np.arange(p99 + 1)
        ax.bar(k_vals + offset, counts / len(cov), width=bar_width,
               color=col, alpha=0.45, label=f"Empirical — {lbl}",
               edgecolor="none")
        log_pmf = (-alpha_val
                   + k_vals * np.log(alpha_val + 1e-300)
                   - np.array([math.lgamma(k + 1) for k in k_vals]))
        ax.plot(k_vals, np.exp(log_pmf), color=col, lw=1.8, ls="--",
                label=fr"Poisson($\alpha$={alpha_val:.2f})")

    ax.set_xlabel("Read depth at each position", fontsize=9)
    ax.set_ylabel("Probability", fontsize=9)
    ax.set_title("② Coverage distribution vs Poisson(α) PMF", fontweight="bold")
    ax.legend(fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, ls="--", alpha=0.35)

    # ── [0,2]  ③ Coverage fraction vs α ────────────────────────────────────
    ax          = axes[0, 2]
    alpha_fine  = np.linspace(0, 8, 500)
    theory_comp = (1 - np.exp(-alpha_fine)) * 100

    ax.plot(alpha_fine, theory_comp, color=GRAY, lw=2, ls="--",
            label=r"Theory: $1 - e^{-\alpha}$")
    for pct in [95, 99]:
        ax.axhline(pct, color="#CCCCCC", lw=1, ls=":")
        ax.text(7.85, pct + 0.8, f"{pct}%", fontsize=7.5,
                color="#AAAAAA", va="bottom")

    for s, col, lbl in [(s0, C0, "No repeats"), (s1, C1, "With repeats")]:
        cov     = np.array(s["coverage"])
        emp_pct = float(np.sum(cov > 0)) / len(cov) * 100
        ax.scatter([s["alpha"]], [emp_pct], color=col, s=90, zorder=5,
                   edgecolors="white", lw=0.8, label=f"Empirical — {lbl}")
        ax.annotate(f"{emp_pct:.1f}%",
                    xy=(s["alpha"], emp_pct),
                    xytext=(s["alpha"] + 0.2, emp_pct - 5),
                    fontsize=8, color=col,
                    arrowprops=dict(arrowstyle="->", color=col, lw=0.7,
                                   connectionstyle="arc3,rad=0.2"))

    alpha_99 = -math.log(0.01)
    for s, col in [(s0, C0), (s1, C1)]:
        ax.annotate("", xy=(alpha_99, 99), xytext=(s["alpha"], 99),
                    arrowprops=dict(arrowstyle="->", color=col, lw=0.7,
                                   linestyle="dotted"),
                    annotation_clip=False)

    ax.set_xlabel(r"$\alpha = NL/G$", fontsize=9)
    ax.set_ylabel("Genome covered (%)", fontsize=9)
    ax.set_title("③ Coverage fraction vs α", fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_xlim(0, 8.2);  ax.set_ylim(0, 106)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, ls="--", alpha=0.35)

    # ── [1,0]  ④ Gap-length distribution vs Exponential fit ────────────────
    ax           = axes[1, 0]
    any_gap_data = False

    for s, col, lbl in [(s0, C0, "No repeats"), (s1, C1, "With repeats")]:
        gaps = extract_gaps(s["contigs"])
        if not gaps:
            continue
        any_gap_data = True
        gaps_arr = np.array(gaps, dtype=float)
        mean_gap = gaps_arr.mean()
        lam      = 1.0 / mean_gap
        n_bins   = min(25, max(5, len(gaps)))
        ax.hist(gaps_arr, bins=n_bins, density=True, color=col, alpha=0.40,
                edgecolor="white", lw=0.5, label=f"Empirical — {lbl}")
        x_fit = np.linspace(0, gaps_arr.max() * 1.05, 400)
        ax.plot(x_fit, lam * np.exp(-lam * x_fit), color=col, lw=2, ls="--",
                label=fr"Exp fit ($\lambda$={lam:.4f}, μ={mean_gap:.0f} bp)")
        ax.axvline(mean_gap, color=col, lw=1, ls=":", alpha=0.75)

    if not any_gap_data:
        ax.text(0.5, 0.5, "No gaps detected\n(genome fully covered)",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=10, color=GRAY)

    ax.set_xlabel("Gap length (bp)", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_title("④ Gap-length distribution vs Exponential fit", fontweight="bold")
    ax.legend(fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, ls="--", alpha=0.35)

    # ── [1,1]  ⑤ N50 / N90 / Max / Mean bar chart ──────────────────────────
    ax      = axes[1, 1]
    metrics = ["N50", "N90", "Max length", "Mean length"]
    vals_0  = [
        compute_nx(s0["contig_lengths"], 0.50),
        compute_nx(s0["contig_lengths"], 0.90),
        max(s0["contig_lengths"]) if s0["contig_lengths"] else 0,
        int(sum(s0["contig_lengths"]) / max(len(s0["contig_lengths"]), 1)),
    ]
    vals_1  = [
        compute_nx(s1["contig_lengths"], 0.50),
        compute_nx(s1["contig_lengths"], 0.90),
        max(s1["contig_lengths"]) if s1["contig_lengths"] else 0,
        int(sum(s1["contig_lengths"]) / max(len(s1["contig_lengths"]), 1)),
    ]

    x_pos = np.arange(len(metrics))
    width = 0.35
    b0 = ax.bar(x_pos - width/2, vals_0, width, color=C0,
                label="No repeats",   edgecolor="white", lw=0.8)
    b1 = ax.bar(x_pos + width/2, vals_1, width, color=C1,
                label="With repeats", edgecolor="white", lw=0.8)

    for bars in (b0, b1):
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(bar.get_x() + bar.get_width() / 2,
                        h + max(vals_0 + vals_1) * 0.01,
                        f"{h:,}", ha="center", va="bottom", fontsize=7)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(metrics, fontsize=9)
    ax.set_ylabel("Contig length (bp)", fontsize=9)
    ax.set_title("⑤ N50 / N90 Assembly Quality", fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, axis="y", ls="--", alpha=0.35)

    # ── [1,2]  ⑥ Contig-length CDF ─────────────────────────────────────────
    ax        = axes[1, 2]
    x_max_all = 1

    for s, col, lbl in [(s0, C0, "No repeats"), (s1, C1, "With repeats")]:
        cl = sorted(s["contig_lengths"])
        if not cl:
            continue
        cl_arr   = np.array(cl, dtype=float)
        cdf_vals = np.arange(1, len(cl_arr) + 1) / len(cl_arr) * 100
        ax.step(cl_arr, cdf_vals, color=col, lw=2, where="post",
                label=f"Empirical — {lbl}")
        x_max_all = max(x_max_all, cl_arr.max())

        for nx_val, pct, mk in [
            (compute_nx(cl, 0.50), 50, "o"),
            (compute_nx(cl, 0.90), 90, "^"),
        ]:
            ax.scatter([nx_val], [pct], color=col, s=55, marker=mk,
                       zorder=5, edgecolors="white", lw=0.7)
            ax.annotate(f"N{'50' if pct==50 else '90'}={nx_val:,}",
                        xy=(nx_val, pct),
                        xytext=(nx_val + x_max_all * 0.04, pct - 6),
                        fontsize=7, color=col,
                        arrowprops=dict(arrowstyle="->", color=col, lw=0.6))

    for pct, lbl_txt in [(50, "N50 = 50%"), (90, "N90 = 90%")]:
        ax.axhline(pct, color="#CCCCCC", lw=1, ls=":")
        ax.text(x_max_all * 1.01, pct + 1, lbl_txt,
                fontsize=7, color="#AAAAAA", va="bottom")

    ax.set_xlabel("Contig length (bp)", fontsize=9)
    ax.set_ylabel("Cumulative fraction (%)", fontsize=9)
    ax.set_title("⑥ Contig-length CDF", fontweight="bold")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 106)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, ls="--", alpha=0.35)

    plt.tight_layout()
    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nExtended summary saved to '{output_png}'")

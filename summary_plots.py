"""
summary_plots.py
----------------
Six-panel overview figure comparing the two simulation runs side-by-side.

Panels
------
[0,0]  Contig map comparison        — both genomes drawn to scale
[0,1]  Empirical vs theoretical     — grouped bar: N·e^{-α} vs observed
[0,2]  Mapped-read fraction         — stacked bar showing unmapped reads
[1,0]  Coverage depth profile       — overlaid line + fill
[1,1]  Contig-length distribution   — overlaid histogram with mean lines
[1,2]  Key metrics summary table    — all headline numbers in one view
"""

import numpy as np
import matplotlib.pyplot as plt


def generate_summary_plots(s0: dict, s1: dict,
                            output_png: str = "summary_comparison.png") -> None:
    """Produce a 2 × 3 summary figure comparing both simulations.

    Args:
        s0:         Stats dict returned by shotgun_simulation_without_repeats().
        s1:         Stats dict returned by shotgun_simulation_with_repeats().
        output_png: Output filename for the saved PNG.
    """
    C0   = "#2196F3"   # blue  — no repeats
    C1   = "#F44336"   # red   — with repeats
    GRAY = "#9E9E9E"   # theory / neutral bars

    labels = ["No Repeats", "With Repeats"]
    stats  = [s0, s1]
    colors = [C0, C1]

    fig = plt.figure(figsize=(18, 11))
    fig.suptitle("Shotgun Sequencing Simulation — Summary & Comparison",
                 fontsize=15, fontweight="bold", y=0.99)
    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.38)

    # ------------------------------------------------------------------
    # [0, 0]  Contig map
    # ------------------------------------------------------------------
    ax_map = fig.add_subplot(gs[0, 0])
    G = s0["G"]

    for y, s, col, lbl in zip([1.6, 0.6], stats, colors, ["No repeats", "With repeats"]):
        ax_map.plot([0, G], [y, y], color="black", linewidth=2, zorder=1)
        ax_map.text(-60, y, lbl, va="center", ha="right", fontsize=8)
        for start, end in s["contigs"]:
            ax_map.plot([start, end], [y, y], color=col, linewidth=10,
                        solid_capstyle="butt", zorder=2, alpha=0.85)

    ax_map.set_xlim(-80, G + 20)
    ax_map.set_ylim(0.1, 2.1)
    ax_map.set_xlabel("Genome Position (bp)", fontsize=9)
    ax_map.set_title("Contig Map Comparison", fontweight="bold")
    ax_map.set_yticks([])
    ax_map.spines[["left", "top", "right"]].set_visible(False)

    # ------------------------------------------------------------------
    # [0, 1]  Empirical vs theoretical contig count
    # ------------------------------------------------------------------
    ax_cnt = fig.add_subplot(gs[0, 1])
    x     = np.arange(2)
    width = 0.32
    emp    = [len(s["contigs"])  for s in stats]
    theory = [s["exp_contigs"]   for s in stats]

    bars_emp = ax_cnt.bar(x - width/2, emp,    width, color=colors,
                          label="Empirical", edgecolor="white", linewidth=0.8)
    bars_thy = ax_cnt.bar(x + width/2, theory, width, color=GRAY,
                          label="Theoretical (N·e⁻ᵅ)", edgecolor="white",
                          linewidth=0.8, hatch="//")

    for bar in bars_emp:
        ax_cnt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f"{bar.get_height():.0f}",
                    ha="center", va="bottom", fontsize=8, fontweight="bold")
    for bar in bars_thy:
        ax_cnt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    f"{bar.get_height():.1f}",
                    ha="center", va="bottom", fontsize=8, color="#555")

    ax_cnt.set_xticks(x)
    ax_cnt.set_xticklabels(labels, fontsize=9)
    ax_cnt.set_ylabel("Number of Contigs", fontsize=9)
    ax_cnt.set_title("Empirical vs Theoretical\nContig Count", fontweight="bold")
    ax_cnt.legend(fontsize=8, loc="upper right")
    ax_cnt.spines[["top", "right"]].set_visible(False)

    # ------------------------------------------------------------------
    # [0, 2]  Mapped-read fraction
    # ------------------------------------------------------------------
    ax_mp = fig.add_subplot(gs[0, 2])
    mapped_pct   = [100 * s["mapped_reads"] / s["N"] for s in stats]
    unmapped_pct = [100 - p for p in mapped_pct]
    bottom = np.zeros(2)

    for pct, col, lbl in zip([mapped_pct, unmapped_pct],
                               [colors, ["#BDBDBD", "#BDBDBD"]],
                               ["Mapped", "Unmapped"]):
        ax_mp.bar(labels, pct, bottom=bottom, color=col, label=lbl,
                  edgecolor="white", linewidth=0.8)
        bottom += np.array(pct)

    for i, (p, col) in enumerate(zip(mapped_pct, colors)):
        ax_mp.text(i, p / 2, f"{p:.1f}%", ha="center", va="center",
                   fontsize=10, fontweight="bold", color="white")

    ax_mp.set_ylim(0, 115)
    ax_mp.set_ylabel("Percentage of Reads (%)", fontsize=9)
    ax_mp.set_title("Mapped vs Unmapped Reads", fontweight="bold")
    ax_mp.legend(fontsize=8, loc="upper right")
    ax_mp.spines[["top", "right"]].set_visible(False)

    # ------------------------------------------------------------------
    # [1, 0]  Coverage depth profile
    # ------------------------------------------------------------------
    ax_cov = fig.add_subplot(gs[1, 0])
    for s, col, lbl in zip(stats, colors, labels):
        cov = s["coverage"]
        ax_cov.plot(range(len(cov)), cov, color=col, linewidth=1,
                    alpha=0.85, label=lbl)
        ax_cov.fill_between(range(len(cov)), cov, alpha=0.15, color=col)

    ax_cov.set_xlabel("Genome Position (bp)", fontsize=9)
    ax_cov.set_ylabel("Read Depth", fontsize=9)
    ax_cov.set_title("Coverage Depth Profile", fontweight="bold")
    ax_cov.legend(fontsize=8)
    ax_cov.spines[["top", "right"]].set_visible(False)

    # ------------------------------------------------------------------
    # [1, 1]  Contig-length distribution
    # ------------------------------------------------------------------
    ax_hist = fig.add_subplot(gs[1, 1])
    for s, col, lbl in zip(stats, colors, labels):
        cl = s["contig_lengths"]
        if cl:
            ax_hist.hist(cl, bins=20, color=col, alpha=0.55,
                         edgecolor="white", label=lbl)
            ax_hist.axvline(sum(cl) / len(cl), color=col, linestyle="--",
                            linewidth=1.5, label=f"{lbl} mean")

    ax_hist.set_xlabel("Contig Length (bp)", fontsize=9)
    ax_hist.set_ylabel("Frequency", fontsize=9)
    ax_hist.set_title("Contig-Length Distribution", fontweight="bold")
    ax_hist.legend(fontsize=7)
    ax_hist.spines[["top", "right"]].set_visible(False)

    # ------------------------------------------------------------------
    # [1, 2]  Key metrics table
    # ------------------------------------------------------------------
    ax_tbl = fig.add_subplot(gs[1, 2])
    ax_tbl.axis("off")

    rows = [
        ["Metric",               "No Repeats",                          "With Repeats"],
        ["Genome length (bp)",   f"{s0['G']}",                          f"{s1['G']}"],
        ["Reads (N)",            f"{s0['N']}",                          f"{s1['N']}"],
        ["Read density (α)",     f"{s0['alpha']:.3f}",                  f"{s1['alpha']:.3f}"],
        ["Mapped reads",
         f"{s0['mapped_reads']} ({100*s0['mapped_reads']/s0['N']:.1f}%)",
         f"{s1['mapped_reads']} ({100*s1['mapped_reads']/s1['N']:.1f}%)"],
        ["Contigs (empirical)",  f"{len(s0['contigs'])}",               f"{len(s1['contigs'])}"],
        ["Contigs (theory)",     f"{s0['exp_contigs']:.1f}",            f"{s1['exp_contigs']:.1f}"],
        ["Mean contig len (bp)",
         f"{sum(s0['contig_lengths'])/len(s0['contig_lengths']):.1f}" if s0['contig_lengths'] else "—",
         f"{sum(s1['contig_lengths'])/len(s1['contig_lengths']):.1f}" if s1['contig_lengths'] else "—"],
        ["Max contig len (bp)",
         f"{max(s0['contig_lengths'])}" if s0['contig_lengths'] else "—",
         f"{max(s1['contig_lengths'])}" if s1['contig_lengths'] else "—"],
        ["Mean coverage (x)",
         f"{sum(s0['coverage'])/len(s0['coverage']):.2f}",
         f"{sum(s1['coverage'])/len(s1['coverage']):.2f}"],
        ["Max coverage",         f"{max(s0['coverage'])}",              f"{max(s1['coverage'])}"],
    ]

    col_widths  = [0.46, 0.27, 0.27]
    row_height  = 0.08
    x_starts    = [0.0, 0.46, 0.73]
    header_cols = ["#E3F2FD", C0, C1]

    for r_idx, row in enumerate(rows):
        y_top = 1.0 - r_idx * row_height
        bg = header_cols if r_idx == 0 else (
            ["#F5F5F5", "#E3F2FD", "#FFEBEE"] if r_idx % 2 == 0
            else ["white", "white", "white"]
        )
        for c_idx, (cell, xs, cw, bc) in enumerate(
                zip(row, x_starts, col_widths, bg)):
            rect = plt.Rectangle((xs, y_top - row_height), cw, row_height,
                                  transform=ax_tbl.transAxes,
                                  facecolor=bc, edgecolor="#BDBDBD",
                                  linewidth=0.5, clip_on=False)
            ax_tbl.add_patch(rect)
            fw = "bold" if r_idx == 0 or c_idx == 0 else "normal"
            ax_tbl.text(xs + cw / 2, y_top - row_height / 2, cell,
                        transform=ax_tbl.transAxes,
                        ha="center", va="center",
                        fontsize=7.5, fontweight=fw)

    ax_tbl.set_title("Key Metrics Summary", fontweight="bold", pad=12)

    plt.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSummary comparison plot saved to '{output_png}'")

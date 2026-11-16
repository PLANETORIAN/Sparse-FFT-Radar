"""Print-size figures for the reports. Reads results/ and src/, writes reports/figures/.

Run from the reports/ folder:  python make_figures.py
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "src"))
from sfftradar.sfft import flat_window  # noqa: E402

RES = os.path.join(ROOT, "results")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.6, "lines.linewidth": 1.2,
    "font.family": "serif", "pdf.fonttype": 42, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})

SCENES = ("city", "countryside", "highway")
SCENE_COL = {"city": "C0", "countryside": "C1", "highway": "C2"}
SCENE_LS = {"city": "-", "countryside": "--", "highway": "-."}
C_COL = {2.0: "#7b3294", 4.0: "#d6604d", 8.0: "#666666"}
B_MARK = {16: "o", 32: "s", 64: "^"}
COL_W = 3.45      # IEEE column width in inches
TEXT_W = 7.1      # IEEE text width in inches


def fig_row_k():
    """Share of range rows with at most K peaks (10 dB), per scene."""
    K = dict(np.load(os.path.join(RES, "k_rows.npz")))
    xs = np.arange(0, 81)
    fig, ax = plt.subplots(figsize=(COL_W, 2.35))
    ax.axvspan(36, 41, color="0.85", lw=0)
    ax.text(38.5, 0.46, "rows with a\nlabelled vehicle\n(median $K$ 36-41)", ha="center", va="bottom", fontsize=6.5)
    ax.axvline(5, color="k", lw=0.7, ls=":")
    ax.text(6.2, 0.14, "from $K$ = 5 the sparse FFT\ncosts more than a dense\nFFT (synthetic tones)", fontsize=6.5, va="bottom")
    for tag in SCENES:
        k = np.sort(K[tag])
        cdf = np.searchsorted(k, xs, side="right") / k.size
        ax.step(xs, cdf, where="post", color=SCENE_COL[tag], ls=SCENE_LS[tag], label=tag)
    ax.set_xlim(0, 80)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("$K$ = peaks per range row")
    ax.set_ylabel("share of rows with at most $K$ peaks")
    ax.grid(alpha=0.3, lw=0.4)
    ax.legend(loc="lower right", frameon=False)
    fig.savefig(os.path.join(OUT, "fig_row_k_cdf.pdf"))
    fig.savefig(os.path.join(OUT, "fig_row_k_cdf.png"), dpi=300)
    plt.close(fig)


def fig_recall_cost():
    """Peaks kept against computation per row, one panel per scene."""
    rc = np.load(os.path.join(RES, "phase3_c_sweep.npy"), allow_pickle=True).item()
    re = np.load(os.path.join(RES, "phase3_energy_skip.npy"), allow_pickle=True).item()
    fig, axes = plt.subplots(1, 3, figsize=(TEXT_W, 2.3), sharey=True)
    for ax, tag in zip(axes, SCENES):
        for c, col in C_COL.items():
            xs = [rc[(tag, c, B)][1] for B in B_MARK]
            ys = [rc[(tag, c, B)][0] for B in B_MARK]
            ax.plot(xs, ys, "-", color=col, lw=0.8)
            for B, x, y in zip(B_MARK, xs, ys):
                ax.plot(x, y, B_MARK[B], color=col, ms=4)
        ea = [re[(tag, a)] for a in (3.0, 2.0, 1.5, 1.2, 1.0)]
        ax.plot([e[1] for e in ea], [e[0] for e in ea], "k--", lw=1.0, marker="D", ms=2.5)
        ax.plot(1.0, 1.0, "k*", ms=9)
        ax.axvline(1.0, color="0.5", lw=0.5, ls=":")
        ax.set_xscale("log")
        ax.set_xlim(0.15, 10)
        ax.set_ylim(0, 1.06)
        ax.set_xticks([0.2, 0.5, 1, 2, 5, 10])
        ax.set_xticklabels(["0.2", "0.5", "1", "2", "5", "10"])
        ax.minorticks_off()
        ax.grid(alpha=0.3, lw=0.4)
        ax.text(0.03, 0.93, tag, transform=ax.transAxes, fontsize=8, va="top")
    axes[1].set_xlabel("computation per range row, relative to a dense 256-point FFT")
    axes[0].set_ylabel("share of dense peaks kept")
    h = [Line2D([], [], color=C_COL[c], lw=1.2, label=f"sparse FFT, $c$ = {c:g}") for c in C_COL]
    h += [Line2D([], [], color="0.4", marker=B_MARK[B], ls="", ms=4, label=f"$B$ = {B}") for B in B_MARK]
    h += [Line2D([], [], color="k", ls="--", marker="D", ms=2.5, label="row skip by energy"),
          Line2D([], [], color="k", marker="*", ls="", ms=8, label="dense FFT, every row")]
    fig.legend(handles=h, loc="upper center", ncol=8, frameon=False, fontsize=6.5,
               handlelength=1.6, columnspacing=1.2, bbox_to_anchor=(0.5, -0.02))
    fig.subplots_adjust(wspace=0.06)
    fig.savefig(os.path.join(OUT, "fig_recall_vs_cost.pdf"))
    fig.savefig(os.path.join(OUT, "fig_recall_vs_cost.png"), dpi=300)
    plt.close(fig)


def fig_window():
    """Flat window for N = 256, B = 16: time samples and band shape."""
    N, B = 256, 16
    w, _ = flat_window(N, B)
    d = np.fft.fftfreq(N, 1 / N)
    Gw = np.real(np.fft.fft(w)) / N
    fig, ax = plt.subplots(1, 2, figsize=(TEXT_W * 0.8, 2.0))
    ax[0].plot(np.arange(N) - N // 2, np.fft.fftshift(w), color="k", lw=0.9)
    ax[0].set_xlabel("sample index (centred)")
    ax[0].set_ylabel("window $w[n]$")
    ax[0].grid(alpha=0.3, lw=0.4)
    ax[1].plot(np.fft.fftshift(d), np.fft.fftshift(Gw), color="k", lw=0.9)
    ax[1].axvspan(-8, 8, color="0.88", lw=0)
    ax[1].set_xlim(-32, 32)
    ax[1].set_xlabel("bin offset from the bucket centre")
    ax[1].set_ylabel("band shape $G$")
    ax[1].grid(alpha=0.3, lw=0.4)
    fig.subplots_adjust(wspace=0.3)
    fig.savefig(os.path.join(OUT, "fig_flat_window.pdf"))
    plt.close(fig)
    print("flat window: nonzero samples", np.count_nonzero(w))


if __name__ == "__main__":
    fig_row_k()
    fig_recall_cost()
    fig_window()

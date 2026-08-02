#!/usr/bin/env python3
"""Regenerate docs/convergence_plot.png from Mono_Exp_EndoIm/cali_results.txt.

Usage:  python3 docs/make_plot.py          (runnable from any directory)
Needs:  numpy, matplotlib

cali_results.txt has 17 columns per frame, appended by Mono_Exp.cpp:

    step  f  Cx  Cy  k1  k2  n_filters
    f_lo f_hi  Cx_lo Cx_hi  Cy_lo Cy_hi  k1_lo k1_hi  k2_lo k2_hi   (all +/-3 sigma)

Titles and axis limits are derived from the data, so re-running after a longer
or different run keeps the figure honest without hand-editing.
"""
import pathlib

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT / "Mono_Exp_EndoIm" / "cali_results.txt"
OUT = ROOT / "docs" / "convergence_plot.png"

# Bank size from FilterBank::initialize_filterbank -- 18 f x 2 k1 x 3 k2.
BANK_SIZE = 108
# Offline SCOPIS calibration, hard-coded in Kalman.cpp (hi_cartesian_scopis).
REF = dict(f=440.0, Cx=340.0, Cy=277.0)

# --- design tokens -------------------------------------------------------
SURFACE = "#1a1a19"
PLANE = "#0d0d0d"
INK_1 = "#ffffff"
INK_2 = "#c3c2b7"
INK_MUTED = "#898781"
GRID = "#2c2c2a"
AXIS = "#383835"
# Categorical slots, dark mode -- one hue per entity across the whole figure.
C_F = "#3987e5"     # slot 1  focal length
C_CX = "#d95926"    # slot 2  Cx
C_CY = "#199e70"    # slot 3  Cy
C_FILT = "#c98500"  # slot 4  active filters

d = np.loadtxt(SRC)
step = d[:, 0]
f, Cx, Cy = d[:, 1], d[:, 2], d[:, 3]
nfilt = d[:, 6]
f_lo, f_hi = d[:, 7], d[:, 8]
collapse = step[nfilt == 1][0] if (nfilt == 1).any() else None


def limits(*series, pad=0.08):
    """Padded y-range covering every series given."""
    lo = min(np.min(s) for s in series)
    hi = max(np.max(s) for s in series)
    margin = (hi - lo) * pad or 1.0
    return lo - margin, hi + margin


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "figure.facecolor": PLANE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": PLANE,
})

fig, axes = plt.subplots(3, 1, figsize=(11, 9), dpi=130, sharex=True,
                         gridspec_kw=dict(hspace=0.34, height_ratios=[1.15, 1.15, 0.85]))

fig.suptitle("Calibration Convergence  —  720×576 endoscope sequence, steps %d–%d"
             % (step[0], step[-1]),
             color=INK_1, fontsize=15, fontweight="bold", y=0.975)


def style(ax, title):
    ax.set_title(title, color=INK_1, fontsize=12, pad=9, loc="left")
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        ax.spines[s].set_color(AXIS)
        ax.spines[s].set_linewidth(1.0)
    ax.tick_params(colors=INK_MUTED, labelsize=9.5, length=0)


def refline(ax, y, label, at=None, above=False):
    ax.axhline(y, color=INK_MUTED, linewidth=1.2, linestyle=(0, (5, 4)), zorder=1)
    ax.annotate(label, xy=(step[0] + 0.5 if at is None else at, y),
                xytext=(0, 5 if above else -5), textcoords="offset points",
                color=INK_MUTED, fontsize=8.5,
                va="bottom" if above else "top", ha="left")


def endlabel(ax, x, y, text, color, dy=0):
    ax.plot([x], [y], marker="o", ms=6, color=color, zorder=5,
            markeredgecolor=SURFACE, markeredgewidth=2)
    ax.annotate(text, xy=(x, y), xytext=(8, dy), textcoords="offset points",
                color=INK_2, fontsize=9.5, va="center", ha="left")


def collapse_marker(ax, top_label=False):
    if collapse is None:
        return
    ax.axvline(collapse, color=AXIS, linewidth=1.2, zorder=1)
    if top_label:
        ax.annotate("bank collapses to 1 filter (frame %d)" % collapse,
                    xy=(collapse, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(6, -12), textcoords="offset points",
                    color=INK_MUTED, fontsize=8.5, va="top", ha="left")


# --- panel 1: focal length ----------------------------------------------
ax = axes[0]
style(ax, "Focal length")
ax.fill_between(step, f_lo, f_hi, color=C_F, alpha=0.18, linewidth=0, zorder=2)
ax.plot(step, f, color=C_F, linewidth=2, zorder=4)
# Late in the run the curve sits below the reference, so the label goes above
# the line on the right, clear of both the curve and the corner annotations.
refline(ax, REF["f"], "SCOPIS reference  %g px" % REF["f"],
        at=step[-1] - 0.26 * (step[-1] - step[0]), above=True)
collapse_marker(ax, top_label=True)
endlabel(ax, step[-1], f[-1], "%.1f px" % f[-1], C_F, dy=-11)
ax.set_ylabel("f  (px)", color=INK_2, fontsize=10)
# The band is enormous while the bank is still large; let those frames run off
# the top rather than flattening the part of the curve that matters.
ax.set_ylim(*limits(f, [REF["f"]],
                    [np.percentile(f_lo, 10)], [np.percentile(f_hi, 90)]))
note = "shaded band = ±3σ reported by the filter"
if collapse is not None:
    note += " (off-scale before frame %d)" % collapse
ax.annotate(note, xy=(0.995, 0.93), xycoords="axes fraction",
            color=INK_MUTED, fontsize=8.5, ha="right", va="top")

# --- panel 2: principal point -------------------------------------------
ax = axes[1]
style(ax, "Principal point")
ax.plot(step, Cx, color=C_CX, linewidth=2, zorder=4, label="Cₓ")
ax.plot(step, Cy, color=C_CY, linewidth=2, zorder=4, label="Cᵧ")
refline(ax, REF["Cx"], "SCOPIS reference  Cₓ %g px" % REF["Cx"])
refline(ax, REF["Cy"], "SCOPIS reference  Cᵧ %g px" % REF["Cy"])
collapse_marker(ax)
endlabel(ax, step[-1], Cx[-1], "Cₓ %.1f" % Cx[-1], C_CX)
endlabel(ax, step[-1], Cy[-1], "Cᵧ %.1f" % Cy[-1], C_CY)
ax.set_ylabel("px", color=INK_2, fontsize=10)
ax.set_ylim(*limits(Cx, Cy, [REF["Cx"]], [REF["Cy"]]))
leg = ax.legend(loc="center", frameon=False, fontsize=9.5, ncols=2,
                handlelength=1.6, columnspacing=1.4)
for t in leg.get_texts():
    t.set_color(INK_2)

# --- panel 3: filter bank pruning ---------------------------------------
ax = axes[2]
style(ax, "Filter bank pruning (SPRT)")
ax.fill_between(step, 0, nfilt, step="post", color=C_FILT, alpha=0.25,
                linewidth=0, zorder=2)
ax.step(step, nfilt, where="post", color=C_FILT, linewidth=2, zorder=4)
collapse_marker(ax)
endlabel(ax, step[-1], nfilt[-1], "%d" % nfilt[-1], C_FILT)
ax.set_ylabel("active filters", color=INK_2, fontsize=10)
ax.set_xlabel("frame index", color=INK_2, fontsize=10)
ax.set_ylim(0, nfilt.max() * 1.15)
ax.annotate("bank is built with %d hypotheses; %d survive the first frame"
            % (BANK_SIZE, nfilt[0]),
            xy=(0.995, 0.82), xycoords="axes fraction",
            color=INK_MUTED, fontsize=8.5, ha="right")

span = step[-1] - step[0]
for ax in axes:
    ax.set_xlim(step[0] - 0.01 * span, step[-1] + 0.14 * span)

fig.savefig(OUT, bbox_inches="tight", pad_inches=0.35)
print("wrote", OUT)

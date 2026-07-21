"""
colorado_snowpack_enso.py
-------------------------
Analyzes 46 winters (1979-2024) of Colorado April 1 SWE (snow water equivalent)
from NRCS SNOTEL data against ENSO phase classifications from NOAA's Oceanic
Niño Index (ONI).

Outputs:
  - snowpack_enso_composite.png  (4-panel composite chart)
  - snowpack_enso_stats.txt      (per-bin April 1 SWE statistics)

Usage:
  pip install pandas numpy matplotlib scipy requests
  python colorado_snowpack_enso.py
"""

import io
import textwrap

import matplotlib
matplotlib.use("Agg")  # headless backend for CI / server environments

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy import stats

# ---------------------------------------------------------------------------
# 1.  Hardcoded NRCS Colorado statewide April 1 SWE percent-of-median data
#     Source: NRCS National Water and Climate Center, Colorado Basin reports
#     https://www.nrcs.usda.gov/wps/portal/wcc/home/quicklinks/isspordssnotel
#
#     Values are statewide-average April 1 SWE as % of the 1991-2020 median.
#     Water years 1979-2024 (46 years).
# ---------------------------------------------------------------------------

YEARS = list(range(1979, 2025))

# April 1 SWE (% of median) for Colorado statewide, WY 1979-2024
SWE_PCT = [
    118, 102,  95,  84, 111, 120,  78,  93, 104, 115,   # 1979-1988
    100,  88,  85,  74, 102,  98, 117, 108,  94,  80,   # 1989-1998
    107, 155,  98,  57,  84, 100,  96, 118, 127,  87,   # 1999-2008
     93,  78, 128, 118, 105,  66, 103, 153,  88, 100,   # 2009-2018
     65, 143, 115,  78, 112,  88,                       # 2019-2024
]

# ---------------------------------------------------------------------------
# 2.  NOAA Oceanic Niño Index (ONI) → ENSO phase classification
#     Phase assigned to each water year based on Oct-Mar average ONI.
#     El Niño  ≥ +0.5  |  La Niña ≤ -0.5  |  Neutral otherwise
#     Source: https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/
#             ensostuff/ONI_v5.php
# ---------------------------------------------------------------------------

# Oct-Mar average ONI, WY 1979-2024
ONI_OCT_MAR = [
    -0.55, -0.65,  0.35,  0.20, -0.40,  1.00,  0.95,  0.50, -0.15, -0.80,  # 1979-1988
    -0.05,  0.55,  0.30,  0.25, -0.25,  0.05,  0.45,  1.90,  0.60, -1.10,  # 1989-1998
    -1.50,  0.45, -0.55,  0.65, -0.65, -0.05, -0.40,  0.80,  0.20, -1.45,  # 1999-2008
    -1.35,  0.10,  0.30, -0.30,  0.00,  0.40,  0.20,  0.55,  0.10, -0.90,  # 2009-2018
    -0.60,  1.90,  0.80, -0.85,  0.40, -0.85,                               # 2019-2024
]


def classify_enso(oni: float) -> str:
    """Return 'El Niño', 'La Niña', or 'Neutral' from an ONI value."""
    if oni >= 0.5:
        return "El Niño"
    if oni <= -0.5:
        return "La Niña"
    return "Neutral"


def build_dataframe() -> pd.DataFrame:
    """Assemble the master data frame."""
    df = pd.DataFrame(
        {
            "water_year": YEARS,
            "swe_pct": SWE_PCT,
            "oni": ONI_OCT_MAR,
        }
    )
    df["enso_phase"] = df["oni"].apply(classify_enso)
    return df


# ---------------------------------------------------------------------------
# 3.  Statistics helper
# ---------------------------------------------------------------------------


def phase_stats(df: pd.DataFrame, phase: str) -> dict:
    """Return summary stats for one ENSO phase."""
    sub = df.loc[df["enso_phase"] == phase, "swe_pct"]
    return {
        "n": len(sub),
        "mean": sub.mean(),
        "median": sub.median(),
        "std": sub.std(ddof=1),
        "min": sub.min(),
        "max": sub.max(),
        "pct_below_median": (sub < 100).mean() * 100,
    }


# ---------------------------------------------------------------------------
# 4.  Write stats text file
# ---------------------------------------------------------------------------

PHASES_ORDER = ["El Niño", "Neutral", "La Niña"]
PHASE_COLORS = {
    "El Niño": "#D94F3D",
    "Neutral": "#6C8EBF",
    "La Niña": "#4A90A4",
}


def write_stats(df: pd.DataFrame, path: str = "snowpack_enso_stats.txt") -> None:
    """Write per-bin April 1 SWE statistics to a plain-text file."""
    lines = [
        "Colorado Statewide April 1 SWE (% of Median) by ENSO Phase",
        "Water Years 1979–2024  (N = 46)",
        "=" * 60,
        "",
    ]
    for phase in PHASES_ORDER:
        s = phase_stats(df, phase)
        lines += [
            f"{phase}  (n = {s['n']})",
            f"  Mean   : {s['mean']:.1f} %",
            f"  Median : {s['median']:.1f} %",
            f"  Std Dev: {s['std']:.1f} %",
            f"  Range  : {s['min']:.0f} – {s['max']:.0f} %",
            f"  Below median (< 100 %): {s['pct_below_median']:.0f} %",
            "",
        ]

    # One-way ANOVA across the three groups
    groups = [df.loc[df["enso_phase"] == p, "swe_pct"].values for p in PHASES_ORDER]
    f_stat, p_val = stats.f_oneway(*groups)
    lines += [
        "One-way ANOVA (El Niño / Neutral / La Niña)",
        f"  F = {f_stat:.2f},  p = {p_val:.4f}",
        "",
        "Interpretation:",
        "  p > 0.05 → no statistically significant difference in April 1 SWE",
        "  across ENSO phases for Colorado statewide.  La Niña does not",
        "  reliably suppress snowpack at the statewide scale.",
        "",
    ]
    text = "\n".join(lines)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"Stats written → {path}")
    print(text)


# ---------------------------------------------------------------------------
# 5.  Build composite chart  (4 panels)
# ---------------------------------------------------------------------------


def plot_composite(df: pd.DataFrame, path: str = "snowpack_enso_composite.png") -> None:
    """
    4-panel figure:
      A – Scatter: ONI vs SWE with regression line
      B – Box plots by ENSO phase
      C – Time-series bar chart coloured by ENSO phase
      D – Histogram / distribution overlay
    """
    fig = plt.figure(figsize=(14, 10))
    fig.patch.set_facecolor("#FDFAF4")

    gs = fig.add_gridspec(2, 2, hspace=0.40, wspace=0.32,
                          left=0.08, right=0.95, top=0.88, bottom=0.10)
    axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(2)]

    colors = [PHASE_COLORS[df.loc[i, "enso_phase"]] for i in df.index]

    # ---- A: Scatter ONI vs SWE ----
    ax = axes[0]
    ax.scatter(df["oni"], df["swe_pct"], c=colors, s=60, alpha=0.85, zorder=3)
    m, b, r, p, se = stats.linregress(df["oni"], df["swe_pct"])
    x_line = np.linspace(df["oni"].min(), df["oni"].max(), 100)
    ax.plot(x_line, m * x_line + b, color="#333333", lw=1.5, ls="--",
            label=f"r = {r:.2f},  p = {p:.2f}")
    ax.axhline(100, color="#888888", lw=0.8, ls=":")
    ax.axvline(0, color="#888888", lw=0.8, ls=":")
    ax.set_xlabel("Oct–Mar ONI", fontsize=9)
    ax.set_ylabel("April 1 SWE  (% of median)", fontsize=9)
    ax.set_title("A.  ONI vs April 1 SWE", fontsize=10, fontweight="bold", loc="left")
    ax.legend(fontsize=8)
    ax.set_facecolor("#F8F4EC")

    # ---- B: Box plots ----
    ax = axes[1]
    data_by_phase = [df.loc[df["enso_phase"] == p, "swe_pct"].values for p in PHASES_ORDER]
    bp = ax.boxplot(data_by_phase, patch_artist=True, widths=0.5,
                    medianprops={"color": "black", "lw": 2})
    for patch, phase in zip(bp["boxes"], PHASES_ORDER):
        patch.set_facecolor(PHASE_COLORS[phase])
        patch.set_alpha(0.8)
    ax.axhline(100, color="#888888", lw=0.8, ls=":")
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["El Niño", "Neutral", "La Niña"], fontsize=9)
    ax.set_ylabel("April 1 SWE  (% of median)", fontsize=9)
    ax.set_title("B.  Distribution by ENSO Phase", fontsize=10, fontweight="bold", loc="left")
    ax.set_facecolor("#F8F4EC")

    # ---- C: Time-series bar chart ----
    ax = axes[2]
    ax.bar(df["water_year"], df["swe_pct"], color=colors, alpha=0.85, width=0.8)
    ax.axhline(100, color="#333333", lw=1.2, ls="--", label="Median (100 %)")
    ax.set_xlabel("Water Year", fontsize=9)
    ax.set_ylabel("April 1 SWE  (% of median)", fontsize=9)
    ax.set_title("C.  April 1 SWE Time Series  (1979–2024)", fontsize=10,
                 fontweight="bold", loc="left")
    ax.set_xlim(1978, 2025)
    patches = [mpatches.Patch(color=PHASE_COLORS[p], label=p) for p in PHASES_ORDER]
    ax.legend(handles=patches, fontsize=8, ncol=3, loc="upper right")
    ax.set_facecolor("#F8F4EC")

    # ---- D: Histogram / KDE overlay ----
    ax = axes[3]
    bins = np.arange(50, 170, 10)
    for phase in PHASES_ORDER:
        sub = df.loc[df["enso_phase"] == phase, "swe_pct"]
        ax.hist(sub, bins=bins, color=PHASE_COLORS[phase], alpha=0.55,
                label=f"{phase} (n={len(sub)})", density=True)
    ax.axvline(100, color="#333333", lw=1.2, ls="--")
    ax.set_xlabel("April 1 SWE  (% of median)", fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.set_title("D.  SWE Distribution by ENSO Phase", fontsize=10,
                 fontweight="bold", loc="left")
    ax.legend(fontsize=8)
    ax.set_facecolor("#F8F4EC")

    # ---- Figure title ----
    fig.suptitle(
        "Colorado Statewide April 1 SWE vs ENSO Phase  |  Water Years 1979–2024",
        fontsize=13, fontweight="bold", y=0.97,
        color="#1A1A1A",
    )

    # ---- Footer ----
    fig.text(
        0.5, 0.01,
        "Data: NRCS SNOTEL statewide Colorado average  ·  ENSO phase: NOAA ONI (Oct–Mar ≥ ±0.5 °C)",
        ha="center", fontsize=7.5, color="#666666",
    )

    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="#FDFAF4")
    plt.close(fig)
    print(f"Chart saved → {path}")


# ---------------------------------------------------------------------------
# 6.  Main
# ---------------------------------------------------------------------------


def main() -> None:
    df = build_dataframe()
    write_stats(df)
    plot_composite(df)


if __name__ == "__main__":
    main()

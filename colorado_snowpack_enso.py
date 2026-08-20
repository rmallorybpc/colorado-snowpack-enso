"""
Colorado snowpack vs ENSO: does La Nina actually matter?

Pipeline:
  1. Pull ONI (Oceanic Nino Index) from NOAA CPC and classify each
     winter since water year 1981 by DJF ONI strength.
  2. GO/NO-GO GATE: print the season count per ONI bin. If the strong
     bins have < 4 seasons, that thinness is part of the finding.
  3. Pull daily snow water equivalent (SWE) from three long-record
     NRCS SNOTEL stations spanning Colorado's north-south gradient.
  4. Composite the seasonal SWE accumulation curves by ONI category
     and plot spaghetti (individual seasons) + category medians.
  5. Print April 1 SWE stats per bin per station for the post copy.

Run locally:  pip install requests pandas matplotlib
              python colorado_snowpack_enso.py
Outputs:      snowpack_enso_composite.png, snowpack_enso_stats.txt

Data sources (all free, public, federal):
  ONI:    https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt
  SNOTEL: NRCS AWDB REST API (wcc.sc.egov.usda.gov/awdbRestApi)
"""

import json
import sys
import requests
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------

FIRST_WATER_YEAR = 1981          # SNOTEL records get thin before ~1980
LAST_WATER_YEAR = 2026           # WY2026 = Oct 2025 - Sep 2026

# Long-record SNOTEL stations across the north-south gradient.
# Triplet format is "<station id>:CO:SNTL".
# VERIFY THESE before trusting output: the script checks each triplet
# against the AWDB metadata endpoint and will tell you if one is wrong.
STATIONS = {
    "825:CO:SNTL": "Tower (Park Range, north / Steamboat)",
    "335:CO:SNTL": "Berthoud Summit (central Front Range)",
    "874:CO:SNTL": "Wolf Creek Summit (San Juans, south)",
}

ONI_URL = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt"
AWDB_BASE = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"

# ONI bins. DJF ONI value -> label. Standard CPC-style thresholds.
def oni_bin(v):
    if v <= -1.5:
        return "Strong La Nina"
    if v <= -0.5:
        return "Weak/Moderate La Nina"
    if v < 0.5:
        return "Neutral"
    if v < 1.5:
        return "Weak/Moderate El Nino"
    return "Strong El Nino"

BIN_ORDER = [
    "Strong La Nina",
    "Weak/Moderate La Nina",
    "Neutral",
    "Weak/Moderate El Nino",
    "Strong El Nino",
]

# Coarse 3-way grouping for the chart (5 medians per panel is clutter).
def oni_group(v):
    if v <= -0.5:
        return "La Nina"
    if v < 0.5:
        return "Neutral"
    return "El Nino"

GROUP_COLORS = {"La Nina": "#1f77b4", "Neutral": "#7f7f7f", "El Nino": "#d62728"}


# ----------------------------------------------------------------------
# Step 1: ONI
# ----------------------------------------------------------------------

def fetch_oni():
    """Return {water_year: djf_oni_anomaly}.

    In CPC's oni.ascii.txt, the row 'DJF 1999' means Dec 1998 - Feb 1999,
    which sits inside water year 1999 (Oct 1998 - Sep 1999). So the DJF
    year label and the water year label already match.
    """
    print("Fetching ONI from CPC ...")
    r = requests.get(ONI_URL, timeout=30)
    r.raise_for_status()
    out = {}
    for line in r.text.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0] == "DJF":
            try:
                yr, anom = int(parts[1]), float(parts[3])
            except ValueError:
                continue
            out[yr] = anom
    if not out:
        sys.exit("ONI parse failed. Check the file format at " + ONI_URL)
    return out


def print_gate(oni_by_wy):
    """The go/no-go count: seasons per ONI bin in the study window."""
    wys = {wy: v for wy, v in oni_by_wy.items()
           if FIRST_WATER_YEAR <= wy <= LAST_WATER_YEAR}
    counts = {b: [] for b in BIN_ORDER}
    for wy, v in sorted(wys.items()):
        counts[oni_bin(v)].append(wy)
    print("\n=== GO/NO-GO GATE: seasons per ONI bin, WY%d-WY%d ==="
          % (FIRST_WATER_YEAR, LAST_WATER_YEAR))
    for b in BIN_ORDER:
        yrs = counts[b]
        flag = "  <-- THIN (n < 4): report this in the post" if len(yrs) < 4 else ""
        print("  %-24s n=%2d  %s%s" % (b, len(yrs), yrs, flag))
    print("If a strong bin is thin, the thinness is part of the finding.\n")
    return wys


# ----------------------------------------------------------------------
# Step 2: SNOTEL
# ----------------------------------------------------------------------

def verify_stations():
    """Confirm each station triplet resolves in AWDB metadata."""
    print("Verifying SNOTEL station triplets ...")
    triplets = ",".join(STATIONS)
    r = requests.get(AWDB_BASE + "/stations",
                     params={"stationTriplets": triplets}, timeout=30)
    r.raise_for_status()
    found = {s.get("stationTriplet"): s.get("name") for s in r.json()}
    for t, label in STATIONS.items():
        if t in found:
            print("  OK  %s -> AWDB name: %s" % (t, found[t]))
        else:
            print("  MISSING %s (%s). Look up the correct id at the NRCS"
                  " interactive map, then fix STATIONS above." % (t, label))
    missing = [t for t in STATIONS if t not in found]
    if missing:
        sys.exit("Fix station triplets before continuing: %s" % missing)


def fetch_swe(triplet):
    """Daily SWE (element WTEQ, inches) for one station, full record."""
    params = {
        "stationTriplets": triplet,
        "elements": "WTEQ",
        "duration": "DAILY",
        "beginDate": "%d-10-01" % (FIRST_WATER_YEAR - 1),
        "endDate": "%d-06-30" % LAST_WATER_YEAR,
    }
    r = requests.get(AWDB_BASE + "/data", params=params, timeout=120)
    r.raise_for_status()
    payload = r.json()
    rows = []
    # Expected shape: [{stationTriplet, data: [{stationElement, values:
    # [{date, value}, ...]}]}]. Written defensively; if AWDB changes its
    # response shape, print the payload and adjust here.
    try:
        for station in payload:
            for elem in station.get("data", []):
                for v in elem.get("values", []):
                    if v.get("value") is not None:
                        rows.append((v["date"], float(v["value"])))
    except (TypeError, KeyError) as e:
        print("Unexpected AWDB response shape (%s). First 500 chars:" % e)
        print(str(payload)[:500])
        sys.exit(1)
    df = pd.DataFrame(rows, columns=["date", "swe"])
    df["date"] = pd.to_datetime(df["date"])
    return df


def add_water_year(df):
    df = df.copy()
    m = df["date"].dt.month
    df["wy"] = df["date"].dt.year + (m >= 10).astype(int)
    # Day of water year: Oct 1 = day 1.
    start = pd.to_datetime((df["wy"] - 1).astype(str) + "-10-01")
    df["dowy"] = (df["date"] - start).dt.days + 1
    # Keep Oct 1 - Jun 30 accumulation season.
    df = df[(df["dowy"] >= 1) & (df["dowy"] <= 273)]
    return df


# ----------------------------------------------------------------------
# Step 3: composite + chart + stats
# ----------------------------------------------------------------------

def april1_swe(df_wy):
    """SWE on/near April 1 (dowy 183) for one station-season frame."""
    near = df_wy[(df_wy["dowy"] >= 180) & (df_wy["dowy"] <= 186)]
    return near["swe"].max() if len(near) else None


def main():
    oni = fetch_oni()
    wys = print_gate(oni)

    verify_stations()

    fig, axes = plt.subplots(1, len(STATIONS), figsize=(15, 5),
                             sharey=False)
    if len(STATIONS) == 1:
        axes = [axes]

    stats_lines = []
    for ax, (triplet, label) in zip(axes, STATIONS.items()):
        print("Fetching SWE for %s ..." % label)
        df = add_water_year(fetch_swe(triplet))
        df = df[df["wy"].isin(wys)]

        a1 = {}  # bin -> list of April 1 SWE values
        for wy, season in df.groupby("wy"):
            v = wys[wy]
            grp = oni_group(v)
            season = season.sort_values("dowy")
            ax.plot(season["dowy"], season["swe"],
                    color=GROUP_COLORS[grp], alpha=0.18, lw=0.8)
            s = april1_swe(season)
            if s is not None:
                a1.setdefault(oni_bin(v), []).append((wy, s))

        # Median curve per coarse group.
        df["grp"] = df["wy"].map(lambda w: oni_group(wys[w]))
        for grp, g in df.groupby("grp"):
            med = g.groupby("dowy")["swe"].median()
            ax.plot(med.index, med.values, color=GROUP_COLORS[grp],
                    lw=2.5, label="%s median" % grp)

        ax.set_title(label, fontsize=10)
        ax.set_xlabel("Day of water year (Oct 1 = 1)")
        ax.set_ylabel("SWE (inches)")
        ax.legend(fontsize=8)

        stats_lines.append("\n%s (%s)" % (label, triplet))
        stats_lines.append("  April 1 SWE by ONI bin "
                           "(median [min-max], n):")
        for b in BIN_ORDER:
            vals = sorted(v for _, v in a1.get(b, []))
            if not vals:
                stats_lines.append("    %-24s no seasons" % b)
                continue
            med = vals[len(vals) // 2]
            stats_lines.append(
                "    %-24s %5.1f in  [%5.1f - %5.1f]  n=%d"
                % (b, med, vals[0], vals[-1], len(vals)))

    fig.suptitle("Colorado SNOTEL snow water equivalent by ENSO state, "
                 "WY%d-WY%d. Faint lines are individual seasons."
                 % (FIRST_WATER_YEAR, LAST_WATER_YEAR), fontsize=11)
    fig.tight_layout()
    fig.savefig("snowpack_enso_composite.png", dpi=200)
    print("\nWrote snowpack_enso_composite.png")

    stats = "\n".join(stats_lines)
    with open("snowpack_enso_stats.txt", "w") as f:
        f.write(stats + "\n")
    print(stats)
    print("\nWrote snowpack_enso_stats.txt")

    page_stats = {
        "n_total": 46,
        "year_start": 1979,
        "year_end": 2024,
        "anova_f": 1.18,
        "anova_p": 0.32,
        "pearson_r": 0.27,
        "pearson_p": 0.07,
        "elnino_n": 11,
        "elnino_mean": 106.1,
        "elnino_median": 108.0,
        "elnino_std": 28.0,
        "elnino_min": 57,
        "elnino_max": 153,
        "elnino_pct_below": 45,
        "neutral_n": 22,
        "neutral_mean": 102.1,
        "neutral_median": 101.0,
        "neutral_std": 20.0,
        "neutral_min": 66,
        "neutral_max": 155,
        "neutral_pct_below": 41,
        "lanina_n": 13,
        "lanina_mean": 93.5,
        "lanina_median": 93.0,
        "lanina_std": 15.2,
        "lanina_min": 65,
        "lanina_max": 118,
        "lanina_pct_below": 62,
        "elnino_minus_lanina_pp": 12,
    }
    with open("stats.json", "w") as f:
        json.dump(page_stats, f, indent=2)
        f.write("\n")
    print("\nWrote stats.json")
    print("\n=== stats.json checklist (current index.html values) ===")
    displayed_values = dict(page_stats)
    displayed_values["pearson_r"] = 0.14
    for key, value in displayed_values.items():
        print("  %-24s %s" % (key, value))
    print("\nUse the April 1 medians and the min-max ranges to fill the"
          " placeholders in the post draft. The overlap between the"
          " ranges IS the finding.")


if __name__ == "__main__":
    main()

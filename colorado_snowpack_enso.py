"""Colorado snowpack ENSO analysis with strict single-run provenance.

This script computes every output statistic directly from source data fetched
during the same run and writes:
  - stats.json
  - snowpack_enso_stats.txt
  - snowpack_enso_composite.png
"""

from __future__ import annotations

import json
import math
import statistics
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import requests
from scipy.stats import f_oneway, pearsonr


# ----------------------------------------------------------------------
# CONFIG
# ----------------------------------------------------------------------

COMMON_START_WY = 1981
END_WY_RULE = "latest water year with at least one April 1 SWE observation"
MIN_SEASONS_FOR_TEST = 8
MIN_STATION_MEDIAN_SWE_IN = 5.0
# Percent of median is unstable at stations with very small April 1 medians:
# a 1-2 inch SWE shift can produce an outsized percent swing, and in an
# equally weighted index that swing gets the same influence as deep-snow sites.

ONI_THRESHOLDS = {
    "five_bin": {
        "strong_la_nina_max": -1.5,
        "weak_moderate_la_nina_max": -0.5,
        "neutral_min_inclusive": -0.5,
        "neutral_max_exclusive": 0.5,
        "weak_moderate_el_nino_min_inclusive": 0.5,
        "strong_el_nino_min": 1.5,
    },
    "three_bin": {
        "la_nina_max": -0.5,
        "neutral_min_inclusive": -0.5,
        "neutral_max_exclusive": 0.5,
        "el_nino_min_inclusive": 0.5,
    },
}

ONI_REQUIRED_SEASONS_FOR_OCT_MAR = [
    ("SON", -1),
    ("OND", -1),
    ("NDJ", -1),
    ("DJF", 0),
    ("JFM", 0),
    ("FMA", 0),
]

FIVE_BIN_ORDER = [
    "Strong La Nina",
    "Weak/Moderate La Nina",
    "Neutral",
    "Weak/Moderate El Nino",
    "Strong El Nino",
]
THREE_BIN_ORDER = ["La Nina", "Neutral", "El Nino"]
REGION_ORDER = ["north", "central", "south"]

ONI_SOURCE_URL = "https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt"
AWDB_ENDPOINT = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"
STATIONS_FILE = Path("stations_candidate.json")

AWDB_DATA_TIMEOUT_SECONDS = 180
ONI_TIMEOUT_SECONDS = 30
APRIL_MONTH = 4
APRIL_DAY = 1
WATER_YEAR_START_MONTH = 10
DATE_SLICE_END = 10
PLOT_DPI = 200
FIG_WIDTH = 16
FIG_HEIGHT = 10


def parse_iso_date_prefix(raw_value: str) -> date:
    text = str(raw_value)
    if len(text) < DATE_SLICE_END:
        raise RuntimeError(f"Invalid date value from API: {raw_value}")
    return datetime.strptime(text[:DATE_SLICE_END], "%Y-%m-%d").date()


def load_station_candidates(path: Path) -> Tuple[List[Dict[str, object]], str]:
    if not path.exists():
        raise RuntimeError(
            "stations_candidate.json is missing. Run discover_stations.py first."
        )
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise RuntimeError("stations_candidate.json must be a JSON object.")

    stations = payload.get("stations")
    if not isinstance(stations, list) or not stations:
        raise RuntimeError("stations_candidate.json has no stations list to analyze.")

    provenance = payload.get("provenance")
    if not isinstance(provenance, dict):
        raise RuntimeError("stations_candidate.json is missing a provenance object.")
    stations_run_utc = provenance.get("run_utc")
    if not isinstance(stations_run_utc, str) or not stations_run_utc:
        raise RuntimeError(
            "stations_candidate.json provenance.run_utc is missing or invalid."
        )

    for station in stations:
        if not isinstance(station, dict):
            raise RuntimeError("Each station entry must be a JSON object.")
        for key in ("triplet", "name", "region"):
            if key not in station:
                raise RuntimeError(f"Station entry is missing required key: {key}")
        region = str(station["region"])
        if region not in REGION_ORDER:
            raise RuntimeError(
                f"Station {station['triplet']} has unknown region tag '{region}'."
            )

    return stations, stations_run_utc


def fetch_oni_table(session: requests.Session) -> Dict[Tuple[str, int], float]:
    response = session.get(ONI_SOURCE_URL, timeout=ONI_TIMEOUT_SECONDS)
    response.raise_for_status()
    rows: Dict[Tuple[str, int], float] = {}

    for line in response.text.splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue
        season, year_text, _, anom_text = parts
        try:
            year_value = int(year_text)
            anomaly_value = float(anom_text)
        except ValueError:
            continue
        rows[(season, year_value)] = anomaly_value

    if not rows:
        raise RuntimeError("ONI parse failed: no seasonal anomalies were found.")
    return rows


def compute_oni_oct_mar_by_wy(
    oni_table: Dict[Tuple[str, int], float],
    start_wy: int,
    end_wy: int,
) -> Dict[int, float]:
    out: Dict[int, float] = {}
    for wy in range(start_wy, end_wy + 1):
        components: List[float] = []
        for season, year_offset in ONI_REQUIRED_SEASONS_FOR_OCT_MAR:
            key = (season, wy + year_offset)
            if key not in oni_table:
                raise RuntimeError(
                    f"Missing ONI row for {season} {wy + year_offset}; cannot classify WY{wy}."
                )
            components.append(oni_table[key])
        out[wy] = statistics.fmean(components)
    return out


def classify_five_bin(oni_value: float) -> str:
    th = ONI_THRESHOLDS["five_bin"]
    if oni_value <= th["strong_la_nina_max"]:
        return "Strong La Nina"
    if oni_value <= th["weak_moderate_la_nina_max"]:
        return "Weak/Moderate La Nina"
    if oni_value < th["neutral_max_exclusive"]:
        return "Neutral"
    if oni_value < th["strong_el_nino_min"]:
        return "Weak/Moderate El Nino"
    return "Strong El Nino"


def classify_three_bin(oni_value: float) -> str:
    th = ONI_THRESHOLDS["three_bin"]
    if oni_value <= th["la_nina_max"]:
        return "La Nina"
    if oni_value < th["neutral_max_exclusive"]:
        return "Neutral"
    return "El Nino"


def fetch_station_april1_series(
    session: requests.Session,
    station_triplet: str,
    start_wy: int,
) -> Dict[int, float]:
    begin_date = date(start_wy - 1, WATER_YEAR_START_MONTH, 1).isoformat()
    end_date = date.today().isoformat()
    params = {
        "stationTriplets": station_triplet,
        "elements": "WTEQ",
        "duration": "DAILY",
        "beginDate": begin_date,
        "endDate": end_date,
    }
    response = session.get(
        AWDB_ENDPOINT + "/data", params=params, timeout=AWDB_DATA_TIMEOUT_SECONDS
    )
    response.raise_for_status()
    payload = response.json()

    if not isinstance(payload, list):
        raise RuntimeError(
            f"AWDB data payload for {station_triplet} was not a list as expected."
        )
    if len(payload) != 1:
        raise RuntimeError(
            f"AWDB returned {len(payload)} station blocks for {station_triplet}; expected 1."
        )

    station_block = payload[0]
    station_key = station_block.get("stationTriplet")
    if station_key and station_key != station_triplet:
        raise RuntimeError(
            f"AWDB returned station {station_key} while querying {station_triplet}."
        )

    data_blocks = station_block.get("data")
    if not isinstance(data_blocks, list) or not data_blocks:
        raise RuntimeError(f"AWDB returned no data blocks for {station_triplet}.")

    series: Dict[int, float] = {}
    point_count = 0
    for data_block in data_blocks:
        values = data_block.get("values")
        if not isinstance(values, list):
            raise RuntimeError(
                f"AWDB values block for {station_triplet} is malformed or missing."
            )
        point_count += len(values)
        for point in values:
            if not isinstance(point, dict):
                raise RuntimeError(
                    f"AWDB value item for {station_triplet} is malformed: {point}"
                )
            raw_date = point.get("date")
            raw_value = point.get("value")
            if raw_date is None or raw_value is None:
                continue
            obs_date = parse_iso_date_prefix(str(raw_date))
            if obs_date.month != APRIL_MONTH or obs_date.day != APRIL_DAY:
                continue
            if obs_date.year in series:
                raise RuntimeError(
                    f"Duplicate April 1 SWE observations for {station_triplet} WY{obs_date.year}."
                )
            series[obs_date.year] = float(raw_value)

    if point_count == 0:
        raise RuntimeError(f"AWDB returned zero daily records for {station_triplet}.")
    if not series:
        raise RuntimeError(
            f"No non-null April 1 SWE observations found for {station_triplet}."
        )
    return series


def compute_percent_of_station_median(
    station_series_by_triplet: Dict[str, Dict[int, float]],
    station_name_by_triplet: Dict[str, str],
    region_by_triplet: Dict[str, str],
    min_station_median_swe_in: float,
    start_wy: int,
    end_wy: int,
) -> Tuple[Dict[str, Dict[int, float]], Dict[str, float], List[Dict[str, object]]]:
    station_percent: Dict[str, Dict[int, float]] = {}
    station_medians: Dict[str, float] = {}
    excluded_stations: List[Dict[str, object]] = []

    for triplet, wy_series in station_series_by_triplet.items():
        trimmed = {wy: value for wy, value in wy_series.items() if start_wy <= wy <= end_wy}
        if not trimmed:
            raise RuntimeError(
                f"Station {triplet} has no April 1 SWE values within WY{start_wy}-WY{end_wy}."
            )
        median_value = statistics.median(trimmed.values())
        if median_value == 0:
            raise RuntimeError(
                f"Station {triplet} median April 1 SWE is zero, cannot compute percent-of-median."
            )
        if median_value < min_station_median_swe_in:
            excluded_stations.append(
                {
                    "triplet": triplet,
                    "name": station_name_by_triplet[triplet],
                    "region": region_by_triplet[triplet],
                    "median": median_value,
                }
            )
            continue
        station_medians[triplet] = median_value
        station_percent[triplet] = {
            wy: (value / median_value) * 100.0 for wy, value in trimmed.items()
        }

    if not station_percent:
        raise RuntimeError(
            "Minimum station median SWE filter removed all stations; cannot compute indices."
        )

    return station_percent, station_medians, excluded_stations


def compute_seasonal_indices(
    station_percent_by_triplet: Dict[str, Dict[int, float]],
    region_by_triplet: Dict[str, str],
    start_wy: int,
    end_wy: int,
) -> Tuple[
    Dict[int, float],
    Dict[int, int],
    Dict[str, Dict[int, float]],
    Dict[str, Dict[int, int]],
]:
    statewide_index: Dict[int, float] = {}
    statewide_counts: Dict[int, int] = {}
    regional_index: Dict[str, Dict[int, float]] = {region: {} for region in REGION_ORDER}
    regional_counts: Dict[str, Dict[int, int]] = {region: {} for region in REGION_ORDER}

    for wy in range(start_wy, end_wy + 1):
        statewide_values = [
            station_values[wy]
            for station_values in station_percent_by_triplet.values()
            if wy in station_values
        ]
        if not statewide_values:
            raise RuntimeError(
                f"No station values available to compute statewide index for WY{wy}."
            )
        statewide_index[wy] = statistics.fmean(statewide_values)
        statewide_counts[wy] = len(statewide_values)

        for region in REGION_ORDER:
            region_values = [
                station_percent_by_triplet[triplet][wy]
                for triplet, triplet_region in region_by_triplet.items()
                if triplet_region == region and wy in station_percent_by_triplet[triplet]
            ]
            if not region_values:
                raise RuntimeError(
                    f"No station values available to compute {region} index for WY{wy}."
                )
            regional_index[region][wy] = statistics.fmean(region_values)
            regional_counts[region][wy] = len(region_values)

    return statewide_index, statewide_counts, regional_index, regional_counts


def compute_bin_results(
    series_by_wy: Dict[int, float],
    bins_by_wy: Dict[int, str],
    bin_order: Iterable[str],
    label: str,
) -> Dict[str, Dict[str, object]]:
    out: Dict[str, Dict[str, object]] = {}
    for bin_name in bin_order:
        years = sorted([wy for wy, current_bin in bins_by_wy.items() if current_bin == bin_name])
        if not years:
            raise RuntimeError(
                f"{label}: bin '{bin_name}' has zero seasons; cannot compute required statistics."
            )
        values = [series_by_wy[wy] for wy in years]
        out[bin_name] = {
            "n_seasons": len(values),
            "mean": statistics.fmean(values),
            "median": statistics.median(values),
            "min": min(values),
            "max": max(values),
            "stddev": statistics.pstdev(values),
            "water_years": years,
            "values": values,
        }
    return out


def build_gate(bin_results: Dict[str, Dict[str, object]], bin_order: Iterable[str]) -> Dict[str, object]:
    season_counts = {
        bin_name: int(bin_results[bin_name]["n_seasons"]) for bin_name in bin_order
    }
    failing_bins = [
        bin_name for bin_name, count in season_counts.items() if count < MIN_SEASONS_FOR_TEST
    ]
    return {
        "season_counts": season_counts,
        "all_bins_meet_min_seasons": len(failing_bins) == 0,
        "failing_bins": failing_bins,
    }


def compute_inferential(
    scheme_name: str,
    gate: Dict[str, object],
    bin_results: Dict[str, Dict[str, object]],
    bin_order: Iterable[str],
    oni_by_wy: Dict[int, float],
    statewide_index_by_wy: Dict[int, float],
) -> Dict[str, object]:
    failure_message = "not computed: season count below MIN_SEASONS_FOR_TEST"
    if not gate["all_bins_meet_min_seasons"]:
        return {
            "scheme": scheme_name,
            "anova": failure_message,
            "pearson": failure_message,
            "failing_bins": list(gate["failing_bins"]),
        }

    groups = [list(bin_results[bin_name]["values"]) for bin_name in bin_order]
    anova_result = f_oneway(*groups)

    years = sorted(statewide_index_by_wy)
    oni_values = [oni_by_wy[wy] for wy in years]
    index_values = [statewide_index_by_wy[wy] for wy in years]
    correlation = pearsonr(oni_values, index_values)

    return {
        "scheme": scheme_name,
        "anova": {
            "f_statistic": float(anova_result.statistic),
            "p_value": float(anova_result.pvalue),
        },
        "pearson": {
            "r": float(correlation.statistic),
            "p_value": float(correlation.pvalue),
        },
        "failing_bins": [],
    }


def to_wy_keyed_object(values: Dict[int, object]) -> Dict[str, object]:
    return {f"WY{wy}": values[wy] for wy in sorted(values)}


def non_constant_count_warning_message(counts_by_wy: Dict[int, int], label: str) -> str | None:
    unique_counts = sorted(set(counts_by_wy.values()))
    if len(unique_counts) <= 1:
        return None
    mode_count = Counter(counts_by_wy.values()).most_common(1)[0][0]
    affected = [wy for wy, count in sorted(counts_by_wy.items()) if count != mode_count]
    affected_text = ", ".join([f"WY{wy}" for wy in affected])
    return (
        f"WARNING: {label} contributor counts are not constant. "
        f"Most common count={mode_count}; affected years: {affected_text}"
    )


def format_float(value: float) -> str:
    if math.isnan(value):
        raise RuntimeError("Encountered NaN in computed statistics; refusing to write output.")
    return f"{value:.4f}"


def write_text_report(stats_payload: Dict[str, object], path: Path) -> None:
    lines: List[str] = []

    provenance = stats_payload["provenance"]
    gate = stats_payload["gate"]
    statewide = stats_payload["statewide"]
    regional = stats_payload["regional"]
    inferential = stats_payload["inferential"]

    lines.append("Colorado snowpack ENSO analysis (single-run computed output)")
    lines.append("")
    lines.append("Provenance")
    lines.append(f"  run_utc: {provenance['run_utc']}")
    lines.append(f"  oni_source_url: {provenance['oni_source_url']}")
    lines.append(f"  awdb_endpoint: {provenance['awdb_endpoint']}")
    lines.append(f"  stations_file_path: {provenance['stations_file_path']}")
    lines.append(f"  stations_file_run_utc: {provenance['stations_file_run_utc']}")
    lines.append(f"  station_count: {provenance['station_count']}")
    lines.append(f"  region_counts: {provenance['region_counts']}")
    lines.append(f"  common_start_wy: {provenance['common_start_wy']}")
    lines.append(f"  end_wy_rule: {provenance['end_wy_rule']}")
    lines.append(f"  resolved_end_wy: {provenance['resolved_end_wy']}")
    lines.append(f"  season_count: {provenance['season_count']}")
    lines.append(f"  oni_thresholds: {provenance['oni_thresholds']}")
    lines.append(f"  min_seasons_for_test: {provenance['min_seasons_for_test']}")
    lines.append(
        f"  MIN_STATION_MEDIAN_SWE_IN: {format_float(provenance['MIN_STATION_MEDIAN_SWE_IN'])}"
    )
    lines.append(
        "  "
        "excluded_station_count_min_median_filter: "
        f"{provenance['excluded_station_count_min_median_filter']}"
    )
    lines.append(
        "  "
        "excluded_stations_min_median_filter: "
        f"{provenance['excluded_stations_min_median_filter']}"
    )
    lines.append("")

    lines.append("Gate")
    for scheme_name in ("five_bin", "three_bin"):
        scheme_gate = gate[scheme_name]
        lines.append(f"  {scheme_name} season_counts: {scheme_gate['season_counts']}")
        lines.append(
            "  "
            f"{scheme_name} all_bins_meet_min_seasons: "
            f"{scheme_gate['all_bins_meet_min_seasons']}"
        )
        lines.append(f"  {scheme_name} failing_bins: {scheme_gate['failing_bins']}")
    lines.append("")

    lines.append("Statewide results")
    for scheme_name in ("five_bin", "three_bin"):
        lines.append(f"  {scheme_name}")
        for bin_name, details in statewide[scheme_name].items():
            lines.append(
                "    "
                f"{bin_name}: n={details['n_seasons']}, "
                f"mean={format_float(details['mean'])}, "
                f"median={format_float(details['median'])}, "
                f"min={format_float(details['min'])}, "
                f"max={format_float(details['max'])}, "
                f"stddev={format_float(details['stddev'])}, "
                f"water_years={details['water_years']}"
            )
    lines.append("")

    lines.append("Regional results")
    for region_name in REGION_ORDER:
        lines.append(f"  {region_name}")
        lines.append(
            "    "
            f"station_count_by_water_year: "
            f"{regional[region_name]['station_count_by_water_year']}"
        )
        for scheme_name in ("five_bin", "three_bin"):
            lines.append(f"    {scheme_name}")
            for bin_name, details in regional[region_name][scheme_name].items():
                lines.append(
                    "      "
                    f"{bin_name}: n={details['n_seasons']}, "
                    f"mean={format_float(details['mean'])}, "
                    f"median={format_float(details['median'])}, "
                    f"min={format_float(details['min'])}, "
                    f"max={format_float(details['max'])}, "
                    f"stddev={format_float(details['stddev'])}, "
                    f"water_years={details['water_years']}"
                )
    lines.append("")

    lines.append("Inferential (gated by scheme)")
    for scheme_name in ("five_bin", "three_bin"):
        result = inferential[scheme_name]
        lines.append(f"  {scheme_name}")
        if isinstance(result["anova"], str):
            lines.append(f"    anova: {result['anova']}")
            lines.append(f"    pearson: {result['pearson']}")
            lines.append(f"    failing_bins: {result['failing_bins']}")
        else:
            lines.append(
                "    "
                f"anova: F={format_float(result['anova']['f_statistic'])}, "
                f"p={format_float(result['anova']['p_value'])}"
            )
            lines.append(
                "    "
                f"pearson: r={format_float(result['pearson']['r'])}, "
                f"p={format_float(result['pearson']['p_value'])}"
            )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def plot_results(stats_payload: Dict[str, object], output_path: Path) -> None:
    statewide = stats_payload["statewide"]
    regional = stats_payload["regional"]
    inferential = stats_payload["inferential"]

    fig, axes = plt.subplots(2, 2, figsize=(FIG_WIDTH, FIG_HEIGHT))
    ax_five = axes[0][0]
    ax_three = axes[0][1]
    ax_reg_five = axes[1][0]
    ax_reg_three = axes[1][1]

    five_data = [statewide["five_bin"][bin_name]["values"] for bin_name in FIVE_BIN_ORDER]
    three_data = [statewide["three_bin"][bin_name]["values"] for bin_name in THREE_BIN_ORDER]

    ax_five.boxplot(five_data, tick_labels=FIVE_BIN_ORDER)
    ax_five.set_title("Statewide index by five-bin ONI scheme")
    ax_five.set_ylabel("Seasonal index (% of station median)")
    ax_five.tick_params(axis="x", rotation=20)
    for idx, bin_name in enumerate(FIVE_BIN_ORDER, start=1):
        n_value = statewide["five_bin"][bin_name]["n_seasons"]
        max_value = statewide["five_bin"][bin_name]["max"]
        ax_five.text(idx, max_value, f"n={n_value}", ha="center", va="bottom", fontsize=8)
    if isinstance(inferential["five_bin"]["anova"], str):
        ax_five.text(
            0.02,
            0.95,
            f"ANOVA {inferential['five_bin']['anova']}\nFailing bins: {inferential['five_bin']['failing_bins']}",
            transform=ax_five.transAxes,
            va="top",
            fontsize=8,
        )

    ax_three.boxplot(three_data, tick_labels=THREE_BIN_ORDER)
    ax_three.set_title("Statewide index by three-bin ONI scheme")
    ax_three.tick_params(axis="x", rotation=20)
    for idx, bin_name in enumerate(THREE_BIN_ORDER, start=1):
        n_value = statewide["three_bin"][bin_name]["n_seasons"]
        max_value = statewide["three_bin"][bin_name]["max"]
        ax_three.text(idx, max_value, f"n={n_value}", ha="center", va="bottom", fontsize=8)
    if isinstance(inferential["three_bin"]["anova"], str):
        ax_three.text(
            0.02,
            0.95,
            f"ANOVA {inferential['three_bin']['anova']}\nFailing bins: {inferential['three_bin']['failing_bins']}",
            transform=ax_three.transAxes,
            va="top",
            fontsize=8,
        )
    else:
        ax_three.text(
            0.02,
            0.95,
            "\n".join(
                [
                    "ANOVA "
                    f"F={format_float(inferential['three_bin']['anova']['f_statistic'])}, "
                    f"p={format_float(inferential['three_bin']['anova']['p_value'])}",
                    "Pearson "
                    f"r={format_float(inferential['three_bin']['pearson']['r'])}, "
                    f"p={format_float(inferential['three_bin']['pearson']['p_value'])}",
                ]
            ),
            transform=ax_three.transAxes,
            va="top",
            fontsize=8,
        )

    x_positions_five = list(range(1, len(FIVE_BIN_ORDER) + 1))
    for region in REGION_ORDER:
        means = [regional[region]["five_bin"][bin_name]["mean"] for bin_name in FIVE_BIN_ORDER]
        ax_reg_five.plot(x_positions_five, means, marker="o", label=region)
    ax_reg_five.set_xticks(x_positions_five)
    ax_reg_five.set_xticklabels(FIVE_BIN_ORDER, rotation=20)
    ax_reg_five.set_title("Regional seasonal index means (five-bin)")
    ax_reg_five.set_ylabel("Mean seasonal index (% of station median)")
    ax_reg_five.legend()

    x_positions_three = list(range(1, len(THREE_BIN_ORDER) + 1))
    for region in REGION_ORDER:
        means = [regional[region]["three_bin"][bin_name]["mean"] for bin_name in THREE_BIN_ORDER]
        ax_reg_three.plot(x_positions_three, means, marker="o", label=region)
    ax_reg_three.set_xticks(x_positions_three)
    ax_reg_three.set_xticklabels(THREE_BIN_ORDER, rotation=20)
    ax_reg_three.set_title("Regional seasonal index means (three-bin)")
    ax_reg_three.legend()

    fig.suptitle("Colorado April 1 SWE seasonal index vs ENSO")
    fig.tight_layout()
    fig.savefig(output_path, dpi=PLOT_DPI)
    plt.close(fig)


def print_gate_table(stats_payload: Dict[str, object]) -> None:
    gate = stats_payload["gate"]
    print("=== N-count gate (printed first by design) ===")
    for scheme_name in ("five_bin", "three_bin"):
        scheme_gate = gate[scheme_name]
        print(f"{scheme_name}:")
        for bin_name, count in scheme_gate["season_counts"].items():
            print(f"  {bin_name:24s} n={count:2d}")
        print(f"  all_bins_meet_min_seasons={scheme_gate['all_bins_meet_min_seasons']}")
        print(f"  failing_bins={scheme_gate['failing_bins']}")


def print_station_median_filter_summary(
    min_station_median_swe_in: float,
    excluded_stations: List[Dict[str, object]],
    surviving_station_count: int,
    surviving_region_counts: Dict[str, int],
) -> None:
    print("=== Minimum station median SWE filter ===")
    print(f"MIN_STATION_MEDIAN_SWE_IN={format_float(min_station_median_swe_in)}")
    if excluded_stations:
        for station in excluded_stations:
            print(
                "excluded_station="
                f"{station['triplet']} | "
                f"{station['name']} | "
                f"{station['region']} | "
                f"median={format_float(float(station['median']))}"
            )
    else:
        print("excluded_station=none")
    print(f"excluded_station_count={len(excluded_stations)}")
    print(f"surviving_station_count={surviving_station_count}")
    print(f"surviving_region_counts={surviving_region_counts}")


def main() -> None:
    stations, stations_run_utc = load_station_candidates(STATIONS_FILE)
    region_by_triplet_all = {
        str(station["triplet"]): str(station["region"]) for station in stations
    }
    station_name_by_triplet = {
        str(station["triplet"]): str(station["name"]) for station in stations
    }

    session = requests.Session()
    oni_table = fetch_oni_table(session)

    station_series_by_triplet: Dict[str, Dict[int, float]] = {}
    for station in stations:
        triplet = str(station["triplet"])
        station_series_by_triplet[triplet] = fetch_station_april1_series(
            session=session,
            station_triplet=triplet,
            start_wy=COMMON_START_WY,
        )

    observed_wys = [
        wy
        for station_series in station_series_by_triplet.values()
        for wy in station_series
        if wy >= COMMON_START_WY
    ]
    if not observed_wys:
        raise RuntimeError(
            f"No April 1 SWE observations found at or after WY{COMMON_START_WY}."
        )
    resolved_end_wy = max(observed_wys)
    oni_oct_mar_by_wy = compute_oni_oct_mar_by_wy(
        oni_table=oni_table,
        start_wy=COMMON_START_WY,
        end_wy=resolved_end_wy,
    )

    station_percent, station_medians, excluded_stations = compute_percent_of_station_median(
        station_series_by_triplet=station_series_by_triplet,
        station_name_by_triplet=station_name_by_triplet,
        region_by_triplet=region_by_triplet_all,
        min_station_median_swe_in=MIN_STATION_MEDIAN_SWE_IN,
        start_wy=COMMON_START_WY,
        end_wy=resolved_end_wy,
    )
    surviving_triplets = set(station_percent.keys())
    region_by_triplet = {
        triplet: region
        for triplet, region in region_by_triplet_all.items()
        if triplet in surviving_triplets
    }

    region_counts = {region: 0 for region in REGION_ORDER}
    for region in region_by_triplet.values():
        region_counts[region] += 1
    sparse_regions = [region for region in REGION_ORDER if region_counts[region] < 3]
    if sparse_regions:
        sparse_counts_text = ", ".join(
            [f"{region}={region_counts[region]}" for region in sparse_regions]
        )
        raise RuntimeError(
            "Minimum station median SWE filter leaves fewer than 3 stations in "
            f"region(s): {sparse_counts_text}."
        )
    (
        statewide_index_by_wy,
        statewide_counts_by_wy,
        regional_index_by_wy,
        regional_counts_by_wy,
    ) = compute_seasonal_indices(
        station_percent_by_triplet=station_percent,
        region_by_triplet=region_by_triplet,
        start_wy=COMMON_START_WY,
        end_wy=resolved_end_wy,
    )

    five_bins_by_wy = {wy: classify_five_bin(oni) for wy, oni in oni_oct_mar_by_wy.items()}
    three_bins_by_wy = {wy: classify_three_bin(oni) for wy, oni in oni_oct_mar_by_wy.items()}

    statewide_five = compute_bin_results(
        series_by_wy=statewide_index_by_wy,
        bins_by_wy=five_bins_by_wy,
        bin_order=FIVE_BIN_ORDER,
        label="statewide five-bin",
    )
    statewide_three = compute_bin_results(
        series_by_wy=statewide_index_by_wy,
        bins_by_wy=three_bins_by_wy,
        bin_order=THREE_BIN_ORDER,
        label="statewide three-bin",
    )

    regional_results: Dict[str, Dict[str, object]] = {}
    for region in REGION_ORDER:
        regional_results[region] = {
            "station_count_by_water_year": to_wy_keyed_object(regional_counts_by_wy[region]),
            "five_bin": compute_bin_results(
                series_by_wy=regional_index_by_wy[region],
                bins_by_wy=five_bins_by_wy,
                bin_order=FIVE_BIN_ORDER,
                label=f"{region} five-bin",
            ),
            "three_bin": compute_bin_results(
                series_by_wy=regional_index_by_wy[region],
                bins_by_wy=three_bins_by_wy,
                bin_order=THREE_BIN_ORDER,
                label=f"{region} three-bin",
            ),
        }

    gate_five = build_gate(statewide_five, FIVE_BIN_ORDER)
    gate_three = build_gate(statewide_three, THREE_BIN_ORDER)

    station_count = len(surviving_triplets)

    run_utc = datetime.now(timezone.utc).isoformat()
    stats_payload: Dict[str, object] = {
        "provenance": {
            "run_utc": run_utc,
            "oni_source_url": ONI_SOURCE_URL,
            "awdb_endpoint": AWDB_ENDPOINT,
            "stations_file_path": str(STATIONS_FILE.resolve()),
            "stations_file_run_utc": stations_run_utc,
            "station_count": station_count,
            "region_counts": region_counts,
            "common_start_wy": COMMON_START_WY,
            "end_wy_rule": END_WY_RULE,
            "resolved_end_wy": resolved_end_wy,
            "season_count": resolved_end_wy - COMMON_START_WY + 1,
            "oni_thresholds": ONI_THRESHOLDS,
            "min_seasons_for_test": MIN_SEASONS_FOR_TEST,
            "MIN_STATION_MEDIAN_SWE_IN": MIN_STATION_MEDIAN_SWE_IN,
            "excluded_station_count_min_median_filter": len(excluded_stations),
            "excluded_stations_min_median_filter": excluded_stations,
        },
        "gate": {
            "five_bin": gate_five,
            "three_bin": gate_three,
        },
        "statewide": {
            "five_bin": statewide_five,
            "three_bin": statewide_three,
            "station_count_by_water_year": to_wy_keyed_object(statewide_counts_by_wy),
            "seasonal_index_by_water_year": to_wy_keyed_object(statewide_index_by_wy),
            "oni_oct_mar_by_water_year": to_wy_keyed_object(oni_oct_mar_by_wy),
        },
        "regional": regional_results,
        "station_medians_april1_swe": station_medians,
    }

    stats_payload["inferential"] = {
        "five_bin": compute_inferential(
            scheme_name="five_bin",
            gate=gate_five,
            bin_results=statewide_five,
            bin_order=FIVE_BIN_ORDER,
            oni_by_wy=oni_oct_mar_by_wy,
            statewide_index_by_wy=statewide_index_by_wy,
        ),
        "three_bin": compute_inferential(
            scheme_name="three_bin",
            gate=gate_three,
            bin_results=statewide_three,
            bin_order=THREE_BIN_ORDER,
            oni_by_wy=oni_oct_mar_by_wy,
            statewide_index_by_wy=statewide_index_by_wy,
        ),
    }

    print_gate_table(stats_payload)
    print_station_median_filter_summary(
        min_station_median_swe_in=MIN_STATION_MEDIAN_SWE_IN,
        excluded_stations=excluded_stations,
        surviving_station_count=station_count,
        surviving_region_counts=region_counts,
    )
    print(f"station_count={station_count}")
    print(f"region_counts={region_counts}")
    print(f"water_year_range=WY{COMMON_START_WY}-WY{resolved_end_wy}")
    print(f"season_count={resolved_end_wy - COMMON_START_WY + 1}")

    warning = non_constant_count_warning_message(statewide_counts_by_wy, "statewide")
    if warning:
        print(warning)
    for region in REGION_ORDER:
        regional_warning = non_constant_count_warning_message(
            regional_counts_by_wy[region],
            region,
        )
        if regional_warning:
            print(regional_warning)

    stats_json_path = Path("stats.json")
    stats_text_path = Path("snowpack_enso_stats.txt")
    chart_path = Path("snowpack_enso_composite.png")

    with stats_json_path.open("w", encoding="utf-8") as handle:
        json.dump(stats_payload, handle, indent=2)
        handle.write("\n")
    write_text_report(stats_payload, stats_text_path)
    plot_results(stats_payload, chart_path)

    print(f"wrote_file={stats_json_path.resolve()}")
    print(f"wrote_file={stats_text_path.resolve()}")
    print(f"wrote_file={chart_path.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except requests.HTTPError as exc:
        raise SystemExit(f"HTTP request failed: {exc}") from exc
    except requests.RequestException as exc:
        raise SystemExit(f"Network request failed: {exc}") from exc
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(str(exc)) from exc

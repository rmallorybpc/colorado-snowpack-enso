#!/usr/bin/env python3
"""Discover candidate Colorado SNOTEL stations from AWDB metadata + data coverage.

This script does not run ENSO analysis. It only discovers and filters
stations, then writes a reproducible candidate list to stations_candidate.json.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import date, datetime, timezone
from typing import Dict, Iterable, List, Optional, Set

import requests


AWDB_BASE = "https://wcc.sc.egov.usda.gov/awdbRestApi/services/v1"
INVENTORY_ENDPOINT = "/stations"
DATA_ENDPOINT = "/data"

# Station sanity check only (these are the three stations used by the analysis).
LEGACY_CHECK_TRIPLETS = {"825:CO:SNTL", "335:CO:SNTL", "874:CO:SNTL"}

TARGET_START_WY = 1981
MAX_MISSING_APRIL1 = 2
DATA_ELEMENT = "WTEQ"


@dataclass
class StationCoverage:
    station: Dict[str, object]
    available_wys: Set[int]
    first_wy: Optional[int]
    last_wy: Optional[int]
    target_missing_count: int
    include: bool
    exclusion_reasons: List[str]


def latest_complete_water_year(today: date) -> int:
    # Water year ends on Sep 30. Before Oct 1, current WY is incomplete.
    return today.year if today.month >= 10 else today.year - 1


def fetch_inventory(session: requests.Session) -> List[Dict[str, object]]:
    """Fetch Colorado SNOTEL station metadata from AWDB."""
    url = AWDB_BASE + INVENTORY_ENDPOINT
    params = {"stationTriplets": "*:CO:SNTL"}
    print(f"Querying AWDB inventory: {url} params={params}")
    response = session.get(url, params=params, timeout=120)
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, list):
        raise RuntimeError(
            "Inventory query returned an unexpected response shape; expected a list."
        )
    if not payload:
        raise RuntimeError(
            "Inventory query returned zero stations for Colorado SNOTEL; cannot continue."
        )
    return payload


def parse_date_yyyy_mm_dd(value: str) -> date:
    return datetime.strptime(value[:10], "%Y-%m-%d").date()


def fetch_april1_available_wys(
    session: requests.Session,
    triplet: str,
    begin_date: date,
    end_date: date,
) -> Set[int]:
    """Return water years with non-null April 1 SWE values for a station."""
    url = AWDB_BASE + DATA_ENDPOINT
    params = {
        "stationTriplets": triplet,
        "elements": DATA_ELEMENT,
        "duration": "DAILY",
        "beginDate": begin_date.isoformat(),
        "endDate": end_date.isoformat(),
    }
    response = session.get(url, params=params, timeout=180)
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, list):
        raise RuntimeError(
            f"Unexpected data response shape for {triplet}; expected list payload."
        )

    available: Set[int] = set()
    for station in payload:
        for element in station.get("data", []):
            for point in element.get("values", []):
                raw_date = point.get("date")
                raw_value = point.get("value")
                if raw_date is None or raw_value is None:
                    continue
                obs_date = parse_date_yyyy_mm_dd(str(raw_date))
                if obs_date.month == 4 and obs_date.day == 1:
                    available.add(obs_date.year)
    return available


def evaluate_station(
    station: Dict[str, object],
    available_wys: Set[int],
    latest_complete_wy: int,
) -> StationCoverage:
    target_wys = set(range(TARGET_START_WY, latest_complete_wy + 1))
    first_wy = min(available_wys) if available_wys else None
    last_wy = max(available_wys) if available_wys else None
    missing_count = len(target_wys - available_wys)

    reasons: List[str] = []
    if first_wy is None:
        reasons.append("no non-null April 1 SWE values found")
    else:
        if first_wy > TARGET_START_WY:
            reasons.append(
                f"first available April 1 water year is {first_wy} (> {TARGET_START_WY})"
            )
        if last_wy < latest_complete_wy:
            reasons.append(
                f"last available April 1 water year is {last_wy} (< {latest_complete_wy})"
            )
    if missing_count > MAX_MISSING_APRIL1:
        reasons.append(
            f"missing April 1 values in WY{TARGET_START_WY}-WY{latest_complete_wy}: "
            f"{missing_count} (> {MAX_MISSING_APRIL1})"
        )

    return StationCoverage(
        station=station,
        available_wys=available_wys,
        first_wy=first_wy,
        last_wy=last_wy,
        target_missing_count=missing_count,
        include=not reasons,
        exclusion_reasons=reasons,
    )


def assign_regions_by_latitude(survivors: List[StationCoverage], latitudes: Iterable[float]) -> None:
    lat_values = list(latitudes)
    lat_min = min(lat_values)
    lat_max = max(lat_values)
    span = lat_max - lat_min
    lower_boundary = lat_min + span / 3.0
    upper_boundary = lat_min + (2.0 * span / 3.0)

    print("\nLatitude band boundaries (computed):")
    print(f"  south: lat < {lower_boundary:.6f}")
    print(f"  central: {lower_boundary:.6f} <= lat < {upper_boundary:.6f}")
    print(f"  north: lat >= {upper_boundary:.6f}")

    band_counts = {"north": 0, "central": 0, "south": 0}
    for cov in survivors:
        lat = float(cov.station["latitude"])
        if lat >= upper_boundary:
            region = "north"
        elif lat >= lower_boundary:
            region = "central"
        else:
            region = "south"
        cov.station["region"] = region
        band_counts[region] += 1

    print("\nSurvivor count by latitude band:")
    for band in ("north", "central", "south"):
        print(f"  {band}: {band_counts[band]}")


def print_exclusion(cov: StationCoverage) -> None:
    station = cov.station
    triplet = station.get("stationTriplet", "<unknown>")
    name = station.get("name", "<unknown>")
    first_wy = cov.first_wy if cov.first_wy is not None else "n/a"
    last_wy = cov.last_wy if cov.last_wy is not None else "n/a"
    print(
        f"EXCLUDED {triplet} | {name} | coverage WY{first_wy}-WY{last_wy} "
        f"| available Apr1 count={len(cov.available_wys)}"
    )
    for reason in cov.exclusion_reasons:
        print(f"  - {reason}")


def print_summary_table(survivors: List[StationCoverage]) -> None:
    print("\nSummary table (north to south):")
    print(
        "{:<15}  {:<34}  {:<7}  {:>8}  {:<14}".format(
            "triplet", "name", "region", "elev_ft", "coverage"
        )
    )
    print("-" * 90)
    for cov in sorted(survivors, key=lambda s: float(s.station["latitude"]), reverse=True):
        station = cov.station
        coverage = f"WY{cov.first_wy}-WY{cov.last_wy}"
        print(
            "{:<15}  {:<34}  {:<7}  {:>8.0f}  {:<14}".format(
                str(station["stationTriplet"]),
                str(station["name"])[:34],
                str(station["region"]),
                float(station["elevation"]),
                coverage,
            )
        )


def main() -> None:
    today = date.today()
    latest_complete_wy = latest_complete_water_year(today)
    target_end_date = date(latest_complete_wy, 6, 30)

    print(
        f"Target continuity window: WY{TARGET_START_WY}-WY{latest_complete_wy} "
        f"(through {target_end_date.isoformat()})"
    )

    session = requests.Session()

    inventory = fetch_inventory(session)
    print(f"Inventory returned {len(inventory)} stations.")

    if not inventory:
        raise RuntimeError(
            "AWDB inventory query succeeded but returned no stations; cannot continue."
        )

    seen_triplets = {str(s.get("stationTriplet", "")) for s in inventory}
    for triplet in sorted(LEGACY_CHECK_TRIPLETS):
        if triplet in seen_triplets:
            print(f"Sanity check present: {triplet}")
        else:
            print(f"Sanity check missing from discovery inventory: {triplet}")

    coverage_rows: List[StationCoverage] = []
    total = len(inventory)
    for idx, station in enumerate(inventory, start=1):
        triplet = str(station.get("stationTriplet", ""))
        name = str(station.get("name", ""))
        begin_raw = str(station.get("beginDate", ""))
        begin_date = parse_date_yyyy_mm_dd(begin_raw)
        query_begin = min(begin_date, date(TARGET_START_WY - 1, 10, 1))
        print(f"[{idx:>3}/{total}] Pulling April 1 availability for {triplet} | {name}")
        available_wys = fetch_april1_available_wys(
            session=session,
            triplet=triplet,
            begin_date=query_begin,
            end_date=target_end_date,
        )
        cov = evaluate_station(
            station=station,
            available_wys=available_wys,
            latest_complete_wy=latest_complete_wy,
        )
        coverage_rows.append(cov)
        if not cov.include:
            print_exclusion(cov)

    with_any_april1 = [c for c in coverage_rows if c.first_wy is not None]
    start_ok = [c for c in with_any_april1 if c.first_wy is not None and c.first_wy <= TARGET_START_WY]
    end_ok = [c for c in start_ok if c.last_wy is not None and c.last_wy >= latest_complete_wy]
    missing_ok = [c for c in end_ok if c.target_missing_count <= MAX_MISSING_APRIL1]
    survivors = [c for c in coverage_rows if c.include]

    for triplet in sorted(LEGACY_CHECK_TRIPLETS):
        matches = [c for c in coverage_rows if c.station.get("stationTriplet") == triplet]
        if not matches:
            continue
        cov = matches[0]
        if cov.include:
            print(f"Sanity check station passes filter: {triplet}")
        else:
            print(f"Sanity check station fails filter: {triplet}")

    if not survivors:
        print("\nNo stations survived all filters.")
    else:
        all_inventory_lats = [float(s["latitude"]) for s in inventory]
        assign_regions_by_latitude(survivors, all_inventory_lats)

    candidate_rows = []
    for cov in sorted(survivors, key=lambda s: float(s.station["latitude"]), reverse=True):
        station = cov.station
        candidate_rows.append(
            {
                "triplet": station["stationTriplet"],
                "name": station["name"],
                "latitude": float(station["latitude"]),
                "longitude": float(station["longitude"]),
                "elevation": float(station["elevation"]),
                "region": station["region"],
                "first_water_year": cov.first_wy,
                "last_water_year": cov.last_wy,
                "april1_available_count": len(cov.available_wys),
            }
        )

    output = {
        "provenance": {
            "run_utc": datetime.now(timezone.utc).isoformat(),
            "api_base": AWDB_BASE,
            "inventory_endpoint": INVENTORY_ENDPOINT,
            "data_endpoint": DATA_ENDPOINT,
            "inventory_query": {"stationTriplets": "*:CO:SNTL"},
            "filter_thresholds": {
                "target_start_water_year": TARGET_START_WY,
                "target_end_water_year": latest_complete_wy,
                "max_missing_april1_values": MAX_MISSING_APRIL1,
                "required_element": DATA_ELEMENT,
            },
            "counts": {
                "inventory_returned": len(inventory),
                "with_any_april1": len(with_any_april1),
                "start_year_ok": len(start_ok),
                "end_year_ok": len(end_ok),
                "missing_ok": len(missing_ok),
                "survivors": len(survivors),
            },
        },
        "stations": candidate_rows,
    }

    with open("stations_candidate.json", "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)
        f.write("\n")
    print(f"\nWrote stations_candidate.json with {len(candidate_rows)} stations.")

    if survivors:
        print_summary_table(survivors)


if __name__ == "__main__":
    try:
        main()
    except requests.HTTPError as exc:
        raise SystemExit(f"AWDB request failed: {exc}") from exc
    except requests.RequestException as exc:
        raise SystemExit(f"AWDB connection failed: {exc}") from exc
    except Exception as exc:  # noqa: BLE001
        raise SystemExit(f"Station discovery failed: {exc}") from exc
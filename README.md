# Colorado ENSO Explorer

Does La Niña actually matter for Colorado snowpack? 46 winters of federal
data say barely.

Live page: https://rmallorybpc.github.io/colorado-snowpack-enso/

## The finding

A strong La Niña winter and a strong El Niño winter land within half an
inch of each other on April 1 at Berthoud Summit (20.8 vs 21.2 inches
median). The only visible tilt is at Tower, near Steamboat: about 6
inches of La Niña advantage against a 39 inch year-to-year spread. The
textbook El Niño advantage in the southern mountains did not show at
Wolf Creek. Since 1981 there have been 4 strong La Niña winters and 6
strong El Niño winters. The annual forecast argument rests on that
sample.

## What is in this repo

| File | What it is |
|---|---|
| `index.html` | The findings page (GitHub Pages) |
| `colorado_snowpack_enso.py` | The full pipeline: ONI pull, ONI binning, SNOTEL pull, composite, chart, stats |
| `snowpack_enso_composite.png` | The composite chart |
| `snowpack_enso_stats.txt` | April 1 SWE medians and ranges per ENSO bin per station |
| `tmg.css`, `TMG-BRAND-GUIDE.md` | TMG design system |

## Method

Each winter since water year 1981 is classified by its December to
February Oceanic Niño Index value into five bins, from strong La Niña
to strong El Niño. Daily snow water equivalent from three long-record
SNOTEL stations (Tower, Berthoud Summit, Wolf Creek Summit) is cut into
October to June seasonal curves. Seasons are grouped by bin. The chart
shows every individual season with the group medians on top. This is a
historical analog composite, not a forecast.

## Data sources

All free, public, and federal.

- ONI: NOAA Climate Prediction Center, https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt
- Snowpack: NRCS SNOTEL via the AWDB REST API

## Reproduce it

```
pip install requests pandas matplotlib
python colorado_snowpack_enso.py
```

The script prints the season count per ONI bin first. That count is the
honest constraint on everything downstream: the strong bins hold 4 and
6 winters.

## Honest limits

Four strong La Niña winters is a thin sample. Three stations are not
the whole state. April 1 is one snapshot of a season. A null here does
not rule out ENSO effects elsewhere in the West, where published
signals are stronger.

# TMG Brand Guide

**Version 1.0 · The Media Guild**

---

## Purpose

This document defines the visual and editorial standards for all TMG
research tools and findings pages.  Drop this file in the root of every
TMG tool repository so contributors always have the spec at hand.

---

## Color Palette

| Token | Hex | Usage |
|-------|-----|-------|
| `--tmg-cream` | `#FDFAF4` | Page background |
| `--tmg-cream-dark` | `#F3EFE6` | Alternating rows, stat blocks |
| `--tmg-white` | `#FFFFFF` | Card backgrounds |
| `--tmg-navy` | `#1B2A3B` | Nav bar, headings, footer |
| `--tmg-navy-mid` | `#2E4160` | Hero gradient, hover states |
| `--tmg-accent` | `#C8512A` | CTAs, links, section labels, finding border |
| `--tmg-accent-hover` | `#A33E1C` | Hover state for accent elements |
| `--tmg-text` | `#1A1A1A` | Body copy |
| `--tmg-text-muted` | `#5A5A5A` | Captions, helper text |
| `--tmg-border` | `#DDD8CE` | Card borders, table rules |

---

## Typography

| Role | Family | Weight | Size |
|------|--------|--------|------|
| Body | Inter / Segoe UI / Arial | 400 | 16 px |
| Nav / labels | Inter | 700 | 0.875 rem |
| Headings | Georgia / Times New Roman | 400 | clamp(1.65 rem, 4 vw, 2.75 rem) |
| Card titles | Inter | 700 | 1.2 rem |
| Captions | Inter | 400 | 0.82 rem |

---

## Spacing

| Scale | Value |
|-------|-------|
| Base unit | 8 px |
| Card padding | 2 rem × 2.25 rem |
| Section gap | 2 rem (bottom margin on cards) |
| Container max-width | 1100 px |
| Nav height | 56 px |

---

## Components

### Navigation (`tmg-nav`)
- Sticky, navy background, `z-index: 1000`
- Logo left: `TMG Research` in all-caps, bold
- Links right: plain text links + Tool Suite dropdown
- Dropdown uses `data-open` attribute toggled by JS; arrow rotates 180°

### Hero (`tmg-hero`)
- Navy-to-navy-mid diagonal gradient
- Kicker in uppercase, 60 % opacity white
- Headline in serif, `clamp` sizing
- Dek at 78 % opacity, 640 px max-width

### Cards (`tmg-card`)
- White background, 1 px border `--tmg-border`
- 14 px border-radius
- `box-shadow: 0 2px 12px rgba(0,0,0,.08)`

### Key Finding Block (`tmg-finding`)
- 4 px left border in `--tmg-accent`
- Background `#FFF7F4`
- Use `<strong>` for the accent-colored emphasis phrase

### Phase Pills (`tmg-phase-pill`)
- Full-bleed colored pill per ENSO phase
- El Niño: `#C0473A` · Neutral: `#5578A8` · La Niña: `#3A7F92`

### Stat Boxes (`tmg-stat`)
- Cream-dark background, centered
- Large bold value, small muted label

### Data Table (`tmg-table`)
- Full-width, cream-dark header row
- 2 px bottom border on `<th>`, 1 px on body rows
- Hover: cream-dark row highlight

---

## Chart / Figure Standards

- Save charts at **150 dpi** PNG with cream (`#FDFAF4`) background
- Panel backgrounds: `#F8F4EC` (slightly darker than page cream)
- Caption below figure, centered, 0.82 rem, `--tmg-text-muted`
- Include a `role="img"` alt description on every `<img>`

---

## File Naming Convention

| File | Name |
|------|------|
| Findings page | `index.html` |
| Stylesheet | `tmg.css` |
| Brand guide | `TMG-BRAND-GUIDE.md` |
| Deploy workflow | `.github/workflows/deploy.yml` |
| Python pipeline | `<tool-name>.py` (snake_case) |
| Chart output | `<tool-name>_composite.png` |
| Stats output | `<tool-name>_stats.txt` |

---

## Editorial Voice

- Lead with the finding, not the method
- Use plain language; avoid jargon without definition
- Quantify uncertainty (p-values, confidence intervals, sample size)
- Credit all data sources inline and in the footer
- Sentence case for headlines; Title Case for section labels

---

## Accessibility

- Minimum contrast ratio 4.5 : 1 for body text on backgrounds
- All charts need a descriptive `alt` attribute
- Dropdown buttons must have `aria-expanded` and `aria-haspopup`
- Navigation landmark must have `aria-label`

---

## Repository Conventions

Every TMG tool repository should have this flat root structure:

```
<tool-name>/
├─ .github/
│  └─ workflows/
│     └─ deploy.yml
├─ index.html
├─ tmg.css
├─ TMG-BRAND-GUIDE.md
├─ README.md
├─ <tool-name>.py
├─ <tool-name>_composite.png
└─ <tool-name>_stats.txt
```

---

*This guide is maintained in the root of each TMG tool repository.
For updates, open an issue or PR in the relevant repository.*

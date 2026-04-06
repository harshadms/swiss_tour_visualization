# UI Redesign — Swiss Tour Visualization
**Date:** 2026-04-06  
**Status:** Approved

---

## Summary

Redesign `index.html` from a bare Folium-generated full-screen map into a polished single-page visualization with dark/light themes, a collapsible bottom drawer, and ride chip navigation. All changes are implemented as a post-processing step that wraps/augments the Folium-generated output — `read_plot.py` remains the data layer, the new UI layer is injected into or layered on top of the generated HTML.

---

## Layout

**Layout A — Full-screen map + bottom drawer.**

- The Leaflet map fills the entire viewport (`position: fixed; inset: 0`).
- A bottom drawer sits above the map as an overlay (`position: fixed; bottom: 0; left: 0; right: 0`).
- Two states:
  1. **Collapsed** — chips bar (~48px tall): scrollable row of ride chips, one per ride.
  2. **Expanded** — full drawer (~260px): elevation chart (full width) + icon stat row below it.
- Clicking a ride chip OR clicking a polyline on the map selects that ride and expands the drawer.
- Clicking the drawer handle (or pressing Escape) collapses it.
- The map bottom padding adjusts dynamically so the drawer never permanently obscures Switzerland.

---

## Theme

Two themes ship in the same HTML, toggled by a button in the top-left badge area.

| Token | Dark | Light |
|---|---|---|
| Map tile | CartoDB Dark Matter | CartoDB Positron |
| Background | `#0a1422` | `#f5f0e8` |
| Accent / rides | `#ff6b35` (orange) | `#c0392b` (red) |
| Text primary | `#ffffff` | `#3a2a1a` |
| Text secondary | `#8ab4d4` | `#6b5a3e` |
| Drawer bg | `#0a1422` | `#f5f0e8` |
| Drawer border-top | `2px solid accent` | `2px solid accent` |
| Canton completed fill | `rgba(255,107,53,0.18)` | `rgba(192,57,43,0.12)` |
| Canton incomplete fill | `rgba(138,180,212,0.1)` | `rgba(100,150,80,0.1)` |

Theme is stored in `localStorage` and applied as a `data-theme` attribute on `<html>`. All colors use CSS variables — a single attribute swap changes everything.

The Leaflet tile layer is swapped in JS when the theme toggles.

---

## Floating Overlays (replace Leaflet's default controls)

| Element | Position | Contents |
|---|---|---|
| Title badge | top-left | "🚴 SWISS TOUR" + theme toggle icon |
| Layer toggles | top-right | Checkboxes: Cantons, Districts, Rides, Stops (replaces Leaflet's collapsed layer control) |

Both styled as frosted-glass cards matching the active theme.

---

## Collapsed Drawer — Ride Chips

- One chip per ride, scrollable horizontally.
- Each chip: `● Ride N · Canton` (colored dot + short label).
- Active chip highlighted with accent color + border.
- A faint "← select a ride" hint fades in when no chip is active.
- Chips are generated from `rides_info.csv` data already embedded in the page.

---

## Expanded Drawer — Ride Detail (Layout B)

```
┌─────────────────────────────────────────────────────────┐
│ RIDE 8 — GRAUBÜNDEN          Zuoz → Chur   142 km   ↓  │  ← handle (click to collapse)
├─────────────────────────────────────────────────────────┤
│  [full-width elevation chart with gradient fill]        │
├──────────┬──────────┬──────────┬──────────┬─────────────┤
│ 📏 142km │ ⛰ 2840m │ 🚴 28kph │ ⚡ 55kph │ 📈 4.2%     │
└──────────┴──────────┴──────────┴──────────┴─────────────┘
```

- Elevation chart: inline SVG generated from the same elevation data used for the PNG plots. No base64 PNG embed needed in the drawer — the SVG is lightweight and theme-aware (line/fill color follow CSS variables).
- Stat row: 5–6 icon+value+label cells, evenly spaced, separated by 1px dividers.
- The existing `generate_stats_from_gps()` output (max elev, min elev, avg/max/min speed, avg/max gradient) maps directly to these cells.

---

## Implementation Approach

The redesign is done as a **post-processing step in `read_plot.py`**:

1. Folium generates the map as before and saves raw HTML to a temp string (not directly to `index.html`).
2. A new `inject_ui(raw_html)` function:
   - Extracts the `<head>` and `<body>` content from Folium's output.
   - Injects CSS variables, theme styles, drawer HTML, and the JS controller into the `<head>` / end of `<body>`.
   - Replaces the Leaflet default layer control with the custom overlay.
   - Writes the final combined HTML to `index.html`.
3. Ride data (names, stats, elevation points) is serialised to a `const RIDES = [...]` JSON block injected into the `<script>` section, so the drawer JS can render chips and charts without parsing the Leaflet layer data.

This keeps `read_plot.py` as the single source of truth for data, and the UI layer is cleanly separable.

---

## Files Changed

| File | Change |
|---|---|
| `read_plot.py` | Add `inject_ui()`, serialise ride data to JS, call inject at end instead of `map1.save()` |
| `index.html` | Regenerated output — not hand-edited |

No new dependencies. No new data files.

---

## Out of Scope

- Mobile responsiveness (drawer behaviour on small screens)
- Ride filtering / search
- Progress tracker / canton completion percentage bar
- Fixing the elevation plot filename bug (separate issue)

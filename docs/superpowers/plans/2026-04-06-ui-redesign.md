# UI Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Transform the Folium-generated `index.html` into a polished single-page visualization with dark/light themes, a collapsible bottom drawer with ride chips and elevation chart, and custom floating overlays replacing Leaflet's default controls.

**Architecture:** A new `inject_ui(raw_html, rides_data)` function in `read_plot.py` post-processes Folium's HTML output — injecting CSS variables, the drawer HTML, and a `const RIDES` JSON block. Folium generates map data as before; `inject_ui` adds the UI shell on top. The final HTML is written to `index.html` via `inject_ui` instead of `map1.save()`.

**Tech Stack:** Python (string manipulation, json, re), vanilla JS + CSS (no new dependencies), Leaflet (already present via Folium), CartoDB tile layers (already available).

---

## File Structure

| File | Role |
|---|---|
| `read_plot.py` | Add `build_rides_data()` to collect per-ride stats/elevation during the ride loop, add `inject_ui(raw_html, rides_data)` to post-process HTML, change final save to call `inject_ui`. |
| `ui_template.py` | New file. Contains `CSS_TEMPLATE`, `HTML_TEMPLATE`, `JS_TEMPLATE` as Python string constants. Keeps `read_plot.py` readable — large HTML/CSS/JS blocks live here. |
| `index.html` | Regenerated output — never hand-edited. |

---

## Task 1: Add `build_rides_data()` to collect per-ride data

**Files:**
- Modify: `read_plot.py` (the ride loop starting at line 375)

This function replaces the ad-hoc per-ride variable handling with a structured dict that `inject_ui` can later serialise to JS.

- [ ] **Step 1: Add a `ride_data_list` collector before the ride loop**

In `read_plot.py`, find the line `for ride_id in rides_info.id.unique():` (line 375). Add this line immediately before it:

```python
ride_data_list = []
```

- [ ] **Step 2: Build a per-ride dict inside the loop and append it**

Inside the ride loop, after the existing line:
```python
stats_html = generate_stats_from_gps(ride_coords, ride_ele, ride_timestamps)
```

Add:

```python
    # Build structured data for UI injection
    ride_towns_list = list(ride.towns.values)
    # Derive canton from bezirk column (first row's bezirk maps to canton name)
    ride_bezirke = list(ride.bezirk.values)
    ride_canton = rides_info.query(f"id == {ride_id}").bezirk.values[0]
    # Map bezirk → canton using the completed_cantons list
    canton_label = next(
        (c for c in completed_cantons if c in ride_canton or ride_canton in c),
        ride_canton
    )

    # Elevation SVG data: downsample to 200 points max for small SVG
    ele_step = max(1, len(ride_ele) // 200)
    ele_sampled = ride_ele[::ele_step]
    coords_sampled = ride_coords[::ele_step]
    distances_sampled = [0.0]
    for i in range(1, len(coords_sampled)):
        d = geodesic(coords_sampled[i - 1], coords_sampled[i]).km
        distances_sampled.append(distances_sampled[-1] + d)

    # Stats dict (reuse logic from generate_stats_from_gps but return raw values)
    ride_stats = _compute_stats(ride_coords, ride_ele, ride_timestamps)

    ride_data_list.append({
        "id": int(ride_id),
        "label": f"Ride {ride_id}",
        "canton": canton_label,
        "towns": ride_towns_list,
        "distances_km": [round(d, 2) for d in distances_sampled],
        "elevations_m": [round(e, 1) if e is not None else 0 for e in ele_sampled],
        "stats": ride_stats,
    })
```

- [ ] **Step 3: Add `_compute_stats()` helper above the ride loop**

This extracts the numeric stats from GPS data (reusing the same logic as `generate_stats_from_gps` but returning a plain dict instead of HTML). Add this function definition after `generate_stats_from_gps` (around line 271):

```python
def _compute_stats(coords, elevations, timestamps):
    """Return a plain dict of ride stats (no HTML). Used for JS data injection."""
    import numpy as np
    raw_speeds = []
    raw_gradients = []
    window_size = 5
    min_dist = 2
    min_time = 1
    ele_copy = list(elevations)

    for i in range(1, len(coords)):
        p1, p2 = coords[i - 1], coords[i]
        lat1, lon1 = p1
        lat2, lon2 = p2
        ele1 = ele_copy[i - 1] if i - 1 < len(ele_copy) else None
        ele2 = ele_copy[i] if i < len(ele_copy) else None
        t1 = timestamps[i - 1]
        t2 = timestamps[i]
        if None in (lat1, lon1, ele1, t1, lat2, lon2, ele2, t2):
            continue
        dist_m = geodesic((lat1, lon1), (lat2, lon2)).meters
        time_s = (t2 - t1).total_seconds()
        if dist_m < min_dist or time_s < min_time:
            continue
        speed_kmh = (dist_m / time_s) * 3.6
        gradient = (ele2 - ele1) / dist_m * 100
        if abs(gradient) < 100:
            raw_gradients.append(gradient)
        raw_speeds.append(speed_kmh)

    smooth_speeds = moving_average(raw_speeds, window_size)
    smooth_gradients = moving_average(raw_gradients, window_size)

    # Total distance
    total_km = sum(
        geodesic(coords[i - 1], coords[i]).km
        for i in range(1, len(coords))
    )

    valid_ele = [e for e in elevations if e is not None]
    return {
        "distance_km": round(total_km, 1),
        "max_elevation_m": round(max(valid_ele), 0) if valid_ele else None,
        "min_elevation_m": round(min(valid_ele), 0) if valid_ele else None,
        "avg_speed_kmh": round(float(np.mean(smooth_speeds)), 1) if len(smooth_speeds) else None,
        "max_speed_kmh": round(float(max(smooth_speeds)), 1) if len(smooth_speeds) else None,
        "avg_gradient_pct": round(float(np.mean(smooth_gradients)), 1) if len(smooth_gradients) else None,
        "max_gradient_pct": round(float(max(smooth_gradients)), 1) if len(smooth_gradients) else None,
    }
```

- [ ] **Step 4: Manual smoke test**

Run:
```bash
cd /Users/harshad/P_Projects/swiss_tour_visualization
source env/bin/activate
python - <<'EOF'
import sys
sys.path.insert(0, '.')
# Quickly verify _compute_stats is importable and ride_data_list structure
# by running just the helpers (not the full script)
from read_plot import _compute_stats, moving_average
print("_compute_stats importable OK")
EOF
```

Expected output: `_compute_stats importable OK`

> Note: `read_plot.py` is a script not a module, so direct import runs the whole thing. The step above uses a path trick — if it doesn't work, just run `python read_plot.py` and confirm no errors instead.

- [ ] **Step 5: Commit**

```bash
git add read_plot.py
git commit -m "feat: collect per-ride structured data for UI injection"
```

---

## Task 2: Create `ui_template.py` — CSS

**Files:**
- Create: `ui_template.py`

All CSS lives here as a Python string. No logic — just the style sheet.

- [ ] **Step 1: Create `ui_template.py` with the CSS constant**

```python
# ui_template.py
# CSS, HTML shell, and JS for the swiss-tour UI overlay.
# Injected into Folium's output by inject_ui() in read_plot.py.

CSS_TEMPLATE = """
<style>
/* ── CSS variables ─────────────────────────────────────── */
:root, [data-theme="dark"] {
  --bg:        #0a1422;
  --bg2:       #0f1e30;
  --accent:    #ff6b35;
  --text1:     #ffffff;
  --text2:     #8ab4d4;
  --border:    #1e2d4a;
  --chip-bg:   #0f1e30;
  --overlay-bg: rgba(10,20,35,0.88);
}
[data-theme="light"] {
  --bg:        #f5f0e8;
  --bg2:       #ece6da;
  --accent:    #c0392b;
  --text1:     #3a2a1a;
  --text2:     #6b5a3e;
  --border:    #d4c9b8;
  --chip-bg:   #ece6da;
  --overlay-bg: rgba(245,240,232,0.92);
}

/* ── Reset map to fill viewport ────────────────────────── */
html, body { margin:0; padding:0; width:100%; height:100%; overflow:hidden; }
.folium-map { position:fixed !important; inset:0 !important; z-index:1; }

/* ── Title badge (top-left) ─────────────────────────────── */
#sw-title {
  position:fixed; top:12px; left:12px; z-index:1000;
  background:var(--overlay-bg);
  border:1px solid var(--border);
  border-radius:8px;
  padding:6px 12px;
  display:flex; align-items:center; gap:10px;
  backdrop-filter:blur(6px);
  -webkit-backdrop-filter:blur(6px);
}
#sw-title .title-text {
  font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
  font-size:12px; font-weight:700; letter-spacing:.8px;
  color:var(--accent);
}
#sw-theme-toggle {
  background:none; border:1px solid var(--border);
  border-radius:20px; padding:2px 8px;
  font-size:11px; cursor:pointer; color:var(--text2);
  transition:border-color .2s;
}
#sw-theme-toggle:hover { border-color:var(--accent); }

/* ── Layer control (top-right) ──────────────────────────── */
#sw-layers {
  position:fixed; top:12px; right:12px; z-index:1000;
  background:var(--overlay-bg);
  border:1px solid var(--border);
  border-radius:8px;
  padding:8px 12px;
  backdrop-filter:blur(6px);
  -webkit-backdrop-filter:blur(6px);
  display:flex; flex-direction:column; gap:5px;
}
#sw-layers label {
  font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
  font-size:11px; color:var(--text2);
  display:flex; align-items:center; gap:6px; cursor:pointer;
  white-space:nowrap;
}
#sw-layers input[type=checkbox] { accent-color:var(--accent); cursor:pointer; }

/* ── Bottom drawer ──────────────────────────────────────── */
#sw-drawer {
  position:fixed; bottom:0; left:0; right:0; z-index:1000;
  background:var(--bg);
  border-top:2px solid var(--accent);
  font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;
  transition:height .25s ease;
  overflow:hidden;
}
#sw-drawer.collapsed { height:48px; }
#sw-drawer.expanded  { height:260px; }

/* chips bar (always visible at bottom of drawer) */
#sw-chips-bar {
  height:48px;
  display:flex; align-items:center; gap:8px;
  padding:0 14px;
  overflow-x:auto;
  scrollbar-width:none;
  flex-shrink:0;
}
#sw-chips-bar::-webkit-scrollbar { display:none; }
#sw-hint {
  font-size:10px; color:var(--border);
  white-space:nowrap; flex-shrink:0;
  font-style:italic;
  transition:opacity .3s;
}

.sw-chip {
  display:flex; align-items:center; gap:5px;
  background:var(--chip-bg);
  border:1px solid var(--border);
  border-radius:20px;
  padding:4px 11px;
  font-size:10px; color:var(--text2);
  white-space:nowrap; flex-shrink:0;
  cursor:pointer;
  transition:border-color .15s, color .15s, background .15s;
}
.sw-chip:hover { border-color:var(--accent); color:var(--text1); }
.sw-chip.active {
  border-color:var(--accent);
  background:rgba(255,107,53,0.12);
  color:var(--accent);
}
[data-theme="light"] .sw-chip.active { background:rgba(192,57,43,0.1); }
.sw-chip .dot {
  width:6px; height:6px; border-radius:50%;
  background:var(--text2); flex-shrink:0;
  transition:background .15s;
}
.sw-chip.active .dot { background:var(--accent); }

/* expanded content */
#sw-drawer-content {
  height:212px; /* 260 - 48 */
  display:flex; flex-direction:column;
  overflow:hidden;
}
#sw-drawer.collapsed #sw-drawer-content { display:none; }

#sw-handle {
  display:flex; align-items:center; justify-content:space-between;
  padding:8px 16px;
  border-bottom:1px solid var(--border);
  cursor:pointer; flex-shrink:0;
}
#sw-handle .ride-name {
  font-size:11px; font-weight:700; letter-spacing:.5px; color:var(--accent);
}
#sw-handle .ride-subtitle { font-size:9px; color:var(--text2); margin-top:2px; }
#sw-handle .handle-right {
  display:flex; align-items:center; gap:10px;
}
#sw-handle .distance-pill {
  font-size:9px; color:var(--text2);
  background:var(--chip-bg);
  border:1px solid var(--border);
  border-radius:10px; padding:2px 8px;
}
#sw-handle .chevron { color:var(--accent); font-size:13px; }

/* elevation chart */
#sw-chart-wrap {
  flex:1; padding:8px 16px 4px;
  min-height:0;
}
#sw-elev-svg {
  width:100%; height:100%;
  display:block;
}
.sw-elev-line { fill:none; stroke:var(--accent); stroke-width:1.8; stroke-linecap:round; stroke-linejoin:round; }
.sw-elev-fill { fill:var(--accent); opacity:.15; stroke:none; }

/* stat row */
#sw-stat-row {
  display:flex;
  border-top:1px solid var(--border);
  flex-shrink:0;
}
.sw-stat {
  flex:1; text-align:center;
  padding:5px 4px;
  border-right:1px solid var(--border);
}
.sw-stat:last-child { border-right:none; }
.sw-stat .s-icon { font-size:11px; }
.sw-stat .s-val  { font-size:12px; font-weight:700; color:var(--text1); margin-top:1px; }
.sw-stat .s-unit { font-size:8px; color:var(--text2); }
.sw-stat .s-lbl  { font-size:7px; color:var(--text2); text-transform:uppercase; letter-spacing:.3px; margin-top:1px; }
</style>
"""
```

- [ ] **Step 2: Commit**

```bash
git add ui_template.py
git commit -m "feat: add CSS template for UI overlay"
```

---

## Task 3: Add HTML shell to `ui_template.py`

**Files:**
- Modify: `ui_template.py`

- [ ] **Step 1: Append `HTML_TEMPLATE` to `ui_template.py`**

```python
HTML_TEMPLATE = """
<!-- ── Title badge ───────────────────────── -->
<div id="sw-title">
  <span class="title-text">🚴 SWISS TOUR</span>
  <button id="sw-theme-toggle" onclick="swToggleTheme()">☀️ Light</button>
</div>

<!-- ── Layer controls ────────────────────── -->
<div id="sw-layers">
  <label><input type="checkbox" id="sw-layer-cantons"   checked onchange="swToggleLayer('Cantons',   this.checked)">Cantons</label>
  <label><input type="checkbox" id="sw-layer-completed" checked onchange="swToggleLayer('Completed', this.checked)">Completed</label>
  <label><input type="checkbox" id="sw-layer-districts" checked onchange="swToggleLayer('Districts', this.checked)">Districts</label>
  <label><input type="checkbox" id="sw-layer-rides"     checked onchange="swToggleLayer('Rides',     this.checked)">Rides</label>
  <label><input type="checkbox" id="sw-layer-stops"     checked onchange="swToggleLayer('Stops',     this.checked)">Stops</label>
</div>

<!-- ── Bottom drawer ─────────────────────── -->
<div id="sw-drawer" class="collapsed">

  <!-- Chips bar (collapsed state) -->
  <div id="sw-chips-bar">
    <span id="sw-hint">← select a ride</span>
    <!-- chips injected by JS -->
  </div>

  <!-- Expanded content -->
  <div id="sw-drawer-content">
    <div id="sw-handle" onclick="swCollapseDrawer()">
      <div>
        <div class="ride-name"  id="sw-ride-name">—</div>
        <div class="ride-subtitle" id="sw-ride-subtitle">—</div>
      </div>
      <div class="handle-right">
        <span class="distance-pill" id="sw-distance">—</span>
        <span class="chevron">↓</span>
      </div>
    </div>
    <div id="sw-chart-wrap">
      <svg id="sw-elev-svg" viewBox="0 0 800 80" preserveAspectRatio="none">
        <polygon class="sw-elev-fill" id="sw-elev-fill" points=""/>
        <polyline class="sw-elev-line" id="sw-elev-line" points=""/>
      </svg>
    </div>
    <div id="sw-stat-row">
      <div class="sw-stat"><div class="s-icon">📏</div><div class="s-val" id="sw-s-dist">—</div><div class="s-lbl">Distance</div></div>
      <div class="sw-stat"><div class="s-icon">⛰</div><div class="s-val" id="sw-s-maxele">—</div><div class="s-lbl">Max Elev</div></div>
      <div class="sw-stat"><div class="s-icon">🚴</div><div class="s-val" id="sw-s-avgspd">—</div><div class="s-lbl">Avg Speed</div></div>
      <div class="sw-stat"><div class="s-icon">⚡</div><div class="s-val" id="sw-s-maxspd">—</div><div class="s-lbl">Max Speed</div></div>
      <div class="sw-stat"><div class="s-icon">📈</div><div class="s-val" id="sw-s-avggrad">—</div><div class="s-lbl">Avg Grade</div></div>
      <div class="sw-stat"><div class="s-icon">⬆</div><div class="s-val" id="sw-s-maxgrad">—</div><div class="s-lbl">Max Grade</div></div>
    </div>
  </div>
</div>
"""
```

- [ ] **Step 2: Commit**

```bash
git add ui_template.py
git commit -m "feat: add HTML shell template for drawer and overlays"
```

---

## Task 4: Add JS controller to `ui_template.py`

**Files:**
- Modify: `ui_template.py`

- [ ] **Step 1: Append `JS_TEMPLATE` to `ui_template.py`**

```python
JS_TEMPLATE = """
<script>
(function() {
  // ── Theme ──────────────────────────────────────────────────
  var _darkTile  = null;
  var _lightTile = null;
  var _map       = null;

  function swFindMap() {
    // Folium stores the map on window under a generated name; find it
    for (var k in window) {
      if (window[k] && window[k]._leaflet_id !== undefined && typeof window[k].addLayer === 'function') {
        return window[k];
      }
    }
    return null;
  }

  function swInitTiles() {
    if (_darkTile) return;
    _map = swFindMap();
    if (!_map) return;
    _darkTile  = L.tileLayer('https://{s}.basemaps.cartocdn.com/dark_all/{z}/{x}/{y}{r}.png', {
      attribution: '&copy; <a href="https://carto.com/">CARTO</a>',
      subdomains: 'abcd', maxZoom: 19
    });
    _lightTile = L.tileLayer('https://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}{r}.png', {
      attribution: '&copy; <a href="https://carto.com/">CARTO</a>',
      subdomains: 'abcd', maxZoom: 19
    });
  }

  window.swToggleTheme = function() {
    var html  = document.documentElement;
    var isDark = html.getAttribute('data-theme') !== 'light';
    var next  = isDark ? 'light' : 'dark';
    html.setAttribute('data-theme', next);
    localStorage.setItem('sw-theme', next);
    document.getElementById('sw-theme-toggle').textContent = next === 'dark' ? '☀️ Light' : '🌙 Dark';

    swInitTiles();
    if (_map) {
      if (next === 'dark') {
        _map.eachLayer(function(l) { if (l._url && l._url.includes('light_all')) _map.removeLayer(l); });
        _darkTile.addTo(_map);
      } else {
        _map.eachLayer(function(l) { if (l._url && l._url.includes('dark_all')) _map.removeLayer(l); });
        _lightTile.addTo(_map);
      }
    }
  };

  // Apply saved theme on load
  var savedTheme = localStorage.getItem('sw-theme') || 'dark';
  document.documentElement.setAttribute('data-theme', savedTheme);
  if (savedTheme === 'light') {
    document.addEventListener('DOMContentLoaded', function() {
      var btn = document.getElementById('sw-theme-toggle');
      if (btn) btn.textContent = '🌙 Dark';
    });
  }

  // ── Layer control ──────────────────────────────────────────
  // Folium stores feature groups by name; we find them in the layer control
  function swGetLayerByName(name) {
    var lc = swFindLayerControl();
    if (!lc) return null;
    var overlays = lc._layers;
    for (var i = 0; i < overlays.length; i++) {
      if (overlays[i].name === name) return overlays[i].layer;
    }
    return null;
  }

  function swFindLayerControl() {
    if (!_map) _map = swFindMap();
    if (!_map) return null;
    var controls = _map._controlContainer;
    // Access Leaflet's internal control list
    for (var k in _map._controls || {}) {
      if (_map._controls[k] && _map._controls[k]._layers) return _map._controls[k];
    }
    // fallback: iterate _leaflet_events
    var found = null;
    _map.eachLayer && _map.eachLayer(function() {});
    // Leaflet stores controls in _controlCorners, access via global
    for (var gk in window) {
      if (window[gk] && window[gk]._layers && window[gk].addTo) { found = window[gk]; }
    }
    return found;
  }

  window.swToggleLayer = function(name, visible) {
    if (!_map) _map = swFindMap();
    if (!_map) return;
    // Walk all layers on the map to find the feature group by Folium name
    _map.eachLayer(function(layer) {
      if (layer.options && layer.options.name === name) {
        if (visible) { _map.addLayer(layer); } else { _map.removeLayer(layer); }
      }
    });
    // Folium uses a LayerControl; toggle via it
    for (var k in window) {
      var obj = window[k];
      if (obj && obj._layers && Array.isArray(obj._layers)) {
        obj._layers.forEach(function(entry) {
          if (entry.name === name) {
            if (visible) { _map.addLayer(entry.layer); }
            else         { _map.removeLayer(entry.layer); }
          }
        });
      }
    }
  };

  // Also hide Leaflet's default layer control (we replace it)
  document.addEventListener('DOMContentLoaded', function() {
    var lc = document.querySelector('.leaflet-control-layers');
    if (lc) lc.style.display = 'none';
  });

  // ── Ride chips ─────────────────────────────────────────────
  var _activeRideId = null;

  function swBuildChips() {
    if (typeof RIDES === 'undefined') return;
    var bar = document.getElementById('sw-chips-bar');
    if (!bar) return;
    RIDES.forEach(function(ride) {
      var chip = document.createElement('div');
      chip.className = 'sw-chip';
      chip.dataset.rideId = ride.id;
      chip.innerHTML = '<span class="dot"></span>' + ride.label + ' &middot; ' + ride.canton;
      chip.addEventListener('click', function() { swSelectRide(ride.id); });
      bar.appendChild(chip);
    });
  }

  window.swSelectRide = function(rideId) {
    if (typeof RIDES === 'undefined') return;
    var ride = RIDES.find(function(r) { return r.id === rideId; });
    if (!ride) return;

    _activeRideId = rideId;

    // Update chips
    document.querySelectorAll('.sw-chip').forEach(function(c) {
      c.classList.toggle('active', parseInt(c.dataset.rideId) === rideId);
    });

    // Hide hint
    var hint = document.getElementById('sw-hint');
    if (hint) hint.style.opacity = '0';

    // Populate handle
    document.getElementById('sw-ride-name').textContent =
      ride.label.toUpperCase() + ' \u2014 ' + ride.canton.toUpperCase();
    document.getElementById('sw-ride-subtitle').textContent =
      ride.towns.length > 1 ? ride.towns[0] + ' \u2192 ' + ride.towns[ride.towns.length - 1] : ride.towns[0] || '';
    var dist = ride.stats.distance_km;
    document.getElementById('sw-distance').textContent = dist ? dist + ' km' : '—';

    // Populate stat row
    function fmt(v, unit) { return v != null ? v + '<span class="s-unit"> ' + unit + '</span>' : '—'; }
    document.getElementById('sw-s-dist').innerHTML    = fmt(ride.stats.distance_km, 'km');
    document.getElementById('sw-s-maxele').innerHTML  = fmt(ride.stats.max_elevation_m, 'm');
    document.getElementById('sw-s-avgspd').innerHTML  = fmt(ride.stats.avg_speed_kmh, 'km/h');
    document.getElementById('sw-s-maxspd').innerHTML  = fmt(ride.stats.max_speed_kmh, 'km/h');
    document.getElementById('sw-s-avggrad').innerHTML = fmt(ride.stats.avg_gradient_pct, '%');
    document.getElementById('sw-s-maxgrad').innerHTML = fmt(ride.stats.max_gradient_pct, '%');

    // Draw elevation SVG
    swDrawElevation(ride);

    // Expand drawer
    document.getElementById('sw-drawer').className = 'expanded';
  };

  function swDrawElevation(ride) {
    var dists = ride.distances_km;
    var eles  = ride.elevations_m;
    if (!dists || !dists.length) return;

    var W = 800, H = 80;
    var minE = Math.min.apply(null, eles);
    var maxE = Math.max.apply(null, eles);
    var rangeE = maxE - minE || 1;
    var maxD = dists[dists.length - 1] || 1;

    var pts = [];
    for (var i = 0; i < dists.length; i++) {
      var x = (dists[i] / maxD) * W;
      var y = H - ((eles[i] - minE) / rangeE) * (H - 8);
      pts.push(x.toFixed(1) + ',' + y.toFixed(1));
    }
    var lineStr = pts.join(' ');
    var fillStr = lineStr + ' ' + W + ',' + H + ' 0,' + H;

    document.getElementById('sw-elev-line').setAttribute('points', lineStr);
    document.getElementById('sw-elev-fill').setAttribute('points', fillStr);
  }

  window.swCollapseDrawer = function() {
    document.getElementById('sw-drawer').className = 'collapsed';
    var hint = document.getElementById('sw-hint');
    if (hint) hint.style.opacity = '1';
    _activeRideId = null;
    document.querySelectorAll('.sw-chip').forEach(function(c) { c.classList.remove('active'); });
  };

  // Keyboard shortcut: Escape to collapse
  document.addEventListener('keydown', function(e) {
    if (e.key === 'Escape') swCollapseDrawer();
  });

  // Init on DOM ready
  document.addEventListener('DOMContentLoaded', function() {
    swBuildChips();
    // Also init tiles so theme toggle works immediately
    setTimeout(swInitTiles, 500);
  });

})();
</script>
"""
```

- [ ] **Step 2: Commit**

```bash
git add ui_template.py
git commit -m "feat: add JS controller template for theme/layers/drawer"
```

---

## Task 5: Add `inject_ui()` to `read_plot.py` and wire it up

**Files:**
- Modify: `read_plot.py`

- [ ] **Step 1: Import `ui_template` at the top of `read_plot.py`**

Add after the existing imports (around line 16, after `from folium import IFrame`):

```python
import json
from ui_template import CSS_TEMPLATE, HTML_TEMPLATE, JS_TEMPLATE
```

- [ ] **Step 2: Add `inject_ui()` function after the imports block (before any global code)**

Add this function definition right after the imports:

```python
def inject_ui(raw_html: str, rides_data: list) -> str:
    """Post-process Folium's HTML output to inject the custom UI shell."""
    rides_json = json.dumps(rides_data, ensure_ascii=False, indent=None)
    rides_script = f"<script>const RIDES = {rides_json};</script>"

    # 1. Inject CSS into <head>
    raw_html = raw_html.replace("</head>", CSS_TEMPLATE + "\n</head>", 1)

    # 2. Inject RIDES data + HTML shell + JS before </body>
    injection = f"\n{rides_script}\n{HTML_TEMPLATE}\n{JS_TEMPLATE}\n"
    raw_html = raw_html.replace("</body>", injection + "</body>", 1)

    return raw_html
```

- [ ] **Step 3: Replace `map1.save("index.html")` at the bottom of `read_plot.py`**

Find the last two lines:
```python
map1.add_child(folium.LayerControl())
map1.save("index.html")
```

Replace with:
```python
map1.add_child(folium.LayerControl())
raw_html = map1._repr_html_()
final_html = inject_ui(raw_html, ride_data_list)
with open("index.html", "w", encoding="utf-8") as f:
    f.write(final_html)
print("index.html written.")
```

- [ ] **Step 4: Run the full script to regenerate `index.html`**

```bash
cd /Users/harshad/P_Projects/swiss_tour_visualization
source env/bin/activate
python read_plot.py
```

Expected: script runs to completion, prints `index.html written.`, file size is ~15MB+.

- [ ] **Step 5: Spot-check the output**

```bash
grep -c "sw-drawer\|sw-chips-bar\|RIDES =" index.html
```

Expected: output shows `3` (all three injected markers are present).

- [ ] **Step 6: Commit**

```bash
git add read_plot.py index.html
git commit -m "feat: inject custom UI shell into Folium output"
```

---

## Task 6: Wire polyline clicks to `swSelectRide()`

**Files:**
- Modify: `read_plot.py` (the `folium.PolyLine` section, around line 433)

Currently polylines open an IFrame popup. We replace this with a lightweight JS onclick that calls `swSelectRide()`.

- [ ] **Step 1: Remove the IFrame popup from the polyline and add an onclick**

Find this block (around line 405–438):
```python
    iframe = IFrame(full_html, width=720, height=400)
    popup = folium.Popup(iframe, max_width=750)
    ...
    folium.PolyLine(
        [(lat, lon) for lat, lon, *_ in ride_coords],
        color="purple",
        popup=popup,
        weight=3,
    ).add_to(rides)
```

Replace the `iframe`, `popup`, and `folium.PolyLine` call with:

```python
    polyline_js = f"swSelectRide({int(ride_id)})"
    folium.PolyLine(
        [(lat, lon) for lat, lon, *_ in ride_coords],
        color="var(--accent, #ff6b35)",
        weight=4,
        opacity=0.85,
        tooltip=f"Ride {ride_id} — {canton_label}",
    ).add_child(folium.ClickForMarker(popup=None)).add_to(rides)
```

Wait — `ClickForMarker` won't call our JS. Use a `folium.Popup` with a one-liner JS redirect instead:

```python
    click_popup = folium.Popup(
        f'<script>swSelectRide({int(ride_id)})</script>',
        max_width=1
    )
    folium.PolyLine(
        [(lat, lon) for lat, lon, *_ in ride_coords],
        color="#ff6b35",
        weight=4,
        opacity=0.85,
        tooltip=f"Ride {ride_id} · {canton_label}",
        popup=click_popup,
    ).add_to(rides)
```

Also remove the now-unused variables `encoded`, `img_html`, `full_html`, `iframe` from the ride loop (they were only used for the old popup). Remove these lines:

```python
    # Create popup with image or iframe
    encoded = base64.b64encode(open(elev_plot_path, "rb").read()).decode()
    img_html = f"""
    <div style="width: 100%; height: 100%;">
        <img src="data:image/png;base64,{encoded}" style="width: 100%; height: auto;">
    </div>
    """

    # Combined HTML
    full_html = f"""
    <div style="font-family: Arial; font-size:12px;">
    <img src="data:image/png;base64,{encoded}" width="700" height="270" style="margin:0; padding:0; display:inline;" />
    {stats_html}
    </div>
    """

    iframe = IFrame(full_html, width=720, height=400)
    popup = folium.Popup(iframe, max_width=750)
```

Also remove `from folium import IFrame` from the imports since it's no longer needed.

- [ ] **Step 2: Regenerate and verify**

```bash
cd /Users/harshad/P_Projects/swiss_tour_visualization
source env/bin/activate
python read_plot.py
```

Expected: completes without error, `index.html` written. File size will be significantly smaller (~3–5MB) since base64 PNG embeds are removed.

```bash
grep -c "swSelectRide" index.html
```

Expected: `10` (one per ride).

- [ ] **Step 3: Commit**

```bash
git add read_plot.py index.html
git commit -m "feat: replace IFrame popups with swSelectRide() onclick on polylines"
```

---

## Task 7: Manual browser test and polish

**Files:**
- Modify: `read_plot.py` and/or `ui_template.py` as needed for fixes found during testing.

- [ ] **Step 1: Open `index.html` in a browser**

```bash
open index.html
```

Verify these work:
1. Map fills the full viewport, no scroll bars.
2. Title badge visible top-left with theme toggle button.
3. Layer checkboxes visible top-right; unchecking "Rides" hides the polylines.
4. Chips bar visible at bottom with one chip per ride (10 chips).
5. Clicking a chip expands the drawer, shows ride name, elevation SVG, and stat row.
6. Clicking drawer handle collapses it; pressing Escape also collapses.
7. Theme toggle switches to light mode — map tiles, drawer background, accent color all change.
8. Clicking a polyline on the map triggers `swSelectRide()` and expands the drawer.

- [ ] **Step 2: Fix layer toggle if needed**

Folium's `FeatureGroup` stores the layer name in `layer.options.name`. If the `swToggleLayer` function doesn't find layers by name, add this debug snippet to the browser console:

```javascript
let m = null;
for (let k in window) {
  if (window[k] && window[k]._leaflet_id !== undefined) { m = window[k]; break; }
}
m.eachLayer(l => console.log(l.options && l.options.name));
```

If the names printed don't match exactly, update the `name` values in `HTML_TEMPLATE`'s checkboxes to match.

- [ ] **Step 3: Fix map bottom padding**

The drawer overlaps the bottom of the map when expanded. Add map padding adjustment to the JS. In `JS_TEMPLATE`, inside `swSelectRide`, add after `document.getElementById('sw-drawer').className = 'expanded';`:

```javascript
    if (_map) _map.invalidateSize();
```

And in `swCollapseDrawer`, add after setting `className = 'collapsed'`:

```javascript
    if (_map) _map.invalidateSize();
```

Update `ui_template.py` with these additions, then regenerate:

```bash
python read_plot.py
```

- [ ] **Step 4: Final commit**

```bash
git add read_plot.py ui_template.py index.html
git commit -m "feat: complete UI redesign — dark/light theme, drawer, ride chips"
```

---

## Self-Review

**Spec coverage check:**

| Spec requirement | Task |
|---|---|
| Full-screen map | Task 2 (CSS: `.folium-map { position:fixed; inset:0 }`) |
| Bottom drawer collapsed/expanded | Tasks 2, 3, 4 |
| Dark/light theme with CSS variables | Task 2 |
| Theme toggle persisted to localStorage | Task 4 |
| Tile layer swap on theme change | Task 4 |
| Title badge top-left | Task 3 |
| Layer control top-right (replaces Leaflet default) | Tasks 3, 4 |
| Ride chips scrollable | Tasks 3, 4 |
| Active chip highlighted | Task 4 |
| Hint text when no ride selected | Tasks 3, 4 |
| Expanded drawer: elevation SVG | Tasks 3, 4 |
| Expanded drawer: stat row (6 stats) | Tasks 3, 4 |
| Click polyline → select ride | Task 6 |
| Escape key collapses drawer | Task 4 |
| `const RIDES` data injection | Tasks 1, 5 |
| Remove base64 PNG embeds (smaller file) | Task 6 |
| `inject_ui()` post-processing approach | Task 5 |

All spec requirements covered. No gaps found.

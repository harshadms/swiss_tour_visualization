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

/* ── Leave Folium's map layout completely alone ─────────── */
/* Only the overlay elements are positioned by us */

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
#sw-hint {
  font-size:10px; color:var(--text2); font-style:italic;
  transition:opacity .3s;
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
#sw-drawer.collapsed { height:0; }
#sw-drawer.expanded  { height:280px; }

/* drawer content: column layout — handle / stats / chart */
#sw-drawer-content {
  height:280px;
  display:flex; flex-direction:column;
  overflow:hidden;
}

/* ── Handle row ── */
#sw-handle {
  display:flex; align-items:center; justify-content:space-between;
  padding:8px 16px 6px;
  border-bottom:1px solid var(--border);
  cursor:pointer; flex-shrink:0;
}
#sw-handle .ride-name {
  font-size:15px; font-weight:700; letter-spacing:.5px; color:var(--accent);
}
#sw-handle .ride-subtitle {
  font-size:11px; color:var(--text2); margin-top:2px;
  white-space:nowrap; overflow:hidden; text-overflow:ellipsis;
  max-width:70vw;
}
#sw-handle .chevron { color:var(--accent); font-size:13px; margin-left:12px; flex-shrink:0; }

/* ── Stat row (Komoot-style) ── */
#sw-stat-row {
  display:flex;
  border-bottom:1px solid var(--border);
  flex-shrink:0;
}
.sw-stat {
  flex:1; text-align:center;
  padding:7px 4px 6px;
  border-right:1px solid var(--border);
}
.sw-stat:last-child { border-right:none; }
.sw-stat .s-val  { font-size:16px; font-weight:700; color:var(--text1); line-height:1; }
.sw-stat .s-unit { font-size:10px; color:var(--text2); margin-left:1px; }
.sw-stat .s-lbl  { font-size:9px; color:var(--text2); text-transform:uppercase; letter-spacing:.4px; margin-top:3px; }

/* ── Elevation chart area ── */
#sw-chart-wrap {
  flex:1; position:relative; min-height:0;
  padding:4px 16px 8px 44px; /* left pad for Y-axis labels, bottom for baseline */
}
#sw-elev-svg {
  width:100%; height:100%;
  display:block; overflow:visible;
}
.sw-elev-line { fill:none; stroke:var(--accent); stroke-width:1.8; stroke-linecap:round; stroke-linejoin:round; }
.sw-elev-fill { fill:var(--accent); opacity:.15; stroke:none; }

/* Y-axis labels */
.sw-y-label {
  position:absolute; left:0; width:40px;
  font-size:11px; color:var(--text2);
  text-align:right; transform:translateY(-50%);
  pointer-events:none; line-height:1;
}
/* Y-axis grid lines rendered in SVG via JS */
.sw-grid-line { stroke:var(--border); stroke-width:0.5; stroke-dasharray:3,3; }

/* Town waypoint markers — now fully SVG, no DOM divs needed */
.sw-wp-line  { stroke:var(--border); stroke-width:0.8; stroke-dasharray:2,2; }
.sw-wp-dot   { fill:var(--accent); stroke:var(--bg); stroke-width:1.5; }
.sw-wp-label { fill:var(--text2); font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif; }

/* Hover cursor line on chart */
#sw-hover-line {
  stroke:var(--text1); stroke-width:1; stroke-dasharray:3,2;
  pointer-events:none; opacity:0;
}
/* Map hover dot */
#sw-map-dot {
  pointer-events:none;
}
</style>
"""

HTML_TEMPLATE = """
<!-- ── Title badge ───────────────────────── -->
<div id="sw-title">
  <span class="title-text">🚴 SWISS TOUR</span>
  <span id="sw-hint">click a ride</span>
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
  <div id="sw-drawer-content">

    <!-- handle: click to collapse -->
    <div id="sw-handle" onclick="swCollapseDrawer()">
      <div style="min-width:0">
        <div class="ride-name" id="sw-ride-name">—</div>
        <div class="ride-subtitle" id="sw-ride-subtitle">—</div>
      </div>
      <span class="chevron">↓</span>
    </div>

    <!-- stat row -->
    <div id="sw-stat-row">
      <div class="sw-stat"><div class="s-val" id="sw-s-dist">—</div><div class="s-lbl">Distance</div></div>
      <div class="sw-stat"><div class="s-val" id="sw-s-maxele">—</div><div class="s-lbl">Max Elev</div></div>
      <div class="sw-stat"><div class="s-val" id="sw-s-minele">—</div><div class="s-lbl">Min Elev</div></div>
      <div class="sw-stat"><div class="s-val" id="sw-s-avgspd">—</div><div class="s-lbl">Avg Speed</div></div>
      <div class="sw-stat"><div class="s-val" id="sw-s-maxspd">—</div><div class="s-lbl">Max Speed</div></div>
      <div class="sw-stat"><div class="s-val" id="sw-s-avggrad">—</div><div class="s-lbl">Avg Grade</div></div>
    </div>

    <!-- elevation chart -->
    <div id="sw-chart-wrap">
      <svg id="sw-elev-svg" viewBox="0 0 800 100" preserveAspectRatio="none">
        <g id="sw-grid-g"></g>
        <polygon class="sw-elev-fill" id="sw-elev-fill" points=""/>
        <polyline class="sw-elev-line" id="sw-elev-line" points=""/>
        <g id="sw-waypoint-lines-g"></g>
        <line id="sw-hover-line" x1="0" y1="0" x2="0" y2="100"/>
      </svg>
      <!-- Y-axis labels and waypoint markers injected by JS into #sw-chart-wrap -->
    </div>

  </div>
</div>
"""

JS_TEMPLATE = """
<script>
(function() {
  // ── Theme ──────────────────────────────────────────────────
  var _darkTile  = null;
  var _lightTile = null;
  var _map       = null;

  function swFindMap() {
    if (window._swMap) return window._swMap;
    for (var k in window) {
      try {
        if (window[k] && window[k]._leaflet_id !== undefined && typeof window[k].addLayer === 'function') {
          return window[k];
        }
      } catch(e) {}
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
    // Apply the saved theme's tile on initial load
    var theme = document.documentElement.getAttribute('data-theme') || 'dark';
    if (theme === 'dark') {
      _map.eachLayer(function(l) { if (l._url && l._url.indexOf('light_all') !== -1) _map.removeLayer(l); });
      _darkTile.addTo(_map);
    }
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
        _map.eachLayer(function(l) { if (l._url && l._url.indexOf('light_all') !== -1) _map.removeLayer(l); });
        _darkTile.addTo(_map);
      } else {
        _map.eachLayer(function(l) { if (l._url && l._url.indexOf('dark_all') !== -1) _map.removeLayer(l); });
        _lightTile.addTo(_map);
      }
    }
  };

  // Apply saved theme on load
  var savedTheme = localStorage.getItem('sw-theme') || 'dark';
  document.documentElement.setAttribute('data-theme', savedTheme);

  // ── Layer control ──────────────────────────────────────────
  window.swToggleLayer = function(name, visible) {
    if (!_map) _map = swFindMap();
    if (!_map) return;
    for (var k in window) {
      try {
        var obj = window[k];
        if (obj && obj._leaflet_id !== undefined && obj.options && obj.options.name === name) {
          if (visible) { _map.addLayer(obj); } else { _map.removeLayer(obj); }
        }
      } catch(e) {}
    }
  };

  // Hide Leaflet's default layer control
  function swHideDefaultControls() {
    var lc = document.querySelector('.leaflet-control-layers');
    if (lc) lc.style.display = 'none';
    var btn = document.getElementById('sw-theme-toggle');
    if (btn && savedTheme === 'light') btn.textContent = '🌙 Dark';
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', swHideDefaultControls);
  } else {
    swHideDefaultControls();
  }

  // ── Ride selection ─────────────────────────────────────────
  var _activeRideId = null;

  window.swSelectRide = function(rideId) {
    if (typeof RIDES === 'undefined') return;
    var ride = RIDES.find(function(r) { return r.id === rideId; });
    if (!ride) return;

    _activeRideId = rideId;

    // Hide hint
    var hint = document.getElementById('sw-hint');
    if (hint) hint.style.opacity = '0';

    // Handle
    document.getElementById('sw-ride-name').textContent =
      ride.label.toUpperCase() + ' \\u2014 ' + ride.canton.toUpperCase();
    document.getElementById('sw-ride-subtitle').textContent =
      ride.towns.join(' \\u2192 ');

    // Stats
    function fmt(v, unit) {
      return v != null
        ? v + '<span class="s-unit"> ' + unit + '</span>'
        : '\\u2014';
    }
    document.getElementById('sw-s-dist').innerHTML    = fmt(ride.stats.distance_km,      'km');
    document.getElementById('sw-s-maxele').innerHTML  = fmt(ride.stats.max_elevation_m,  'm');
    document.getElementById('sw-s-minele').innerHTML  = fmt(ride.stats.min_elevation_m,  'm');
    document.getElementById('sw-s-avgspd').innerHTML  = fmt(ride.stats.avg_speed_kmh,    'km/h');
    document.getElementById('sw-s-maxspd').innerHTML  = fmt(ride.stats.max_speed_kmh,    'km/h');
    document.getElementById('sw-s-avggrad').innerHTML = fmt(ride.stats.avg_gradient_pct, '%');

    // Draw elevation chart with Y-axis and waypoints
    swDrawElevation(ride);
    swSetupHover(ride);

    // Expand drawer
    document.getElementById('sw-drawer').className = 'expanded';
    if (_map) _map.invalidateSize();
  };

  function swDrawElevation(ride) {
    var dists = ride.distances_km;
    var eles  = ride.elevations_m;
    if (!dists || !dists.length) return;

    // SVG coordinate space — use actual pixel width so text doesn't stretch
    // Top pad: 22 (waypoint labels + dot), bottom pad: 18 (baseline gap like Komoot)
    var svg = document.getElementById('sw-elev-svg');
    var W = Math.round(svg.getBoundingClientRect().width) || 800;
    var H = 120;
    var PAD_TOP = 22, PAD_BOT = 18;
    var chartH = H - PAD_TOP - PAD_BOT; // usable height for elevation

    var minE = Math.min.apply(null, eles);
    var maxE = Math.max.apply(null, eles);
    var rangeE = maxE - minE || 1;
    var maxD = dists[dists.length - 1] || 1;

    // Helper: elevation value → SVG y
    function eToY(e) {
      return PAD_TOP + chartH - ((e - minE) / rangeE) * chartH;
    }

    // ── Elevation line + fill ──
    var pts = [];
    for (var i = 0; i < dists.length; i++) {
      var x = (dists[i] / maxD) * W;
      var y = eToY(eles[i]);
      pts.push(x.toFixed(1) + ',' + y.toFixed(1));
    }
    var lineStr = pts.join(' ');
    // Fill closes down to baseline (PAD_BOT above bottom), not the very bottom edge
    var baseY = H - PAD_BOT;
    var fillStr = lineStr + ' ' + W + ',' + baseY + ' 0,' + baseY;
    document.getElementById('sw-elev-line').setAttribute('points', lineStr);
    document.getElementById('sw-elev-fill').setAttribute('points', fillStr);

    // Update SVG viewBox to match new H
    document.getElementById('sw-elev-svg').setAttribute('viewBox', '0 0 ' + W + ' ' + H);
    document.getElementById('sw-hover-line').setAttribute('y2', H);

    // ── Y-axis grid lines + labels ──
    var wrap = document.getElementById('sw-chart-wrap');
    var oldLabels = wrap.querySelectorAll('.sw-y-label');
    for (var j = 0; j < oldLabels.length; j++) oldLabels[j].parentNode.removeChild(oldLabels[j]);

    var gridG = document.getElementById('sw-grid-g');
    while (gridG.firstChild) gridG.removeChild(gridG.firstChild);

    // Pick ~4 nice round elevation ticks
    var rawStep = rangeE / 4;
    var magnitude = Math.pow(10, Math.floor(Math.log(rawStep) / Math.LN10));
    var niceStep = Math.ceil(rawStep / magnitude) * magnitude;
    var tickStart = Math.ceil(minE / niceStep) * niceStep;

    for (var e = tickStart; e <= maxE; e += niceStep) {
      var yPx = eToY(e);

      // Grid line in SVG
      var gline = document.createElementNS('http://www.w3.org/2000/svg', 'line');
      gline.setAttribute('class', 'sw-grid-line');
      gline.setAttribute('x1', '0'); gline.setAttribute('x2', W);
      gline.setAttribute('y1', yPx.toFixed(1)); gline.setAttribute('y2', yPx.toFixed(1));
      gridG.appendChild(gline);

      // Label pinned to left of chart-wrap using %
      // Map SVG y → % of chart-wrap clientHeight
      // We can't know exact px from JS easily so use a ratio: yPx/H as fraction of wrap height
      var lbl = document.createElement('div');
      lbl.className = 'sw-y-label';
      lbl.textContent = Math.round(e) + 'm';
      lbl.style.top = ((yPx / H) * 100).toFixed(1) + '%';
      wrap.appendChild(lbl);
    }

    // ── Town waypoint markers (fully in SVG so they scale) ──
    var wpG = document.getElementById('sw-waypoint-lines-g');
    while (wpG.firstChild) wpG.removeChild(wpG.firstChild);

    var towns     = ride.towns;
    var townDists = ride.town_distances_km;

    if (townDists && townDists.length === towns.length) {
      for (var t = 0; t < towns.length; t++) {
        var frac   = townDists[t] / maxD;
        var wpX    = frac * W;

        // Find the elevation at this town to place the dot on the line
        // Binary search in dists for nearest index
        var lo = 0, hi = dists.length - 1;
        while (lo < hi) { var mid = (lo + hi) >> 1; if (dists[mid] < townDists[t]) lo = mid + 1; else hi = mid; }
        var dotY = eToY(eles[lo]);

        // Dashed vertical line from top pad down to baseline
        var vline = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        vline.setAttribute('class', 'sw-wp-line');
        vline.setAttribute('x1', wpX.toFixed(1)); vline.setAttribute('x2', wpX.toFixed(1));
        vline.setAttribute('y1', PAD_TOP.toFixed(1)); vline.setAttribute('y2', (H - PAD_BOT).toFixed(1));
        wpG.appendChild(vline);

        // Dot on the elevation line
        var dot = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
        dot.setAttribute('class', 'sw-wp-dot');
        dot.setAttribute('cx', wpX.toFixed(1));
        dot.setAttribute('cy', dotY.toFixed(1));
        dot.setAttribute('r', '4');
        wpG.appendChild(dot);

        // Town name label at top, inside PAD_TOP area
        // Use SVG text with a small background rect for legibility
        var labelY = PAD_TOP - 6; // above the chart area

        var rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        var approxW = towns[t].length * 5.5 + 6;
        var rectX = Math.min(Math.max(wpX - approxW / 2, 2), W - approxW - 2);
        rect.setAttribute('x', rectX.toFixed(1));
        rect.setAttribute('y', (labelY - 9).toFixed(1));
        rect.setAttribute('width', approxW.toFixed(1));
        rect.setAttribute('height', '11');
        rect.setAttribute('rx', '2');
        rect.setAttribute('fill', 'var(--bg2)');
        rect.setAttribute('opacity', '0.85');
        wpG.appendChild(rect);

        var txt = document.createElementNS('http://www.w3.org/2000/svg', 'text');
        txt.setAttribute('class', 'sw-wp-label');
        txt.setAttribute('x', (rectX + approxW / 2).toFixed(1));
        txt.setAttribute('y', (labelY).toFixed(1));
        txt.setAttribute('text-anchor', 'middle');
        txt.setAttribute('font-size', '8');
        txt.textContent = towns[t];
        wpG.appendChild(txt);
      }
    }
  }

  // ── Elevation hover: cursor line + map dot ────────────────
  var _hoverMarker = null;
  var _hoverRide   = null;

  function swSetupHover(ride) {
    _hoverRide = ride;
    var svg = document.getElementById('sw-elev-svg');
    var hline = document.getElementById('sw-hover-line');

    svg.addEventListener('mousemove', function(e) {
      if (!_hoverRide) return;
      var rect = svg.getBoundingClientRect();
      var frac = (e.clientX - rect.left) / rect.width;
      frac = Math.max(0, Math.min(1, frac));

      // Show cursor line at SVG x coordinate (viewBox matches actual px width)
      var vb = svg.viewBox.baseVal;
      var svgX = (frac * (vb.width || 800)).toFixed(1);
      hline.setAttribute('x1', svgX);
      hline.setAttribute('x2', svgX);
      hline.style.opacity = '1';

      // Find nearest sample by distance fraction
      var dists = _hoverRide.distances_km;
      var maxD  = dists[dists.length - 1] || 1;
      var target = frac * maxD;
      // Binary search
      var lo = 0, hi = dists.length - 1;
      while (lo < hi) {
        var mid = (lo + hi) >> 1;
        if (dists[mid] < target) lo = mid + 1; else hi = mid;
      }
      var coords = _hoverRide.coords;
      if (!coords || !coords[lo]) return;
      var latlng = L.latLng(coords[lo][0], coords[lo][1]);

      if (!_hoverMarker) {
        _hoverMarker = L.circleMarker(latlng, {
          radius: 7, color: '#ffffff', weight: 2,
          fillColor: '#ff6b35', fillOpacity: 1,
          pane: 'markerPane'
        }).addTo(_map);
      } else {
        _hoverMarker.setLatLng(latlng);
      }
    });

    svg.addEventListener('mouseleave', function() {
      var hline = document.getElementById('sw-hover-line');
      hline.style.opacity = '0';
      if (_hoverMarker && _map) {
        _map.removeLayer(_hoverMarker);
        _hoverMarker = null;
      }
    });
  }

  window.swCollapseDrawer = function() {
    document.getElementById('sw-drawer').className = 'collapsed';
    _activeRideId = null;
    _hoverRide = null;
    if (_hoverMarker && _map) { _map.removeLayer(_hoverMarker); _hoverMarker = null; }
    var hline = document.getElementById('sw-hover-line');
    if (hline) hline.style.opacity = '0';
    var hint = document.getElementById('sw-hint');
    if (hint) hint.style.opacity = '1';
    if (_map) _map.invalidateSize();
  };

  // Escape key collapses drawer
  document.addEventListener('keydown', function(e) {
    if (e.key === 'Escape') swCollapseDrawer();
  });

  // Init
  function swInit() {
    setTimeout(swInitTiles, 500);
  }
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', swInit);
  } else {
    swInit();
  }

})();
</script>
"""

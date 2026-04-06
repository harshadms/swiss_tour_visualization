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

JS_TEMPLATE = """
<script>
(function() {
  // ── Theme ──────────────────────────────────────────────────
  var _darkTile  = null;
  var _lightTile = null;
  var _map       = null;

  function swFindMap() {
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
    // Walk all global vars to find Folium feature groups by name
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
  document.addEventListener('DOMContentLoaded', function() {
    var lc = document.querySelector('.leaflet-control-layers');
    if (lc) lc.style.display = 'none';
    // Apply saved theme button label
    var btn = document.getElementById('sw-theme-toggle');
    if (btn && savedTheme === 'light') btn.textContent = '🌙 Dark';
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
      ride.towns.length > 1
        ? ride.towns[0] + ' \u2192 ' + ride.towns[ride.towns.length - 1]
        : (ride.towns[0] || '');
    var dist = ride.stats.distance_km;
    document.getElementById('sw-distance').textContent = dist ? dist + ' km' : '\u2014';

    // Populate stat row
    function fmt(v, unit) {
      return v != null ? v + '<span class="s-unit"> ' + unit + '</span>' : '\u2014';
    }
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
    if (_map) _map.invalidateSize();
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
    if (_map) _map.invalidateSize();
  };

  // Escape key collapses drawer
  document.addEventListener('keydown', function(e) {
    if (e.key === 'Escape') swCollapseDrawer();
  });

  // Init on DOM ready
  document.addEventListener('DOMContentLoaded', function() {
    swBuildChips();
    setTimeout(swInitTiles, 500);
  });

})();
</script>
"""

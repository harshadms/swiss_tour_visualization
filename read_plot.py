import json
import time
import warnings
from geopy.distance import geodesic

warnings.simplefilter(action="ignore", category=FutureWarning)

import geopandas as gpd
import pandas as pd
import folium
import gpxpy
import traceback
import requests

from ui_template import CSS_TEMPLATE, HTML_TEMPLATE, JS_TEMPLATE


def inject_ui(raw_html: str, rides_data: list) -> str:
    """Post-process Folium's HTML output to inject the custom UI shell."""
    import re
    rides_json = json.dumps(rides_data, ensure_ascii=False)
    rides_script = f"<script>const RIDES = {rides_json};</script>"

    # Extract the Leaflet map variable name so we can call invalidateSize()
    # Folium generates something like: var map_abc123 = L.map(...)
    map_var_match = re.search(r'var (map_[a-f0-9]+) = L\.map\(', raw_html)
    map_var = map_var_match.group(1) if map_var_match else None

    # Build a small init snippet that fixes map size after the drawer is rendered
    map_fix = ""
    if map_var:
        map_fix = f"""
<script>
(function waitForMap() {{
  if (typeof {map_var} === 'undefined') {{ setTimeout(waitForMap, 50); return; }}
  window._swMap = {map_var};
  {map_var}.invalidateSize();
  // Wire polyline clicks: tooltip is "ride-id-N"
  {map_var}.eachLayer(function(layer) {{
    if (layer.getTooltip) {{
      var tt = layer.getTooltip();
      if (tt) {{
        var content = tt.getContent ? tt.getContent() : tt._content;
        var m = content && content.match(/^ride-id-(\\d+)$/);
        if (m) {{
          var rid = parseInt(m[1]);
          layer.setStyle && layer.setStyle({{color: '#ff6b35', weight: 4, opacity: 0.85}});
          layer.on('click', function() {{ swSelectRide(rid); }});
        }}
      }}
    }}
  }});
}})();
</script>"""

    raw_html = raw_html.replace("</head>", CSS_TEMPLATE + "\n</head>", 1)
    # Inject before </html> — Folium puts its JS *after* </body>, so we must
    # go after all of Folium's code to ensure the map variable exists.
    injection = f"\n{rides_script}\n{HTML_TEMPLATE}\n{JS_TEMPLATE}\n{map_fix}\n"
    raw_html = raw_html.replace("</html>", injection + "</html>", 1)
    return raw_html


def read_gpx(ride_name):
    gpx_file = open(ride_name, "r")
    gpx = gpxpy.parse(gpx_file)
    coords = []
    ele = []
    time = []
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                coords.append(
                    (point.latitude, point.longitude)
                )  # -> {point.elevation}')
                ele.append(point.elevation)
                time.append(point.time)
    return coords, ele, time


def listPoints(someGeometry):
    """List the points in a Polygon in a geometry entry - some polygons are more complex than others, so accommodating for that"""
    pointList = []
    try:
        # Note: might miss parts within parts with this
        for part in someGeometry:
            x, y = part.exterior.coords.xy
            pointList.append(list(zip(x, y)))
    except:
        try:
            x, y = someGeometry.exterior.coords.xy
            pointList.append(list(zip(x, y)))
        except:
            # this will return the geometry as is, enabling you to see if special handling is required - then modify the function as need be
            pointList.append(someGeometry)
    return pointList


def convert_to_latlon(coords):
    latlon = []
    for c in coords:
        latlon.append((c[1], c[0]))

    return latlon


def plot_polyline(geom_object, color, layer, fillcolor, popup):
    if "multipolygon.MultiPolygon" in str(type(geom_object)):
        for i in geom_object.geoms:
            coords = convert_to_latlon(list(i.exterior.coords))
            folium.Polygon(
                coords,
                dash_array="5",
                fill_color=fillcolor,
                color=color,
                popup=popup,
                weight=1.5,
            ).add_to(layer)
    else:
        coords = convert_to_latlon(list(geom_object.exterior.coords))
        folium.Polygon(
            coords,
            dash_array="5",
            fill_color=fillcolor,
            color=color,
            popup=popup,
            weight=1.5,
        ).add_to(layer)


def get_station_coordinates(city):
    df = pd.read_csv("./coords_from_api.csv")
    record = df.query(f"town == '{city}'")
    if len(record) >= 1:
        return record.lat.values[0], record.lon.values[0]

    query = f"{city} Bahnhof, Switzerland"
    url = "https://nominatim.openstreetmap.org/search"
    params = {"q": query, "format": "json"}
    headers = {"User-Agent": "PersonalApp/1.0 (usernamesdontmatter@proton.me)"}
    response = requests.get(url, headers=headers, params=params)

    if response.status_code != 200:
        print(f"Received invalid response - {response.status_code}")
        print(response.text)
        raise Exception

    data = response.json()
    if data:
        with open("./coords_from_api.csv", "a") as f:
            f.write(f"{city}, {data[0]['lat']}, {data[0]['lon']}\n")

        return float(data[0]["lat"]), float(data[0]["lon"])

    return None


def generate_elevation_plot(coords, elevations, ride_id):
    import matplotlib.pyplot as plt
    import os

    os.makedirs("./elevation_plots", exist_ok=True)

    distances = [0.0]
    for i in range(1, len(coords)):
        d = geodesic(coords[i - 1], coords[i]).km
        distances.append(distances[-1] + d)

    # Beautify the plot
    plt.figure(figsize=(8, 3))
    plt.plot(distances, elevations, color="#6A0DAD", linewidth=2)
    plt.fill_between(distances, elevations, color="#D8BFD8", alpha=0.4)

    lower_limit = int(round((min(elevations) - 50) / 50.0) * 50)
    upper_limit = int(round((max(elevations) + 50) / 50.0) * 50)

    plt.ylim(lower_limit, upper_limit)
    plt.title(f"Elevation Profile – Ride {ride_id}", fontsize=14, fontweight="bold")
    plt.xlabel("Distance (km)", fontsize=12)
    plt.ylabel("Elevation (m)", fontsize=12)
    plt.grid(True, linestyle="--", alpha=0.5)

    plt.tight_layout()
    path = f"./elevation_plots/ride_{ride_id}.png"
    plt.savefig(path, dpi=100, bbox_inches="tight", pad_inches=0.05)
    plt.close()

    return path


def moving_average(values, window_size=5):
    import numpy as np

    if len(values) < window_size:
        return values
    return np.convolve(values, np.ones(window_size) / window_size, mode="valid")


def generate_stats_from_gps(coords, elevations, timestamps):
    raw_speeds = []
    raw_gradients = []
    import numpy as np

    window_size = 5

    min_dist = 2
    min_time = 1

    for i in range(1, len(coords)):
        p1, p2 = coords[i - 1], coords[i]

        lat1, lon1 = p1
        lat2, lon2 = p2

        ele1 = elevations[i - 1]
        ele2 = elevations[i]

        t1 = timestamps[i - 1]
        t2 = timestamps[i]

        # Skip if missing data
        if None in (lat1, lon1, ele1, t1, lat2, lon2, ele2, t2):
            continue

        dist_m = geodesic((lat1, lon1), (lat2, lon2)).meters
        time_s = (t2 - t1).total_seconds()

        if dist_m < min_dist or time_s < min_time:
            continue

        speed_kmh = (dist_m / time_s) * 3.6
        gradient = (ele2 - ele1) / dist_m * 100

        # Filter extreme gradient noise
        if abs(gradient) < 100:
            raw_gradients.append(gradient)

        raw_speeds.append(speed_kmh)
        elevations.extend([ele1, ele2])

    # Apply moving average smoothing
    smooth_speeds = moving_average(raw_speeds, window_size)
    smooth_gradients = moving_average(raw_gradients, window_size)

    stats = {
        "Max Elevation": round(max(elevations), 1) if elevations else None,
        "Min Elevation": round(min(elevations), 1) if elevations else None,
        "Max Gradient (%)": (
            round(max(smooth_gradients), 2) if len(smooth_gradients) > 0 else None
        ),
        "Min Gradient (%)": (
            round(min(smooth_gradients), 2) if len(smooth_gradients) > 0 else None
        ),
        "Avg Gradient (%)": (
            round(np.mean(smooth_gradients), 2) if len(smooth_gradients) > 0 else None
        ),
        "Max Speed (km/h)": (
            round(max(smooth_speeds), 1) if len(smooth_speeds) > 0 else None
        ),
        "Min Speed (km/h)": (
            round(min(smooth_speeds), 1) if len(smooth_speeds) > 0 else None
        ),
        "Avg Speed (km/h)": (
            round(np.mean(smooth_speeds), 1) if len(smooth_speeds) > 0 else None
        ),
    }
    
    # Use emojis or Font Awesome/Material if embedding full HTML later
    icons = {
        "Max Elevation": '<i class="fas fa-mountain"></i>',
        "Min Elevation": '<i class="fas fa-tree"></i>',
        "Max Gradient (%)": '<i class="fas fa-arrow-up"></i>',
        "Min Gradient (%)": '<i class="fas fa-arrow-down"></i>',
        "Avg Gradient (%)": '<i class="fas fa-chart-line"></i>',
        "Max Speed (km/h)": '<i class="fas fa-tachometer-alt"></i>',
        "Min Speed (km/h)": '<i class="fas fa-walking"></i>',
        "Avg Speed (km/h)": '<i class="fas fa-bicycle"></i>',
    }

    elevation_rows = "".join(
        f"<tr><td>{icons.get(k, '')} {k}</td><td style='text-align:right'>{v}</td></tr>"
        for k, v in stats.items() if "Elevation" in k
    )

    gradient_rows = "".join(
        f"<tr><td>{icons.get(k, '')} {k}</td><td style='text-align:right'>{v}</td></tr>"
        for k, v in stats.items() if "Gradient" in k
    )

    speed_rows = "".join(
        f"<tr><td>{icons.get(k, '')} {k}</td><td style='text-align:right'>{v}</td></tr>"
        for k, v in stats.items() if "Speed" in k
    )

    html = f"""
    <div style="font-family:Arial, sans-serif; font-size:13px; max-width:700px;">
      <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" rel="stylesheet">
      <table style="border-collapse:collapse; width:100%; margin-top:4px;">
        <thead>
          <tr><th colspan="2" style="text-align:left; padding:6px; border-bottom:1px solid #ccc; background:#f8f8f8;">🏔️ Elevation</th></tr>
        </thead>
        <tbody>{elevation_rows}</tbody>
        
        <thead>
          <tr><th colspan="2" style="text-align:left; padding:6px; border-bottom:1px solid #ccc; background:#f8f8f8;">📉 Gradient</th></tr>
        </thead>
        <tbody>{gradient_rows}</tbody>

        <thead>
          <tr><th colspan="2" style="text-align:left; padding:6px; border-bottom:1px solid #ccc; background:#f8f8f8;">🚴 Speed</th></tr>
        </thead>
        <tbody>{speed_rows}</tbody>
      </table>
    </div>
    """
    return html


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


# All towns
towns = pd.read_csv("rides/towns.csv", sep=",")

map1 = folium.Map(location=[46.8182, 8.2275], tiles="CartoDB Positron", zoom_start=8.5)

# Plot kantons
shapefile = gpd.read_file(
    "./shapefiles/swissBOUNDARIES3D_1_5_TLM_KANTONSGEBIET.shp", engine="pyogrio"
)
shapefile.to_crs(crs=4326, inplace=True)

# latlon = shapefile.query("NAME == 'Schweiz'").geometry.apply(lambda x: listPoints(x)).values.tolist()
completed_cantons = ["Zürich", "Zug", "Aargau", "Graubünden"]
cantons = folium.FeatureGroup(name="Cantons")
completed = folium.FeatureGroup(name="Completed")
districts = folium.FeatureGroup(name="Districts")
rides = folium.FeatureGroup(name="Rides")
stops = folium.FeatureGroup(name="Stops")

districts_shp = gpd.read_file(
    "./shapefiles/swissBOUNDARIES3D_1_5_TLM_BEZIRKSGEBIET.shp", engine="pyogrio"
)
districts_shp.to_crs(crs=4326, inplace=True)

kantons = [
    "Zürich",
    "Bern",
    "Luzern",
    "Uri",
    "Schwyz",
    "Obwalden",
    "Nidwalden",
    "Glarus",
    "Zug",
    "Fribourg",
    "Solothurn",
    "Basel-Stadt",
    "Basel-Landschaft",
    "Schaffhausen",
    "Appenzell",
    "Appenzell",
    "St.",
    "Graubünden",
    "Aargau",
    "Thurgau",
    "Ticino",
    "Vaud",
    "Valais",
    "Neuchâtel",
    "Genève",
    "Jura",
]

style1 = {"fillColor": "#e6e6e6", "color": "#003399"}
# Plot rides
rides_info = pd.read_csv("./rides/rides_info.csv")

bezirke_complete_border = "#e65c00"
bezirke_complete_fill = "#ffe0cc"

kanton_incomplete_border = "#1e7802"
kanton_incomplete_fill = "#dff7d7"

kanton_not_started_border = "#003399"
kanton_not_started_fill = "#e6e6e6"

for row in range(len(shapefile)):
    kanton_obj = shapefile.iloc[row]
    kanton_name = kanton_obj.NAME.strip()

    if kanton_name in completed_cantons:
        color = bezirke_complete_border
        fillColor = bezirke_complete_fill
        layer = completed
        kanton_num = kantons.index(kanton_name) + 1
        dist_obj = districts_shp.query(f"KANTONSNUM == {kanton_num}")

        if dist_obj.empty:
            print(f"No districts found for {kanton_name}, using canton geometry.")
            plot_polyline(kanton_obj.geometry, color, cantons, fillColor, kanton_name)
            continue

        for i in range(len(dist_obj.geometry)):
            dist_name = dist_obj.iloc[i].NAME
            if dist_name not in towns.bezirk.values:
                print(f"Bezirke not completed: {dist_name}")
                color = kanton_incomplete_border
                fillColor = kanton_incomplete_fill
            else:
                color = bezirke_complete_border
                fillColor = bezirke_complete_fill

            popup = dist_name
            tooltip = "Click me!"
            plot_polyline(dist_obj.iloc[i].geometry, color, districts, fillColor, popup)
    else:
        layer = cantons
        color = kanton_not_started_border
        fillColor = kanton_not_started_fill
        popup = kanton_name
        plot_polyline(kanton_obj.geometry, color, layer, fillColor, popup)

ride_data_list = []

for ride_id in rides_info.id.unique():
    ride = rides_info.query(f"id == {ride_id}")
    print(ride_id)
    ride_name = f"./rides/ride{ride_id}.xml"
    ride_coords, ride_ele, ride_timestamps = read_gpx(ride_name)
    ride_towns = ride.towns.values

    # Generate elevation plot
    elev_plot_path = generate_elevation_plot(ride_coords, ride_ele, ride_id)
    stats_html = generate_stats_from_gps(ride_coords, ride_ele, ride_timestamps)

    # Build structured data for UI injection
    ride_towns_list = list(ride.towns.values)
    ride_canton = rides_info.query(f"id == {ride_id}").bezirk.values[0]
    canton_label = next(
        (c for c in completed_cantons if c in ride_canton or ride_canton in c),
        ride_canton
    )

    ele_step = max(1, len(ride_ele) // 200)
    ele_sampled = ride_ele[::ele_step]
    coords_sampled = ride_coords[::ele_step]
    distances_sampled = [0.0]
    for i in range(1, len(coords_sampled)):
        d = geodesic(coords_sampled[i - 1], coords_sampled[i]).km
        distances_sampled.append(distances_sampled[-1] + d)

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
    for c in ride_towns:
        town = towns.query(f"town in '{c}'")
        try:
            try:
                st_coord = get_station_coordinates(c)
                time.sleep(1.2)
                if isinstance(st_coord, type(None)):
                    raise Exception
            except Exception as e:
                # traceback.print_exc()
                print(f"Picking from the list: {e} {c}")

                st_coord = [
                    town.lat.values[0],
                    town.lon.values[0],
                ]  # .#town.Geopos.values[0].split(",")
            folium.CircleMarker(
                st_coord, popup=c, fill=True, fill_opacity=1, radius=4
            ).add_to(stops)
        except:
            traceback.print_exc()
            print(c)
            pass

    # Add polyline — ride_id stored in tooltip for JS click wiring
    folium.PolyLine(
        [(lat, lon) for lat, lon in ride_coords],
        color="#ff6b35",
        weight=4,
        opacity=0.85,
        tooltip=f"ride-id-{int(ride_id)}",
    ).add_to(rides)


# Plot switzerland border
shapefile = gpd.read_file(
    "./shapefiles/swissBOUNDARIES3D_1_5_TLM_LANDESGEBIET.shp", engine="pyogrio"
)
shapefile.to_crs(crs=4326, inplace=True)

style1 = {"fillColor": "#FFFFFF", "color": "#006600", "opacity": 1, "weight": 0}
border_layer = folium.FeatureGroup(name="Schweiz")
folium.GeoJson(
    shapefile.query("NAME == 'Schweiz'").geometry, style_function=lambda x: style1
).add_to(border_layer)

map1.add_child(border_layer)
map1.add_child(cantons)
map1.add_child(completed)
map1.add_child(districts)
map1.add_child(rides)
map1.add_child(stops)


map1.add_child(folium.LayerControl())
raw_html = map1.get_root().render()
final_html = inject_ui(raw_html, ride_data_list)
with open("index.html", "w", encoding="utf-8") as f:
    f.write(final_html)
print("index.html written.")

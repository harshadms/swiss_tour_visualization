import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import geopandas as gpd
import pandas as pd
import folium
import gpxpy

def read_gpx(ride_name):
    gpx_file = open(ride_name, 'r')
    gpx = gpxpy.parse(gpx_file)
    coords = []
    for track in gpx.tracks:
        for segment in track.segments:
            for point in segment.points:
                coords.append((point.latitude,point.longitude))# -> {point.elevation}')

    return coords

def listPoints(someGeometry):
    '''List the points in a Polygon in a geometry entry - some polygons are more complex than others, so accommodating for that'''    
    pointList = []
    try:
        #Note: might miss parts within parts with this
        for part in someGeometry:
            x, y = part.exterior.coords.xy
            pointList.append(list(zip(x,y)))
    except:
        try:
            x,y = someGeometry.exterior.coords.xy
            pointList.append(list(zip(x,y)))
        except:
            #this will return the geometry as is, enabling you to see if special handling is required - then modify the function as need be
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
            folium.Polygon(coords, dash_array='5', fill_color=fillcolor, color=color, popup=popup, weight=1.5).add_to(layer)
    else:
        coords = convert_to_latlon(list(geom_object.exterior.coords))
        folium.Polygon(coords, dash_array='5', fill_color=fillcolor, color=color, popup=popup, weight=1.5).add_to(layer)

# All towns
towns = pd.read_csv("rides\haltestelle-offnungszeiten.csv", sep=';')

map1 = folium.Map(location=[46.8182, 8.2275], tiles="CartoDB Positron", zoom_start=8.5)

# Plot kantons
shapefile = gpd.read_file("swissBOUNDARIES3D_1_5_TLM_KANTONSGEBIET.shp", engine='pyogrio')
shapefile.to_crs(crs=4326, inplace=True)

# latlon = shapefile.query("NAME == 'Schweiz'").geometry.apply(lambda x: listPoints(x)).values.tolist()
completed_cantons = ["Zürich"]
cantons = folium.FeatureGroup(name="Cantons")
completed = folium.FeatureGroup(name="Completed")
districts = folium.FeatureGroup(name='Districts')
rides = folium.FeatureGroup(name='Rides')
stops = folium.FeatureGroup(name='Stops')

districts_shp = gpd.read_file("swissBOUNDARIES3D_1_5_TLM_BEZIRKSGEBIET.shp", engine='pyogrio')
districts_shp.to_crs(crs=4326, inplace=True)

style1 = {'fillColor': '#e6e6e6', 'color': '#003399'}

for row in range(len(shapefile)):
    kanton_obj = shapefile.iloc[row]

    if kanton_obj.NAME in completed_cantons:
        color = '#e65c00'
        fillColor='#ffe0cc'
        #layer = completed
        dist_obj = districts_shp.query(f"KANTONSNUM == 1")
    
        for i in range(len(dist_obj.geometry)):
            popup=dist_obj.iloc[i].NAME,
            tooltip="Click me!",
            #coords = convert_to_latlon(list(geom_object.iloc[i].exterior.coords))
            plot_polyline (dist_obj.iloc[i].geometry, color, districts, fillColor, popup[0])
    else:
        layer = cantons
        color = '#003399'
        fillColor = '#e6e6e6'
        popup=kanton_obj.NAME,
        plot_polyline(kanton_obj.geometry, color, layer, fillColor, popup[0])

# Plot rides
rides_info = pd.read_csv("./rides/rides_info.csv")
for i in range(len(rides_info)):
    ride = rides_info.iloc[i]
    
    ride_name = f"./rides/ride{ride.id}.xml"
    ride_coords = read_gpx(ride_name)
    ride_towns = ride.towns

    for c in ride_towns.split("-"):
        town = towns.query(f"stopname in '{c}'")
        try:
            print (f"{c} - {town.Geopos.values[0]}")
            st_coord = town.Geopos.values[0].split(",")
            folium.CircleMarker(st_coord, popup=c, fill=True, fill_opacity=1, radius=4).add_to(stops)
        except:
            pass

    folium.PolyLine(ride_coords, color='purple', popup=popup, weight=3).add_to(rides)

# Plot switzerland border
shapefile = gpd.read_file("swissBOUNDARIES3D_1_5_TLM_LANDESGEBIET.shp", engine='pyogrio')
shapefile.to_crs(crs=4326, inplace=True)

style1 = {'fillColor': '#FFFFFF', 'color': '#006600', 'opacity':1, 'weight':0}
border_layer = folium.FeatureGroup(name="Schweiz")
folium.GeoJson(shapefile.query("NAME == 'Schweiz'").geometry, style_function=lambda x:style1).add_to(border_layer)

map1.add_child(border_layer)
map1.add_child(cantons)
#map1.add_child(completed)
map1.add_child(districts)
map1.add_child(rides)
map1.add_child(stops)


map1.add_child(folium.LayerControl())  
map1.save( 'index.html')

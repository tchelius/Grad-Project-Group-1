import arcpy
import geopandas as gpd
from shapely.geometry import Point
from shapely.ops import unary_union
import numpy as np
import requests
import os
import matplotlib.pyplot as plt
from shapely.affinity import scale, rotate

# Set workspace
arcpy.env.workspace = r"E:\gis programming\gis_project\wildfire"
arcpy.env.overwriteOutput = True

# Define input paths
fuel_data = r"E:\gis_programming\gis_project\Vegetation_-_Western_San_Diego_County_[ds2964]\Vegetation_-_Western_San_Diego_County_[ds2964].shp"
slope_data = r"E:\gis_programming\gis_project\Slopes_CN_shapefile\Slopes_CN.shp"
evac_url = "https://gis-public.sandiegocounty.gov/arcgis/rest/services/Hosted/Zonehaven_Evacuation_Zones___SD_County/FeatureServer/105/query?where=1=1&outFields=*&f=geojson"
evac_zones = gpd.read_file(evac_url)

# Define simulation parameters
wind_speed = 45  # kph
duration = 5040  # minutes (7 days)
time_step = 480  # minutes (8 hours)
humidity = 6     # relative humidity %
start_lat, start_lon = 32.5947, -116.8437

# Ensure 'fuel_type' field exists
fields = [f.name for f in arcpy.ListFields(fuel_data)]
if 'fuel_type' not in fields:
    arcpy.AddField_management(fuel_data, 'fuel_type', 'TEXT', field_length=20)

with arcpy.da.UpdateCursor(fuel_data, ['GROUP_', 'fuel_type']) as cursor:
    for row in cursor:
        group_val = row[0]
        if group_val is None:
            row[1] = 'unknown'
        elif 'grass' in group_val.lower():
            row[1] = 'grassland'
        elif 'shrub' in group_val.lower() or 'scrub' in group_val.lower() or 'chaparral' in group_val.lower():
            row[1] = 'shrub'
        elif 'urban' in group_val.lower() or 'developed' in group_val.lower():
            row[1] = 'urban'
        else:
            row[1] = 'shrub'
        cursor.updateRow(row)

print("✅ 'fuel_type' field created and populated.")

# --- Helper Functions ---
def create_wind_ellipse(buffer_geom, wind_direction):
    stretch_x, stretch_y = 2.0, 1.0
    direction_angle = {"W": 270, "SW": 225, "E": 90, "NE": 45}.get(wind_direction, 0)
    ellipse = scale(buffer_geom, xfact=stretch_x, yfact=stretch_y)
    ellipse = rotate(ellipse, direction_angle, origin='centroid')
    return ellipse

def get_random_point_on_fire(fire_geometry):
    if fire_geometry.geom_type == "Point":
        return fire_geometry
    perimeter = fire_geometry.exterior
    random_distance = np.random.uniform(0, perimeter.length)
    return perimeter.interpolate(random_distance)

def get_slope_at_point(slope_shapefile, point):
    slope_gdf = gpd.read_file(slope_shapefile).to_crs("EPSG:3310")
    point_gdf = gpd.GeoDataFrame(geometry=[point], crs="EPSG:3310")
    joined = gpd.sjoin(point_gdf, slope_gdf, how="left", predicate="intersects")
    if not joined.empty and 'gridcode' in joined.columns:
        return joined['gridcode'].iloc[0]
    return 0

def get_fuel_type(fuel_data, point):
    gdf_fuel = gpd.read_file(fuel_data).to_crs("EPSG:3310")
    point_gdf = gpd.GeoDataFrame(geometry=[point], crs="EPSG:3310")
    joined = gpd.sjoin(point_gdf, gdf_fuel, how="left", predicate="intersects")
    if not joined.empty:
        if "fuel_type" in joined.columns:
            fuel = joined['fuel_type'].iloc[0]
            if fuel is None or str(fuel).lower() == 'nan':
                return "unknown"
            return fuel
    return "unknown"

def calculate_ros(fuel_type, wind_speed, slope, humidity):
    fuel_ros = {'grassland': 0.3, 'shrub': 0.7, 'forest': 1.5, 'urban': 0.1}
    wind_speed_m_per_min = (wind_speed * 1000) / 60
    wind_factor = 0.05 * wind_speed_m_per_min
    slope_factor = 0.03 * slope
    humidity_factor = 1 + (1 - humidity / 100)
    base_ros = fuel_ros.get(fuel_type, 0.1)
    ros = base_ros * (1 + wind_factor + slope_factor) * humidity_factor
    if ros is None or np.isnan(ros) or ros <= 0:
        ros = 0.1
    return ros

def get_wind_direction(t):
    if t < 1440:
        return "W"
    elif t < 2880:
        return "SW"
    elif t < 4320:
        return "E"
    else:
        return "NE"

# --- Simulate Fire Spread ---
def simulate_fire(start_point, fuel_data, slope_data, wind_speed, duration, time_step, humidity):
    current_fire = start_point
    fire_buffers = []
    
    for t in range(0, duration, time_step):
        if t == 0:
            fire_buffers.append(current_fire)
            print(f"Time {t} min: Fire ignition point established.")
            continue

        random_point = get_random_point_on_fire(current_fire)
        slope = get_slope_at_point(slope_data, random_point)
        fuel_type = get_fuel_type(fuel_data, random_point)
        ros = calculate_ros(fuel_type, wind_speed, slope, humidity)
        distance = ros * time_step
        distance = min(distance, 26.5)

        # Slow down if fuel is unknown
        if fuel_type == "unknown":
            distance *= 0.1  # 10% spread on unknown

        if not np.isfinite(distance) or distance <= 0:
            print(f"⚠️ Skipping timestep {t}: Invalid spread {distance}")
            continue

        if distance < 20:
            print(f"⚠️ Spread very small at timestep {t}: {distance:.2f}m — fire crawls slowly.")
            # Keep crawling, no break

        buffer_fire = current_fire.buffer(distance)
        wind_direction = get_wind_direction(t)
        buffer_fire = create_wind_ellipse(buffer_fire, wind_direction)

        if current_fire.is_empty:
            current_fire = buffer_fire
        else:
            current_fire = unary_union([current_fire, buffer_fire])

        fire_buffers.append(buffer_fire)
        print(f"Time {t} min: Wind {wind_direction}, Spread {distance:.2f}m, Fuel {fuel_type}, Slope {slope:.2f}%")

    if not fire_buffers:
        print("No fire spread detected. Check parameters.")
        exit()

    return fire_buffers

# --- Visualization ---
def visualize_fire(fire_buffers, evac_zones):
    evac_zones = evac_zones.to_crs("EPSG:3310")
    fig, ax = plt.subplots(figsize=(12, 12))
    step = max(1, len(fire_buffers) // 10)
    for i, fire in enumerate(fire_buffers):
        if i % step == 0:
            gpd.GeoSeries([fire]).plot(ax=ax, color='red', alpha=0.3)
    evac_zones.plot(ax=ax, color='blue', alpha=0.2)
    plt.legend(['Fire Perimeters', 'Evacuation Zones'])
    plt.title("Wildfire Spread Simulation (Selected Steps)")
    plt.show()

# --- Setup Start Point ---
start_point_geom = Point(start_lon, start_lat)
gdf_start = gpd.GeoDataFrame([{'geometry': start_point_geom}], crs="EPSG:4326").to_crs("EPSG:3310")
start_point = gdf_start.geometry.iloc[0]

# --- Run Simulation ---
fire_buffers = simulate_fire(start_point, fuel_data, slope_data, wind_speed, duration, time_step, humidity)

# --- Visualize Fire Spread ---
visualize_fire(fire_buffers, evac_zones)

# --- Final Acreage and Evac Zones Report ---
fire_final = fire_buffers[-1]
gdf_fire = gpd.GeoSeries([fire_final], crs="EPSG:3310")
burned_acres = gdf_fire.area.iloc[0] / 4046.86
print(f"🔥 Total burned area: {burned_acres:.2f} acres")

gdf_evacs = evac_zones.to_crs("EPSG:3310")
evacs_impacted = gdf_evacs[gdf_evacs.intersects(fire_final)]
print(f"🏡 Evacuation zones impacted: {len(evacs_impacted)} zones")


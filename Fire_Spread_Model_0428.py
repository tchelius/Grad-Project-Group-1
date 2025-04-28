import arcpy
from arcpy import CopyFeatures_management
import os

gbd_path = r"C:\Users\rebec\Documents\ArcGIS\Projects\MSState\GR6363 GIS Programming\Group Project\Group Project.gdb"
arcpy.env.workspace = gbd_path

class Toolbox(object):
    def __init__(self):
        self.label = "Wildfire Toolbox"
        self.alias = "wildfire"
        self.tools = [WildfireSimulation]

class WildfireSimulation(object):
    def __init__(self):
        self.label = "Run Wildfire Simulation"
        self.description = "Simulates wildfire spread given slope, fuel, and weather parameters."

    def getParameterInfo(self):
        params = [        
            arcpy.Parameter(
                displayName="Fuel Raster",
                name="fuel_data",
                datatype="DEFile",
                 parameterType="Required",
                direction="Input"
           ),
            arcpy.Parameter(
                displayName="Slope Shapefile",
                name="slope_data",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Input"
           ),
            arcpy.Parameter(
                displayName="Start Latitude",
                name="start_lat",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           ),
            arcpy.Parameter(
                displayName="Start Longitude",
                name="start_long",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           ),
            arcpy.Parameter(
                displayName="Wind Speed (kph)",
                name="wind_speed",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           ),
            arcpy.Parameter(
                displayName="Humidity (%)",
                name="humidity",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           ),
            arcpy.Parameter(
                displayName="Duration (minutes)",
                name="duration",
                datatype="Long",
                parameterType="Optional",
                direction="Input"
           ),
            arcpy.Parameter(
                displayName="Time Step (minutes)",
                name="time_step",
                datatype="Long",
                parameterType="Optional",
                direction="Input"
           ),
            arcpy.Parameter(
                displayName="Output Shapefile (Fire Perimeter)",
                name="output_fc",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Output"
           )
        ]

        #params[6].value = 10080
        #params[7].value = 720
        return params

    def execute(self, parameters, messages):
        #import packages and libraries
        import geopandas as gpd
        import numpy as np
        from shapely.geometry import Point
        from shapely.ops import unary_union
        from shapely.affinity import scale, rotate
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        import rasterio
        from shapely import wkt
        #import matplotlib_scalebar.scalebar as sb

        # Get parameters
        fuel_data = parameters[0].valueAsText
        slope_data = parameters[1].valueAsText
        start_lat = float(parameters[2].value)
        start_lon = float(parameters[3].value)
        wind_speed = float(parameters[4].value)
        humidity = float(parameters[5].value)
        duration = int(parameters[6].value)
        time_step = int(parameters[7].value)
        output_fc = parameters[8].valueAsText

        duration = 10080
        time_step = 720

        arcpy.AddMessage("📍 Initializing wildfire simulation...")

        # Convert to projected CRS
        start_point_geom = Point(start_lon, start_lat)  #uses shapely.geometry.Point to create a point geodatabase in WGS84 (= EPSG:4326)
        gdf_start = gpd.GeoDataFrame([{'geometry': start_point_geom}], crs="EPSG:4326").to_crs("EPSG:3310") #creates a geopandas dataframe from shapely point and converts to EPSG:3310
        start_point = gdf_start.geometry.iloc[0]        #defines 'start_point' by grapping point object from gdf_start using .iloc[0]

        slope_gdf = gpd.read_file(slope_data).to_crs("EPSG:3310")
        #fuel_gdf = gpd.read_file(fuel_data).to_crs("EPSG:3310")

        def create_wind_ellipse(buffer_geom, wind_direction):
            stretch_x, stretch_y = 2.0, 1.0
            direction_angle = {"W": 270, "SW": 225, "E": 90, "NE": 45}.get(wind_direction, 0)
            ellipse = scale(buffer_geom, xfact=stretch_x, yfact=stretch_y)
            return rotate(ellipse, direction_angle, origin='centroid')

        def get_random_point_on_fire(fire_geometry):
            if fire_geometry.geom_type == "Point":
                return fire_geometry
            perimeter = fire_geometry.exterior
            random_distance = np.random.uniform(0, perimeter.length)
            return perimeter.interpolate(random_distance)

        def get_slope_at_point(point):
            point_gdf = gpd.GeoDataFrame(geometry=[point], crs="EPSG:3310")
            joined = gpd.sjoin(point_gdf, slope_gdf, how="left", predicate="intersects")
            return joined['gridcode'].iloc[0] if not joined.empty and 'gridcode' in joined.columns else 0

        def get_fuel_type(point):
            point_gdf = gpd.GeoDataFrame(geometry=[point], crs="EPSG:3310")
            with rasterio.open(fuel_data) as src:           #need to make sure raster layer is in SR 3310
                sampled_values = [x for x in src.sample([(x.x, x.y) for x in point_gdf.geometry])]
                if not sampled_values or sampled_values[0] is None:
                    arcpy.AddWarning(f"⚠️ No valid fuel type found for point {point}. Returning 'unknown'.")
                    return "unknown"
                fuel_type = sampled_values[0][0]
                if fuel_type == 32767:
                    return 0.1
            print(fuel_type)
            return fuel_type

            #point_gdf = gpd.GeoDataFrame(geometry=[point], crs="EPSG:3310")
            #joined = gpd.sjoin(point_gdf, fuel_gdf, how="left", predicate="intersects")
            #if not joined.empty and "fuel_type" in joined.columns:
            #   fuel = joined['fuel_type'].iloc[0]
            #    return "unknown" if fuel is None or str(fuel).lower() == 'nan' else fuel
            #return "unknown"

        def calculate_ros(fuel_type, wind_speed, slope, humidity):
            fuel_ros = {'2': 0.3, '4': 0.7, '5': 1.5, '7': 0.1, '8': 0.1, '9': 0.1, '10': 0.1, '11': 0.1, '12': 0.1, '91': 0.1, '93': 0.1, '98': 0.1, '99': 0.1}
            wind_factor = 0.05 * (wind_speed * 1000 / 60)
            slope_factor = 0.03 * slope
            humidity_factor = 1 + (1 - humidity / 100)
            base_ros = fuel_ros.get(fuel_type, 0.1)
            ros = base_ros * (1 + wind_factor + slope_factor) * humidity_factor
            return max(ros, 0.1) if np.isfinite(ros) else 0.1

        def get_wind_direction(t):
            return "W" if t < 2880 else "SW" if t < 5760 else "E" if t < 7920 else "NE"

        # Run simulation
        current_fire = start_point
        fire_buffers = []

        for t in range(0, duration, time_step):
            if t == 0:
                fire_buffers.append(current_fire)
                arcpy.AddMessage(f"⏱️ Time {t} min: Fire ignition established.")
                continue

            pt = get_random_point_on_fire(current_fire)
            slope = get_slope_at_point(pt)
            fuel_type = get_fuel_type(pt)
            ros = calculate_ros(fuel_type, wind_speed, slope, humidity)
            distance = min(ros * time_step, 9.25)
            if fuel_type == "unknown": distance *= 0.1
            if not np.isfinite(distance) or distance <= 0:
                arcpy.AddWarning(f"⚠️ Time {t}: Invalid spread distance {distance}.")
                continue

            buffer_fire = current_fire.buffer(distance)
            buffer_fire = create_wind_ellipse(buffer_fire, get_wind_direction(t))
            current_fire = unary_union([current_fire, buffer_fire]) if not current_fire.is_empty else buffer_fire
            fire_buffers.append(buffer_fire)
            arcpy.AddMessage(f"⏱️ Time {t} min: Spread {distance:.2f}m, Fuel {fuel_type}, Slope {slope}%")

        if not fire_buffers:
            raise RuntimeError("🔥 No fire spread occurred. Check input data or parameters.")

        # Export final buffer
        output_fc = "wildfire_result"
        output_fc_path = os.path.join(gbd_path, output_fc)
        spatial_ref = arcpy.SpatialReference(3310)
        arcpy.management.CreateFeatureclass(gbd_path, output_fc, "POLYGON", spatial_reference=spatial_ref)
        last_fire_geo = fire_buffers[-1]
        wkt_str = last_fire_geo.wkt #convert shapely geometry of last_fire_geo to WKT
        arc_geo = arcpy.FromWKT(wkt_str, spatial_ref)
        fields = ['SHAPE@']
        with arcpy.da.InsertCursor(output_fc_path, fields) as cursor:
            cursor.insertRow([arc_geo])

        #gdf_result = gpd.GeoDataFrame(geometry=[fire_buffers[-1]], crs="EPSG:3310")
        #gdf_result.to_file("temp_output.shp")
        #arcpy.CopyFeatures_management("temp_output.shp", output_fc)
        #burned_acres = gdf_result.area.iloc[0] / 4046.86
        #arcpy.AddMessage(f"🔥 Total burned area: {burned_acres:.2f} acres")

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return
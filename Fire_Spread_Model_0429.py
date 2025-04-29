import arcpy
from arcpy import CopyFeatures_management
import os

gbd_path = r"C:\Users\rebec\Documents\ArcGIS\Projects\MSState\GR6363 GIS Programming\Group Project\Group Project.gdb"       #your gbd path here
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
        params = []        
        fuel_data = arcpy.Parameter(
                displayName="Fuel Raster",
                name="fuel_data",
                datatype="DEFile",
                 parameterType="Required",
                direction="Input"
           )
        slope_data = arcpy.Parameter(
                displayName="Slope Shapefile",
                name="slope_data",
                datatype="DEFeatureClass",
                parameterType="Required",
                direction="Input"
           )
        start_lat = arcpy.Parameter(
                displayName="Start Latitude",
                name="start_lat",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           )
        start_long = arcpy.Parameter(
                displayName="Start Longitude",
                name="start_long",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           )
        wind_speed = arcpy.Parameter(
                displayName="Wind Speed (kph)",
                name="wind_speed",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           )
        wind_dir = arcpy.Parameter(
                displayName="Wind Direction (degrees)",
                name="wind_dir",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           )
        rh = arcpy.Parameter(
                displayName="Relative Humidity (%)",
                name="humidity",
                datatype="Double",
                parameterType="Required",
                direction="Input"
           )
        duration = arcpy.Parameter(
                displayName="Duration (minutes)",
                name="duration",
                datatype="Long",
                parameterType="Optional",
                direction="Input"
           )
        time_step = arcpy.Parameter(
                displayName="Time Step (minutes)",
                name="time_step",
                datatype="Long",
                parameterType="Optional",
                direction="Input"
           )
        output_fc = arcpy.Parameter(
                displayName="Output Shapefile (Fire Perimeter)",
                name="output_fc",
                datatype="DEFeatureClass",
                parameterType="Derived",
                direction="Output"
           )   
        
        wind_speed.filter.type = "Range"
        wind_speed.filter.list = [0,100]        #limits wind speed to nonnegative values between 0 and 100 kph
        wind_dir.filter.type = "Range"
        wind_dir.filter.list = [0,360]          #limits wind direction to 360 degrees of a compass rose
        rh.filter.type = "Range" 
        rh.filter.list = [0, 100]       #limit RH values to between 0 and 100
        duration.filter.type = "Range"
        duration.filter.list = [0, 360]     #limit predicted fire duration to 360 mins = 6 hrs

        params = [fuel_data, slope_data, start_lat, start_long, wind_speed, wind_dir, rh, duration, time_step, output_fc]
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
        from shapely.wkt import loads as load_wkt
        import math
        #import matplotlib_scalebar.scalebar as sb

        # Get parameters
        fuel_data = parameters[0].valueAsText
        slope_data = parameters[1].valueAsText
        start_lat = float(parameters[2].value)
        start_lon = float(parameters[3].value)
        wind_speed = float(parameters[4].value)
        wind_dir = float(parameters[5].value)
        humidity = float(parameters[6].value)
        duration = int(parameters[7].value)
        time_step = int(parameters[8].value)
        output_fc = parameters[9].valueAsText

        #set up arcpy crs
        sr = arcpy.SpatialReference(3310)

        arcpy.AddMessage("📍 Initializing wildfire simulation...")

        #import start point using gpd AND arcpy! Shapely version gets fed into model; arcpy version exports as feature layer to final output.
        #shapely channel
        start_point_geom = Point(start_lon, start_lat)  #uses shapely.geometry.Point to create a point geodatabase in WGS84 (= EPSG:4326)
        gdf_start = gpd.GeoDataFrame([{'geometry': start_point_geom}], crs="EPSG:4326").to_crs("EPSG:3310") #creates a geopandas dataframe from shapely point and converts to EPSG:3310
        start_point = gdf_start.geometry.iloc[0]        #defines 'start_point' by grapping point object from gdf_start using .iloc[0]

        #arcpy channel
        pt = arcpy.Point(start_lon, start_lat)          #Create start arcpy point geometry in geographic coords
        pt_sr = arcpy.SpatialReference(4326)
        pt_geo = arcpy.PointGeometry(pt, spatial_reference=pt_sr)
        start_fc_geo="Fire_Origin_Geometric"  #create empty Feature Class which will hold fire origin point
        arcpy.management.CreateFeatureclass(gbd_path, start_fc_geo, "POINT", spatial_reference=pt_sr)
        start_fc_geo = os.path.join(gbd_path, start_fc_geo)         #add point geometry to feature class using search cursor
        fields = ['SHAPE@']     #create fields to be added to feature class
        with arcpy.da.InsertCursor(start_fc_geo, fields) as cursor:     #add start point using InsertCursor
            cursor.insertRow([pt_geo])
        #project into 3310 - not sure if this step is totally necessary
        start_proj_name = "Fire_Origin_Point"
        start_fc = os.path.join(gbd_path, start_proj_name)
        input_feature_class = arcpy.management.Project(start_fc_geo, start_fc, sr)

        #read in slope data via gpd
        slope_gdf = gpd.read_file(slope_data).to_crs("EPSG:3310")

        #define helper functions

        def get_slope_at_point(point):          #intersect slope layer and start point using gpd and extract slope value
            point_gdf = gpd.GeoDataFrame(geometry=[point], crs="EPSG:3310")
            joined = gpd.sjoin(point_gdf, slope_gdf, how="left", predicate="intersects")
            slope_data_value = joined['gridcode'].iloc[0] if not joined.empty and 'gridcode' in joined.columns else 0
            arcpy.AddMessage(f'Slope value from slope shapefile: {slope_data_value}')
            return slope_data_value

        def get_fuel_type(point):
            point_gdf = gpd.GeoDataFrame(geometry=[point], crs="EPSG:3310")
            with rasterio.open(fuel_data) as src:           #need to make sure raster layer is in SR 3310
                if point_gdf.crs != src.crs:                #check if point crs is same as raster crs
                    point_gdf = point_gdf.to_crs(src.crs)   #and convert if not the same - better than reprojecting entire raster
                sampled_values = [x for x in src.sample([(x.x, x.y) for x in point_gdf.geometry])]
                if not sampled_values or sampled_values[0] is None:
                    arcpy.AddWarning(f"⚠️ No valid fuel type found for point {point}. Returning 'unknown'.")
                    return "unknown"
                fuel_type = sampled_values[0][0]
                if fuel_type == 32767:
                    arcpy.AddWarning(f"⚠️ No valid fuel type found for point {point}. Returning default 0.1")
                    return 0.1
            arcpy.AddMessage(f"Fuel value from fuel raster: {fuel_type}")
            return fuel_type

        def calculate_ros(fuel_type, wind_speed, slope, humidity):              #lay out mathematical relationships between all parameters
            fuel_ros = {'2': 0.05, '4': 0.2, '5': 0.3, '7': 0.4, '8': 0.5, '9': 0.6, '10': 0.7, '11': 0.8, '12': 0.9, '91': 1.0, '93': 1.1, '98': 1.2, '99': 1.3}
            base_ros = fuel_ros.get(fuel_type, 0.1)
            wind_factor = 0.05 * (wind_speed * 1000 / 60)
            slope_factor = 0.03 * slope
            humidity_factor = 1 + (1 - humidity / 100)
            ros = base_ros * (1 + wind_factor + slope_factor) * humidity_factor
            arcpy.AddMessage("Rate of spread calculated.")
            return max(ros, 0.1) if np.isfinite(ros) else 0.1
        
        def create_fire_ellipse(wind_speed, wind_dir, out_ellipse_path):              #use arcpy tools (incl. arcpy start point) to create fire spread ellipse
            major = ros
            minor = (ros + 0.1)/((wind_speed * 0.15) +1)
            azimuth = wind_dir

            #calculate projected coordinates for TableToEllipse function
            arcpy.management.CalculateGeometryAttributes(input_feature_class, [['Long_Proj', "POINT_X"],
                                                                            ['Lat_Proj', 'POINT_Y']])
            
            wind_dir_rad = math.radians((wind_dir + 180) % 360)

            #this next section provides logic which allows the tool to add a new field if it is not already created
            #and as doing so, specifies a NUMERIC field type -- this is necessary for the TableToEllipse function!
            fields_to_add = [['Major', "Float"], 
                            ['Minor', "Float"], 
                            ['Azimuth', "Float"],
                            ['x_shifted', "Float"],
                            ['y_shifted', "Float"]]
            
            existing_fields = [f.name for f in arcpy.ListFields(input_feature_class)]

            for field in fields_to_add:
                if field not in existing_fields:
                    arcpy.management.AddField(input_feature_class, field[0], field[1])
            
            #use UpdateCursor to apply mathematical defns of major and minor axes and azimuth direction as specified above
            with arcpy.da.UpdateCursor(input_feature_class, ['Long_Proj', 'Lat_Proj', 'Major', 'Minor', 'Azimuth', 'x_shifted', 'y_shifted']) as cursor:
                for row in cursor:
                    long_proj, lat_proj = row[0], row[1]

                    scale_factor = (wind_speed/100)*1000
                    shift_dist = scale_factor*(major/2)
                    x_shift = math.sin(wind_dir_rad)*shift_dist
                    y_shift = math.cos(wind_dir_rad)*shift_dist

                    row[2] = major
                    row[3] = minor
                    row[4] = azimuth
                    row[5] = long_proj + x_shift
                    row[6] = lat_proj + y_shift

                    cursor.updateRow(row)

            arcpy.management.TableToEllipse(input_feature_class, out_ellipse_path, 'x_shifted', 'y_shifted', 'Major', 'Minor', "KILOMETERS", 'Azimuth', "DEGREES", "", "", "", "POLYGON","PLANAR")
            
            #use SearchCursor to use WKT to return shapely (not arcpy) ellipse geometry
            with arcpy.da.SearchCursor(out_ellipse_path, ['SHAPE@WKT']) as cursor:
                for row in cursor:
                    return load_wkt(row[0])
            return

        def get_random_point_on_fire(fire_geometry):
            if fire_geometry.geom_type == "Point":
                return fire_geometry
            perimeter = fire_geometry.exterior
            random_distance = np.random.uniform(0, perimeter.length)
            return perimeter.interpolate(random_distance)

        # Run simulation
        current_fire = start_point      #set initially 'current_fire' as shapely start_point
        fire_buffers = []               #create empty list for buffers (created in subsequent steps)

        for t in range(0, duration, time_step):         #sets up loop to calculate fire buffer ellipses at different timesteps
            if t == 0:                                  #if the fire is just starting, adds initial location to list of buffer geometries
                fire_buffers.append(current_fire)
                arcpy.AddMessage(f"⏱️ Time {t} min: Fire ignition established.")
                continue

            pt = get_random_point_on_fire(current_fire)     #grabs random point on perimeter of fire ellipse
            slope = get_slope_at_point(pt)                  #grabs slope at that point
            fuel_type = get_fuel_type(pt)                   #fuel, too
            ros = calculate_ros(fuel_type, wind_speed, slope, humidity)         #calculates ros for that point
            distance = min(ros * time_step, 9.25)           #calculates buffer distance
            if fuel_type == "unknown": distance *= 0.1
            if not np.isfinite(distance) or distance <= 0:
                arcpy.AddWarning(f"⚠️ Time {t}: Invalid spread distance {distance}.")
                continue

            buffer_fire = current_fire.buffer(distance)
            out_ellipse_path = os.path.join(gbd_path, f'Fire_spread_ellipse_timestep_{t}')
            buffer_fire = create_fire_ellipse(wind_speed, wind_dir, out_ellipse_path)
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
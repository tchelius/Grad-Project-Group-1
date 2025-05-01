import arcpy
from arcpy import CopyFeatures_management
import os

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
        start_lon = arcpy.Parameter(
                displayName="Start Longitude",
                name="start_lon",
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
                displayName="Workspace (GDB or Folder)",
                name="gbd_path",
                datatype="DEWorkspace",
                parameterType="Required",
                direction="Input"
           )   
        evac_zones = arcpy.Parameter(
                displayName="Evacuation Zones",
                 name="evac_zones",
                datatype="DEFile",
                parameterType="Optional",
                direction="Input"
           )   
        
        
        wind_speed.filter.type = "Range"
        wind_speed.filter.list = [0,100]        #limits wind speed to nonnegative values between 0 and 100 kph
        wind_dir.filter.type = "Range"
        wind_dir.filter.list = [0,360]          #limits wind direction to 360 degrees of a compass rose
        rh.filter.type = "Range" 
        rh.filter.list = [0, 100]             #limit RH values to between 0 and 100
        duration.filter.type = "Range"
        duration.filter.list = [0, 360]         #limit predicted fire duration to 360 mins = 6 hrs

        params = [fuel_data, slope_data, start_lat, start_lon, wind_speed, wind_dir, rh, duration, time_step, output_fc, evac_zones]
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
        gbd_path = parameters[9].valueAsText
        evac_zones = parameters[10].valueAsText

        #set up arcpy crs, workspace
        arcpy.env.workspace = gbd_path
        sr = arcpy.SpatialReference(3310)

        arcpy.AddMessage("📍 Initializing wildfire simulation...")

        #unpack start point using shapely AND arcpy! Shapely version gets fed into model; arcpy version exports as feature layer to final output

        #shapely channel
        start_point_geom = Point(start_lon, start_lat)  
        gdf_start = gpd.GeoDataFrame([{'geometry': start_point_geom}], crs="EPSG:4326").to_crs("EPSG:3310") 
        start_point = gdf_start.geometry.iloc[0]

        #arcpy channel
        #unpack geographic coords of start point and add to feature class
        pt = arcpy.Point(start_lon, start_lat)
        pt_sr = arcpy.SpatialReference(4326)
        pt_geo = arcpy.PointGeometry(pt, spatial_reference=pt_sr)
        start_fc_geo="Fire_Origin_Geometric"          #create empty Feature Class which will hold geometric fire origin point
        arcpy.management.CreateFeatureclass(gbd_path, start_fc_geo, "POINT", spatial_reference=pt_sr)
        start_fc_geo = os.path.join(gbd_path, start_fc_geo)         #add point geometry to feature class using search cursor
        fields = ['SHAPE@']     
        with arcpy.da.InsertCursor(start_fc_geo, fields) as cursor:
            cursor.insertRow([pt_geo])
        #project geographic feature class into 3310
        start_proj_name = "Fire_Origin_Point"
        start_fc = os.path.join(gbd_path, start_proj_name)
        arcpy.management.Project(start_fc_geo, start_fc, sr)

        #add to current map
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        map_obj = aprx.activeMap
        map_obj.addDataFromPath(start_fc)

        #read in shapely slope data with gpd
        slope_gdf = gpd.read_file(slope_data).to_crs("EPSG:3310")

        #define helper functions

        def get_slope_at_point(point):          #intersect slope layer and shapely start point using gpd and extract slope value
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
            fuel_ros = {'2': 2.5, '4': 6, '5': 0.5, '7': 5, '8': 1.0, '9': 7, '10': 8, '11': 4, '12': 0.9, '91': 9, '93': 5, '98': 0, '99': 0}
            base_ros = fuel_ros.get(fuel_type, 0.1)
            wind_factor = 0.05 * (wind_speed * 1000 / 60)
            slope_factor = 0.03 * slope
            humidity_factor = 1 + (1 - humidity / 100)
            ros = base_ros * (1 + wind_factor + slope_factor) * humidity_factor
            arcpy.AddMessage(f"Rate of spread calculated: {ros}")
            return max(ros, 0.1) if np.isfinite(ros) else 0.1
        
        def create_fire_ellipse(start_fc, wind_speed, wind_dir, ros, out_ellipse_path):
            major = ros
            minor = (ros + 0.1)/((wind_speed * 0.15) +1)
            azimuth = wind_dir
    
            #this next section provides logic which allows the tool to add a new field if it is not already created
            #and as doing so, specifies a NUMERIC field type -- this is necessary for the TableToEllipse function!
            fields_to_add = [['Major', "Float"], ['Minor', "Float"], ['Azimuth', "Float"], ['x_shifted', "Float"], ['y_shifted', "Float"]]
            existing_fields = [f.name for f in arcpy.ListFields(start_fc)]
            for field in fields_to_add:
                if field not in existing_fields:
                    arcpy.management.AddField(start_fc, field[0], field[1])

            #calculate projected coordinates for TableToEllipse function
            arcpy.management.CalculateGeometryAttributes(start_fc, [['Long_Proj', "POINT_X"], ['Lat_Proj', 'POINT_Y']])
            wind_dir_rad = math.radians((wind_dir + 180) % 360)

            #use UpdateCursor to apply mathematical defns of major and minor axes and azimuth direction as specified above
            with arcpy.da.UpdateCursor(start_fc, ['Long_Proj', 'Lat_Proj', 'Major', 'Minor', 'Azimuth', 'x_shifted', 'y_shifted']) as cursor:
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
            arcpy.management.TableToEllipse(start_fc, out_ellipse_path, 'x_shifted', 'y_shifted', 'Major', 'Minor', "KILOMETERS", 'Azimuth', "DEGREES", "", "", "", "POLYGON","PLANAR")
            
            #use SearchCursor to use WKT to return shapely (not arcpy) ellipse geometry
            with arcpy.da.SearchCursor(out_ellipse_path, ['SHAPE@WKT']) as cursor:
                for row in cursor:
                    return load_wkt(row[0])
            return

        def get_random_point_on_fire(fire_geometry):
            if fire_geometry.geom_type == "Point":
                return fire_geometry
            elif fire_geometry.geom_type == "Polygon":
                perimeter = fire_geometry.exterior
            elif fire_geometry.geom_type == "MultiPolygon":
                largest = max(fire_geometry.geoms, key=lambda g: g.area)
                perimeter = largest.exterior
            else:
                raise TypeError(f"Unsupported geometry type for fire perimeter: {fire_geometry.geom_type}")
            random_distance = np.random.uniform(0, perimeter.length)
            return perimeter.interpolate(random_distance)

        # Run simulation
        fire_buffers = []               #create empty list for buffers (created in subsequent steps)

        for t in range(0, duration, time_step):         #sets up loop to calculate fire buffer ellipses at different timesteps
            if t == 0:                                  #if the fire is just starting, adds initial location to list of buffer geometries
                arcpy.AddMessage(f"⏱️ Time {t} min: Fire ignition established.")
                pt = start_point
                slope = get_slope_at_point(pt)                  #grabs slope at that point
                fuel_type = get_fuel_type(pt)                   #fuel, too
                ros = calculate_ros(fuel_type, wind_speed, slope, humidity)
                out_ellipse_path = os.path.join(gbd_path, f'Fire_spread_ellipse_timestep_{t}')
                fire_ellipse = create_fire_ellipse(start_fc, wind_speed, wind_dir, ros, out_ellipse_path)
                fire_buffers.append(fire_ellipse)
                arcpy.AddMessage(f"🔥 Time {t+1} min: Initial fire ellipse created.")
                continue
            if t >= 1:
                prev_fire = fire_buffers[-1]
                pt = get_random_point_on_fire(prev_fire)
                slope = get_slope_at_point(pt)                  
                fuel_type = get_fuel_type(pt)                   
                ros = calculate_ros(fuel_type, wind_speed, slope, humidity)         
                distance = max(ros * time_step, 9.25)           
                if fuel_type == "unknown": 
                    distance *= 0.1
                if not np.isfinite(distance) or distance <= 0:
                    arcpy.AddWarning(f"⚠️ Time {t}: Invalid spread distance {distance}.")
                current_fire = prev_fire.buffer(distance)
                fire_buffers.append(current_fire)
                arcpy.AddMessage(f"⏱️ Time {t} min: Spread {distance:.2f}m, Fuel {fuel_type}, Slope {slope}%")
                continue

        if not fire_buffers:
            raise RuntimeError("🔥 No fire spread occurred. Check input data or parameters.")

        #Export buffers and ellipse
        output_fc_path = os.path.join(gbd_path, "Fire_Buffer_Ellipses")
        if arcpy.Exists(output_fc_path):
            arcpy.management.Delete(output_fc_path)
        arcpy.management.CreateFeatureclass(gbd_path, "Fire_Buffer_Ellipses", "POLYGON", spatial_reference=sr)
        fields = ['SHAPE@']
        with arcpy.da.InsertCursor(output_fc_path, fields) as cursor:
            for idx, geom in enumerate(fire_buffers):
                try:
                    wkt_str = geom.wkt
                    arc_geo = arcpy.FromWKT(wkt_str, sr)
                    cursor.insertRow([arc_geo])
                except Exception as e:
                    arcpy.AddWarning(f"Skipped fire buffer geometry at index {idx}:{e}.")
        map_obj.addDataFromPath(output_fc_path)

        #intersect with evacuation zones
        if evac_zones and os.path.exists(evac_zones):
            arcpy.AddMessage("🧭 Calculating intersection with evacuation zones...")
            #Save last fire buffer as a temporary feature class
            last_buffer_geom = fire_buffers[-1]
            last_buffer_wkt = last_buffer_geom.wkt
            last_buffer_arc = arcpy.FromWKT(last_buffer_wkt, sr)
            temp_last_buffer_fc = os.path.join(gbd_path, "Temp_Last_Fire_Buffer")
            if arcpy.Exists(temp_last_buffer_fc):
                arcpy.management.Delete(temp_last_buffer_fc)
            arcpy.management.CreateFeatureclass(gbd_path, "Temp_Last_Fire_Buffer", "POLYGON", spatial_reference=sr)
            with arcpy.da.InsertCursor(temp_last_buffer_fc, ['SHAPE@']) as cursor:
                cursor.insertRow([last_buffer_arc])

            #Set up for select by location - set up output path, make intersection layer (necessary since imported as DEFile), perform selection.
            affected_areas_fc = os.path.join(gbd_path, "Affected_Areas")
            evac_layer = "Evac_Zones_Layer"
            arcpy.management.MakeFeatureLayer(evac_zones, evac_layer)
            arcpy.management.SelectLayerByLocation(in_layer=evac_layer, overlap_type="INTERSECT", select_features=temp_last_buffer_fc, selection_type="NEW_SELECTION")

            #Create and add to map affected areas
            affected_areas_fc = os.path.join(gbd_path, "Affected_Areas")
            arcpy.management.CopyFeatures(evac_layer, affected_areas_fc)
            map_obj.addDataFromPath(affected_areas_fc)
            arcpy.AddMessage("🚨 Affected areas calculated and added to map.")
        else:
            arcpy.AddMessage("ℹ️ No evacuation zones provided or file not found — skipping affected areas step.")

    def postExecute(self, parameters):
        """This method takes place after outputs are processed and
        added to the display."""
        return
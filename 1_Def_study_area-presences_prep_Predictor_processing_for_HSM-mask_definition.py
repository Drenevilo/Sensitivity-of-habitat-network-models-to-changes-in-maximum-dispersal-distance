#-------------------------------------------------------------------------------
# Script name:       Def_study_area_and_presences-Predictor_processing_for_HSM-mask_definition
# From article:      Predicting species occurrences with habitat network models

# Purpose:      Definition of Study area and species records,
#               processing of predictors for Habitat Suitability Modelling,
#               and mask definition

# Author:      Damian O. Ortiz-Rodriguez, Antoine Guisan, Rolf Holderegger, Maarten J. van Strien
# Created:     10/10/2017

#-------------------------------------------------------------------------------

def main():
    pass

if __name__ == '__main__':
    main()


# Import system modules
import arcpy
from arcpy import env
from arcpy.sa import *

#For license usage
arcpy.CheckExtension('Spatial')
arcpy.CheckOutExtension('Spatial')

# ### Study area delimitation ########################
#Get shape of all Mittelland
#Select all regions labeled as 'Mittelland' as defined in "Biogeographic regions of Switzerland" by OFEV (2011))
biogreg = "M:\\people\\damiano\\PhD\\GeoData\\All_layers\\CH_Biogeographic_Regions\\biogreg.shp"
Mittelland = "M:\\people\\damiano\\PhD\\GeoData\\All_layers\\CH_Biogeographic_Regions\\Mittelland.shp"
arcpy.Select_analysis(biogreg, Mittelland, "BIOGREG_C6 IN (2)")

#Get Switzerland (exclude Liechtenstein)
Landesgebiet = "M:\\people\\damiano\\PhD\\GeoData\\All_layers\\swissBOUNDARIES3D\\BOUNDARIES_2016\\DATEN\\swissBOUNDARIES3D\\SHAPEFILE_LV03_LN02\\swissBOUNDARIES3D_1_3_TLM_LANDESGEBIET.shp
Switzerland = "M:\\people\\damiano\\PhD\\GeoData\\All_layers\\CH_Biogeographic_Regions\\Switzerland.shp"
arcpy.Select_analysis(Landesgebiet, Switzerland, "ICC IN (CH)")

#Buffer = 2Km: commonly reported amphibian dispersal distances; Smith & Green 2005, to prevent border effect
CH_border_buffer_2Km = "M:\\people\\damiano\\PhD\\GeoData\\All_layers\\CH_Biogeographic_Regions\\CH_border_buffer_2Km.shp"
arcpy.buffer(Switzerland, CH_border_buffer_2Km, "-2000 Meters", "FULL", "ROUND", "ALL", "", "PLANAR")

#clip Switzerland to its buffer inner (CH_border_buffer_2Km.shp)
#Use 'Erase' on Switzerland.shp, erasing the CH_border_buffer_2Km.shp
#Output:
CH_minus2Kmbuffer ="M:\\people\\damiano\\PhD\\GeoData\\All_layers\\CH_Biogeographic_Regions\\CH_minus2Kmbuffer.shp"
#process
arcpy.Erase_analysis(Switzerland, CH_border_buffer_2Km, CH_minus2Kmbuffer)

#Clip 'Mittelland.shp' with the -2 Km buffer to Swiss border
Core_area_Hy_a_multipart = "M:\\people\\damiano\\PhD\\Automated_Run_BioCHECNET\\GeoData\\New Study Area\\Core_area_Hy_a_multipart.shp"
arcpy.Clip_analysis(Mittelland, CH_minus2Kmbuffer, Core_area_Hy_a_multipart, "")

#Core_area_Hy_a_multipart has a multipart (non-contiguous) feature
#Make this feature its own feature class
##Use Multipart To Singlepart to get rid of the islands/noncontinuous parts.
Hochrhein_Biogreg_multi = "M:\\people\\damiano\\PhD\\Automated_Run_BioCHECNET\\GeoData\\New Study Area\\Hochrhein_Biogreg_multi.shp"
arcpy.MultipartToSinglepart_management(Core_area_Hy_a_multipart, Hochrhein_Biogreg_multi)

#Select - to get only the the contiguous parts of Hochrhein as independent .shp
Hochrhein_Biogreg_main = "M:\\people\\damiano\\PhD\\Automated_Run_BioCHECNET\\GeoData\\New Study Area\\Hochrhein_Biogreg_main.shp"
arcpy.Select_analysis(Hochrhein_Biogreg_multi, Hochrhein_Biogreg_main, "BIOGREG_C6 IN (2)")

#Select by attribute from Core_area_Hy_a_multipart to get whole study area except Hochrhein
Core_area_Hy_a_noHochrhein = "M:\\people\\damiano\\PhD\\Automated_Run_BioCHECNET\\GeoData\\New Study Area\\Core_area_Hy_a_noHochrhein.shp"
arcpy.Select_analysis(Core_area_Hy_a_multipart, Core_area_Hy_a_noHochrhein, "FID in (1,2,3)")
##FID in (1,2,3) - All the Mittelland except Hochrhein

#join the main part (biggest, contiguous) to the original feature class
#Merge
Core_area_Hy_a_class = "M:\\people\\damiano\\PhD\\Automated_Run_BioCHECNET\\GeoData\\New Study Area\\Core_area_Hy_a_class.shp"
arcpy.Merge_management([Core_area_Hy_a_noHochrhein, Hochrhein_Biogreg_main], Core_area_Hy_a_class)

#Dissolve to get 1 polygon of all
##Field on which to aggregate: ORIG_FID
##Output: Core_area_Hy_a_dissolved.shp
Study_area = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\New_Study_Area\\Core_area_Hy_a_dissolved.shp"
##This is the used study area to crop the records and all the variables
arcpy.Dissolve_management(Core_area_Hy_a_class, Study_area, "ORIG_FID")




# ### Species records spatial preparation ########################
#set workspace
#inDir1 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\"
#env.workspace = inDir1

#Generate 'species raster' for snapping all the layers to it
#Out of .shp with all records of group of interest (amphibians) in the study area
#?Presence_points_all? -> Location data of all amphibian observations, provided by InfoSpecies-KARCH. Converted into point shapefile
#Extent = Whole study area without (border) buffer
    #Point to Raster
        # local variables
Presence_points_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Centered_amph_objects_Mittelland_20160414.shp"
Presence_points_all_raster = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Centered_amph_objects_Mittelland.tif"
        # Execute PointToRaster
arcpy.PointToRaster_conversion(Presence_points_all, "ART", Presence_points_all_raster, "MOST_FREQUENT", "NONE", "100")
#meaning of terms: (inFeature, valField, outRaster, assignmentType, priorityField, cellSize)

    #Raster-Clip all records outside of the (larger) study area: All regions labeled as 'Mittelland' as defined by OFEV (2011)
        #  More local variables
#These are environments for most processes
Mittelland_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\New_Study_Area\\Mittelland.shp" #Extent
spp_records_raster = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\amph_objs_Mitt_all.tif" #Snap raster

arcpy.Clip_management (Presence_points_all_raster, "485410 109645 768768 295145", spp_records_raster, Mittelland_all, "", "NONE", "NO_MAINTAIN_EXTENT")
#meaning of terms: (in_raster, rectangle, out_raster, {in_template_dataset}, {nodata_value}, {clipping_geometry}, {maintain_clipping_extent})

#global var, final clip for all processes
Study_area = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\New_Study_Area\\Core_area_Hy_a_dissolved.shp"
Presence_points_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Presence_points_core.shp"
# ####### Do separate shp's for every species #########
#Crop spp shp's to study area
arcpy.Clip_analysis(Presence_points_all, Study_area, Presence_points_core, "")

#Take out salamanders for Mask analysis, as they are not pond-based
outSHP = "E:/PhD/Automated_Run_BioCHECNET/GeoData/Mask/Presence_points_core_no_sal.shp"
arcpy.Select_analysis (Presence_points_core, outSHP, "\"ART\" NOT LIKE 'Salamandra%'")

#Separate the spp records into individual .shp's per spp per year
#Adapted from comment on:http://gis.stackexchange.com/questions/9998/exporting-feature-class-into-multiple-feature-classes-based-on-field-values-usin
#http://gis.stackexchange.com/a/44435
# "my_shapefile" = Presence_points_core
import arcgisscripting
# Starts Geoprocessing
gp = arcgisscripting.create(10.3)
#changed version of ArcGIS
gp.OverWriteOutput = 1
#Set Input Output variables
#inputFile = Presence_points_core
outDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Ind_spp\\" #<-- CHANGE
#Don't forget to also put the double backslashes at the end of every directory,
#otherwise it will just be added as part of the name, and will be stored in the directory which was last followed by the double slashes

# Reads My_shapefile for different values in the attribute
rows = gp.searchcursor(Presence_points_core)
row = rows.next()
attribute_types = set([])

while row:
    attribute_types.add(row.ART) #<-- CHANGE to the name of specific attribute ("ART" = 'kind' (species) in German
    row = rows.next()

# Output a Shapefile for each different attribute
for each_attribute in attribute_types:
    outName = each_attribute.replace(" ", "_")
    outName = outName.replace(".", "")
    outSHP = outDir + outName + u".shp"
    print outSHP
    arcpy.Select_analysis (Presence_points_core, outSHP, "\"ART\" = '" + each_attribute + "'")
    arcpy.TableToTable_conversion(outSHP, outDir, outName+".csv")

del rows, row, attribute_types, gp

#END

# ###Do Shp's for ind. spp. by year #######################
#Set Output variables
inDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Ind_spp\\"
outDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Ind_spp\\By_year\\" #<-- CHANGE

#Import functions
import arcgisscripting
import os
import enum
import dbf

#Define functions
def unique_values(table, field): #From: https://arcpy.wordpress.com/2012/02/01/create-a-list-of-unique-field-values/
    with arcpy.da.SearchCursor(table, [field]) as cursor:
        return sorted({row[0] for row in cursor})

# Starts Geoprocessing
gp = arcgisscripting.create(10.4)
gp.OverWriteOutput = 1

#Set the working directory
os.chdir(inDir)

# do the same for all the shapefiles in the directory
#Set folder from which files are selected and processed
for filename in os.listdir(inDir):
    if filename.endswith(".shp"):
        years = unique_values(filename, "C_YEAR") #Get a vector with all the years in which observations have been found for a certain species
        species = filename.split('.')[0] #Get the species name
        # Output a Shapefile for each different attribute
        for year in years:
            outSHP = outDir + species + "_" + year + u".shp"
            gp.Select_analysis (filename, outSHP, "\"C_YEAR\" = '" + year + "'")

del filename, years, species, outSHP
#END




# ### Predictor development ########################

#Standard Cell size for all final predictors = 100x100 m (maximum resolution of the species raster)
#"Old environments" ("old_envs")- Used for many processes:
    # Snap raster to Species_raster: spp_records_raster
    #Extent: Mittelland_all
#Extent and shape of all final predictors: Study_area


#  ### Predictors based on 'roads' (SwissTLM_Strasse) category of SwissTLM3D  (Swisstopo, 2016) ###
#Density of roads, Density of highways
        #Extract by attribute - Get only the relevant categories (road categories with significant car traffic, according to catalog of SwissTLM3D)
        #select
            #For roads: Set 5 (see Readme file: Data Preparation - Transport)
            #For highways: Set1 (see Readme file: Data Preparation - Transport)

# Local variables:
TLM_STRASSE = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\swissTLM3D_2016_LV03_LN02.gdb\\TLM_STRASSEN\\TLM_STRASSE"
Strasse_Set5_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Strasse_Set5_Mitt_all.shp"
Actual_roads_notunnels_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Actual_roads_notunnels_Mitt_all.shp"
# ### For Density of roads only
Roads_Mitt_all_25_tif = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Roads_Mitt_all_25.tif"
Rec_Roads_Mitt_all_25_tif = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Rec_Roads_Mitt_all_25.tif"
Road_density_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Road_density_25.tif"
Road_density_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Road_density_100_Mitt_all.tif"
Road_density_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Road_density_100_core.tif"

# Select All actual roads
arcpy.Select_analysis(TLM_STRASSE, Strasse_Set5_Mitt_all, "OBJEKTART IN ( 0, 1, 2, 3, 4, 5, 6, 8, 9, 10,11, 12, 20, 21)")
# Select (2) -	Exclude tunnels
arcpy.Select_analysis(Strasse_Set5_Mitt_all, Actual_roads_notunnels_Mitt_all, "KUNSTBAUTE NOT IN (1000)")

#Polyline to raster
            #To 25m (because it's a good resolution and an easy extrapolation/resampling to 100 m, apply old_envs, field value= Objektart
            #I expect it to auto-clip the raster to the extent of Mittelland_all by just specifying the environment
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PolylineToRaster_conversion(Actual_roads_notunnels_Mitt_all, "OBJEKTART", Roads_Mitt_all_25_tif, "MAXIMUM_LENGTH", "NONE", "25")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Reclassify
            #Set fixed looking at rasterizing value field options, as all the values will be 1, no data = 0
arcpy.gp.Reclassify_sa(Roads_Mitt_all_25_tif, "Value", "0 1;1 1;2 1;3 1;4 1;5 1;6 1;8 1;9 1;10 1;11 1;12 1;20 1;21 1;NODATA 0", Rec_Roads_Mitt_all_25_tif, "DATA")
#meaning of terms: (inputraster,valuefield(to be reclasified), "originalvalue1 newvalue1 (1); originalvalue2 newvalue2; originalvalueX newvalueX; originalvalueNoData newvalue0", outputRaster, "ChangemissingvaluestoNoData(Unchecked)")

        #Focal statistics
            #Mean, circular buffer: (2 Km)
arcpy.gp.FocalStatistics_sa(Rec_Roads_Mitt_all_25_tif, Road_density_25, "Circle 2000 MAP", "MEAN", "DATA")

        #Resample
            #To std Cell Size
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Road_density_25, Road_density_100_Mitt_all, "100", "BILINEAR")
#Meaning of terms: (input, output, cell size (same as), Resampling method (nearest neighbor)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

        #Clip raster (geometry)
            #To shape of Sampling_area
arcpy.Clip_management(Road_density_100_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Road_density_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
#Meaning of terms: (in_raster, rectangle, out_raster, {in_template_dataset}, {nodata_value}, {clipping_geometry}, {maintain_clipping_extent})


# ## Density of Highways ######################
#Local vars
Highways_set1_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Highways_set1_Mitt_all.shp"
Highways_Mitt_all_25_tif = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Highways_Mitt_all_25.tif"
Rec_Highways_Mitt_all_25_tif = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Rec_Highways_Mitt_all_25.tif"
Highway_density_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Highway_density_25.tif"
Highway_density_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Highway_density_100_Mitt_all.tif"
Highway_density_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\Highway_density_100_core.tif"

# Select (3) Highways and semi-highways only
arcpy.Select_analysis(Actual_roads_notunnels_Mitt_all, Highways_set1_Mitt_all, "OBJEKTART IN ( 0, 1, 2, 3, 5, 21)")


#Polyline to raster
            #To 25m (apply old_envs, field value= Objektart)
            #I expect it to auto-clip the raster to the extent of Mittelland_all by just specifying the environment
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PolylineToRaster_conversion(Highways_set1_Mitt_all, "OBJEKTART", Highways_Mitt_all_25_tif, "MAXIMUM_LENGTH", "NONE", "25")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Reclassify
            #Try to Set fixed, as all the values will be 1, no data = 0
arcpy.gp.Reclassify_sa(Highways_Mitt_all_25_tif, "VALUE", "0 1;1 1;2 1;3 1;5 1;21 1;NODATA 0", Rec_Highways_Mitt_all_25_tif, "DATA")
#meaning of terms: (inputraster,valuefield(to be reclasified), "originalvalue1 newvalue1 (1); originalvalue2 newvalue2; originalvalueX newvalueX; originalvalueNoData newvalue0", outputRaster, "ChangemissingvaluestoNoData(Unchecked)")

        #Focal statistics
            #Mean, circular buffer: Max. dispersal dist. of the sp. (2 Km)
arcpy.gp.FocalStatistics_sa(Rec_Highways_Mitt_all_25_tif, Highway_density_25, "Circle 2000 MAP", "MEAN", "DATA")

        #Resample
            #To std Cell Size
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Highway_density_25, Highway_density_100_Mitt_all, "100", "BILINEAR")
#Meaning of terms: (input, output, cell size (same as), Resampling method (nearest neighbor)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

        #Clip raster (geometry)
            #To shape of Sampling_area
arcpy.Clip_management(Highway_density_100_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Highway_density_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# #########################################################
# ### Predictors based on 'railways' category of SwissTLM3D #####
    #Same procedures, one less selection
# Local variables:
TLM_RAILS = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\swissTLM3D_2016_LV03_LN02.gdb\\TLM_OEV\\TLM_EISENBAHN"
Rails_notunnels_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_rails\\Rails_notunnels_Mitt_all.shp"
Rails_Mitt_all_25_tif = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_rails\\Rails_Mitt_all_25.tif"
Rec_Rails_Mitt_all_25_tif = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_rails\\Rec_Rails_Mitt_all_25.tif"
Rail_density_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_rails\\Rail_density_25.tif"
Rail_density_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_rails\\Rail_density_100_Mitt_all.tif"
Rail_density_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_rails\\Rail_density_100_core.tif"

# Exclude tunnels
arcpy.Select_analysis(TLM_RAILS, Rails_notunnels_Mitt_all, "KUNSTBAUTE NOT IN (800)")

#Polyline to raster
            #To 25m (because it's a good resolution and an easy extrapolation/resampling to 100 m, apply old_envs, field value= Objektart
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PolylineToRaster_conversion(Rails_notunnels_Mitt_all, "KUNSTBAUTE", Rails_Mitt_all_25_tif, "MAXIMUM_LENGTH", "NONE", "25")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Reclassify
            #all the values will be 1, no data = 0
arcpy.gp.Reclassify_sa(Rails_Mitt_all_25_tif, "VALUE", "100 1;200 1;300 1;400 1;900 1;NODATA 0", Rec_Rails_Mitt_all_25_tif, "DATA")


#Focal statistics
            #Mean, circular buffer: Max. dispersal dist. of the sp. (2 Km)
arcpy.gp.FocalStatistics_sa(Rec_Rails_Mitt_all_25_tif, Rail_density_25, "Circle 2000 MAP", "MEAN", "DATA")

#Resample
            #To std Cell Size
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Rail_density_25, Rail_density_100_Mitt_all, "100", "BILINEAR")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

        #Clip raster (geometry)
            #To shape of Study_area
arcpy.Clip_management(Rail_density_100_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Rail_density_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# ############################################################
# Recreation intensity
# Local variables:
inix_10 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_recreation\\inix10.tif" #Recreation intensity, from Auswertungsprotokoll f?r Parameter 31b Potenzielle Naherholungsgebiete um Siedlungen (BAFU, 2012)
inix10_res100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_recreation\\inix10_res100.tif"
inix10_res100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_recreation\\inix10_res100_core.tif"

#Resample
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(inix_10, inix10_res100, "100", "BILINEAR")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Clip raster (geometry)#To shape of study area
arcpy.Clip_management(inix10_res100, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", inix10_res100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")


# ############################################################
# ### Population density - From STATPOP (BFS 2015)
#  Local variables
#Raw layer input
STATPOP15_whole_CH = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Population_density\\STATPOP2015B_TOT_CORR.shp"
#Rasterized feature (to 100m resolution)
STATPOP15_Mitt_all_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Population_density\\STATPOP15_Mitt_all_100.tif"
STATPOP15_Mitt_all_100_noNoData = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Population_density\\STATPOP15_Mitt_all_100_noNoData.tif"
#focal stats output
Population_Density_15 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Population_density\\Population_Density_15.tif"
#Final predictor
Population_Density_15_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Population_density\\Population_Density_15_core.tif"

#Point to Raster
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PointToRaster_conversion(STATPOP15_whole_CH, "B15BTOT", STATPOP15_Mitt_all_100, "MEAN", "NONE", spp_records_raster)
#Meaning of things inside parenthesis:(inputFeatures, valueField, outputRaster, assignmentType, priorityField, cellSize)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Change all NoData to 0
STATPOP15_Mitt_all_100_noNoData = Con(IsNull(STATPOP15_Mitt_all_100),0,STATPOP15_Mitt_all_100)
# Save the output
STATPOP15_Mitt_all_100_noNoData.save("E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Population_density\\STATPOP15_Mitt_all_100_noNoData.tif")

#Focal Statistics
arcpy.gp.FocalStatistics_sa(STATPOP15_Mitt_all_100_noNoData, Population_Density_15, "Circle 2000 MAP", "MEAN", "DATA")
# Meaning of terms inside parenthesis: (input raster, output raster, "neighborhood_shape radius_of_shape map_units(m)", "calculated statistic", "Ignore NoData in calculations: checked")

#Clip raster (geometry)#To shape of study area
arcpy.Clip_management(Population_Density_15, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Population_Density_15_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
#Meaning of terms: (in_raster, rectangle, out_raster, {in_template_dataset}, {nodata_value}, {clipping_geometry}, {maintain_clipping_extent})


# ############################################################
# Forest-related variables - From Waldmischungsgrad (BFS 2013)
    #Density of forest
    #At-site categorical binary vars: deciduous, mixed and coniferous forest

# Local vars:
Forest_mixture_CH = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\wmg25"
Forest_mit_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Forest_mit_all.tif"

# Clip to all_Mittelland extent
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
arcpy.Clip_management(Forest_mixture_CH, "485410 109645 768768 295145", Forest_mit_all, Mittelland_all, "", "NONE", "NO_MAINTAIN_EXTENT")
arcpy.env.snapRaster = tempEnvironment0

#Reclassify to produce 4 different rasters
#All values with relevant category of forest=1, all others=0
#local vars for the rest of the forest rasters:
Rec_Forest_25_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Forest_nonforest_25_Mitt_all.tif"
Deciduous_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Deciduous_25_Mitt_all.tif"
Coniferous_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Coniferous_25_Mitt_all.tif"
Mixed_forest_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Mixed_forest_25_Mitt_all.tif"

Forest_density_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Forest_density_25.tif"

Forest_density_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Forest_density_100_Mitt_all.tif"
Deciduous_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Deciduous_100_Mitt_all.tif"
Coniferous_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Coniferous_100_Mitt_all.tif"
Mixed_forest_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Mixed_forest_100_Mitt_all.tif"

Forest_density_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Forest_density_100_core.tif"
Deciduous_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Deciduous_100_core.tif"
Coniferous_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Coniferous_100_core.tif"
Mixed_forest_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Mixed_forest_100_core.tif"

    # Reclass 1 - Into Forest & non-forest
arcpy.gp.Reclassify_sa(Forest_mit_all, "Value", "0 0;1 1;2 1;3 1;4 1;9 1;NODATA 0", Rec_Forest_25_Mitt_all, "DATA")
    #meaning of terms: (inputraster,valuefield(to be reclasified), "originalvalueX newvalueX; originalvalueNoData newvalue0", outputRaster, "ChangemissingvaluestoNoData(Unchecked)")
    #Reclassify 2 -	Deciduous
arcpy.gp.Reclassify_sa(Forest_mit_all, "Value", "0 0;1 0;2 0;3 0;4 1;9 0;NODATA 0", Deciduous_25, "DATA")
#Reclassify 3 -	Coniferous
arcpy.gp.Reclassify_sa(Forest_mit_all, "Value", "0 0;1 1;2 0;3 0;4 0;9 0;NODATA 0", Coniferous_25, "DATA")
#Reclassify 4 - Mixed forest (join categories 2 and 3)
arcpy.gp.Reclassify_sa(Forest_mit_all, "Value", "0 0;1 0;2 1;3 1;4 0;9 0;NODATA 0", Mixed_forest_25, "DATA")

# Focal statistics (Only to Rec_Forest_25_Mitt_all)
arcpy.gp.FocalStatistics_sa(Rec_Forest_25_Mitt_all, Forest_density_25, "Circle 2000 MAP", "MEAN", "DATA")


# Resample (All 4 forest layers) To std Cell Size
# At-site categorical rasters with nearest neighbor, density with bilinear interpolation)
    #Resample 1: Forest_density
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Forest_density_25, Forest_density_100_Mitt_all, "100", "BILINEAR")
#Meaning of terms: (input, output, cell size (same as), Resampling method)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1
    #Resample 2: Deciduous
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Deciduous_25, Deciduous_100_Mitt_all, "100", "NEAREST")
#Meaning of terms: (input, output, cell size (same as), Resampling method)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1
    #Resample 3: Coniferous
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Coniferous_25, Coniferous_100_Mitt_all, "100", "NEAREST")
#Meaning of terms: (input, output, cell size (same as), Resampling method)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1
    #Resample 4: Mixed
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Mixed_forest_25, Mixed_forest_100_Mitt_all, "100", "NEAREST")
#Meaning of terms: (input, output, cell size (same as), Resampling method)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

# Clip geometry - To shape of Study area
    #Clip1: forest density
arcpy.Clip_management(Forest_density_100_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Forest_density_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
#Meaning of terms: (in_raster, rectangle, out_raster, {in_template_dataset}, {nodata_value}, {clipping_geometry}, {maintain_clipping_extent})
    #Clip2: At site Deciduous status (0/1)
arcpy.Clip_management(Deciduous_100_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Deciduous_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Clip3: At site Coniferous status (0/1)
arcpy.Clip_management(Coniferous_100_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Coniferous_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
 #Clip4: At site Mixed_forest status (0/1)
arcpy.Clip_management(Mixed_forest_100_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Mixed_forest_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")


# ##  Distance to forest edge #####
#New vars for this
forest_edge_polygons_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\forest_edge_polygons_Mitt_all.shp"
forest_edge_polyline_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\forest_edge_polyline_Mitt_all.shp"
Dist_to_forest_edge_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Dist_to_forest_edge_100_Mitt_all.tif"
Dist_to_forest_edge_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Dist_to_forest_edge_100_core.tif"
DirectionRaster_Dist_to_forest_edge = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_forest\\Direction_forest_edge_100.tif"


#Raster to Polygon
# raster must be integer type
arcpy.RasterToPolygon_conversion (Rec_Forest_25_Mitt_all, forest_edge_polygons_Mitt_all, "NO_SIMPLIFY", "Value")

#Polygon to Polyline
arcpy.PolygonToLine_management (forest_edge_polygons_Mitt_all,forest_edge_polyline_Mitt_all)

#Euclidean distance (from the lines) - Set to Old envs
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.gp.EucDistance_sa(forest_edge_polyline_Mitt_all, Dist_to_forest_edge_Mitt_all, "", spp_records_raster, DirectionRaster_Dist_to_forest_edge)
#(input, outputdistance, maxDistance, cellSize, outDirectionRaster)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Clip geometry
arcpy.Clip_management(Dist_to_forest_edge_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Dist_to_forest_edge_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# ############################################################
#Climate variables
    #Mean annual precipitation
    #Mean summer precipitation
    #Mean annual direct solar radiation
    #Mean annual temperature
#local vars
Temp_avg_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Taveyy_03.tif"
Precip_year_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Precyy_03.tif"
Precip_sum_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Precsu_03.tif"
Solardir_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Sdiryy_03.tif"

Temp_avg_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Temp_avg_100.tif"
Precip_year_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Precip_year_100.tif"
Precip_sum_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Precip_sum_100.tif"
Solardir_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Solardir_100.tif"

Temp_avg_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Temp_avg_100_core.tif"
Precip_year_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Precip_year_100_core.tif"
Precip_sum_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Precip_sum_100_core.tif"
Solardir_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_climate\\Solardir_100_core.tif"

#Resample (all the climatic layers)
#With old envs
    #Resample 1 - temp_avg_yy
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Temp_avg_25, Temp_avg_100, "100", "BILINEAR")
#Meaning of terms: (input, output, cell size (same as), Resampling method)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

    #Resample 2 - Precip_year
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Precip_year_25, Precip_year_100, "100", "BILINEAR")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

    #Resample 3 - Precip_summer
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Precip_sum_25, Precip_sum_100, "100", "BILINEAR")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

    #Resample 4 - direct solar radiation
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Solardir_25, Solardir_100, "100", "BILINEAR")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Clip geometry (all the climatic layers)
    #Clip geometry 1 - temp_avg_yy
arcpy.Clip_management(Temp_avg_100, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Temp_avg_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Clip geometry 2 - Precip_year
arcpy.Clip_management(Precip_year_100, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Precip_year_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Clip geometry 3 - Precip_summer
arcpy.Clip_management(Precip_sum_100, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Precip_sum_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Clip geometry 4 - direct solar radiation
arcpy.Clip_management(Solardir_100, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Solardir_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")



# ############################################################
# Density of traffic
#To get raster 'Traffic_CH_raster10' out of the NPVM shp with all the traffic data of CH:
    #1.Open 'DTV_SZ_MIV_2010_KGM_1250_NEU_GV+Netz2010_reduziert_VSt_Basis_Kal_link', that's the layer of MIV_2010 based on 2010 calibration.
    #2.	Define projection, to CH1903_LV03
    #3.	Start editing, add new field: TOT_TRAFF, short integer, precision: 5 (Settings for the Traffic columns (VOLVEHPR~1, R_VOLVEH~7: Gesamtbelastung Fz/Tag (Total traffic load in each direction)
    #4.	Populate the field with field calculator Sum TrafficDir1+TrafficDir2

#Take out tunnels
        #Defined categories in NPVM shp properties (TYPE_NO field): '0', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '40', '41', '42', '57', '58', '60', '61', '62', '63', '80', '81', '90', '91'

#Local vars:
Traffic_CH_bothdirections = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\DTV_SZ_MIV_2010_KGM_1250_NEU_GV+Netz2010_reduziert_VSt_Basis_Kal_link.shp"
Traffic_CH_bothdirs_notunnels = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_CH_bothdirs_notunnels.shp"

Traffic_CH_raster10 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_CH_raster10_notunnels.tif"

Traffic_CH_raster10_IsNull = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_CH_raster10_IsNull.tif"
Traffic_CH_raster10_noNoData = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_CH_raster10_noNoData.tif"

Traffic_Density_10 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_Density_10_notunnels.tif"
Traffic_Density_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_Density_100_notunnels.tif"
Traffic_Density_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_Density_100_core_notunnels.tif"

# Select  -	Exclude tunnels
arcpy.Select_analysis(Traffic_CH_bothdirections, Traffic_CH_bothdirs_notunnels, "TYPE_NO IN ('0', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '40', '41', '42', '57', '58', '60', '61', '62', '63', '80', '81', '90', '91')")

# Perform Polyline to Raster
        #Output: Traffic_CH_raster10.tif
            #To 10m
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PolylineToRaster_conversion(Traffic_CH_bothdirs_notunnels, "TOT_TRAFF", Traffic_CH_raster10, "MAXIMUM_LENGTH", "NONE", "10")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Change all NoData to 0
#Embed IsNull in a Con function
# outCon = Con(IsNull(inputRaster),0,inputRaster)
Traffic_CH_raster10_noNoData= Con(IsNull(Traffic_CH_raster10),0,Traffic_CH_raster10)
# Save the output
Traffic_CH_raster10_noNoData.save("E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_CH_raster10_noNoData.tif")

# Focal Statistics
arcpy.gp.FocalStatistics_sa(Traffic_CH_raster10_noNoData, Traffic_Density_10, "Circle 2000 MAP", "MEAN", "DATA")

# Resample
    # old envs
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Traffic_Density_10, Traffic_Density_100, "100", "BILINEAR")
#Meaning of terms: (input, output, cell size (same as), Resampling method (nearest neighbor)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

# Clip geometry
arcpy.Clip_management(Traffic_Density_100, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Traffic_Density_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")


# ############################################################
# Noise - From EMPA (2011)
#Raw rasters at 10x10 resolution
#Get one raster for Total noise at daytime and another for total noise at nighttime

#local vars
Strassenlaerm_Tag = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Strassenlaerm_Tag.tif"
Strassenlaerm_Nacht = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Strassenlaerm_Nacht.tif"
Eisenbahnlaerm_Tag = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Eisenbahnlaerm_Tag.tif"
Eisenbahnlaerm_Nacht = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Eisenbahnlaerm_Nacht.tif"

Totalnoise_daytime_10 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_daytime_10.tif"
Totalnoise_nighttime_10 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_nighttime_10.tif"

Totalnoise_daytime_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_daytime_100.tif"
Totalnoise_nighttime_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_nighttime_100.tif"

Totalnoise_daytime_100_noNoData = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_daytime_100_noNoData.tif"
Totalnoise_nighttime_100_noNoData = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_daytime_100_noNoData.tif"

Totalnoise_daytime_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_daytime_100_core.tif"
Totalnoise_nighttime_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_nighttime_100_core.tif"

#Cell statistics Train+road noise
    #Day
arcpy.gp.CellStatistics_sa("E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Strassenlaerm_Tag.tif;E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Eisenbahnlaerm_Tag.tif", Totalnoise_daytime_10, "SUM", "DATA")
#(in_raster_or_constant1;in_rasters_or_constant2, outRaster, {statistics_type}, {ignore_nodata})
    #Night
arcpy.gp.CellStatistics_sa("E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Strassenlaerm_Nacht.tif;E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Eisenbahnlaerm_Nacht.tif", Totalnoise_nighttime_10, "SUM", "DATA")

#Resample
#Use Old envs
    #Day
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Totalnoise_daytime_10, Totalnoise_daytime_100, "100", "BILINEAR")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1
    #Night
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Totalnoise_nighttime_10, Totalnoise_nighttime_100, "100", "BILINEAR")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Change all NoData to 0
    #Day
Totalnoise_daytime_100_noNoData= Con(IsNull(Totalnoise_daytime_100),0,Totalnoise_daytime_100)
# Save the output
Totalnoise_daytime_100_noNoData.save("E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_daytime_100_noNoData.tif")
    #Night
Totalnoise_nighttime_100_noNoData= Con(IsNull(Totalnoise_nighttime_100),0,Totalnoise_nighttime_100)
# Save the output
Totalnoise_nighttime_100_noNoData.save("E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_noise\\Totalnoise_nighttime_100_noNoData.tif")

# Clip geometry
    #Day
arcpy.Clip_management(Totalnoise_daytime_100_noNoData, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Totalnoise_daytime_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Night
arcpy.Clip_management(Totalnoise_nighttime_100_noNoData, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Totalnoise_nighttime_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# ############################################################
#Land cover (Base categories) - From Arealstatistik (OFS 2010)
# transformed into:
# Agriculture density
# Presence of Agriculture categories
# Presence of Human settlement categories

# Local variables:
Arealstatistik_points = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\XY_AREA_NOAS04_72_85_09.shp"
Arealstatistik_raster_CH = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Arealstatistik_raster_CH.tif"
#Agro Classified by all 72 categories of 2004-2009
Agriculture_only_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Agriculture_only_Mitt_all.tif"
Orchard_vineyard_hort_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Orchard_vineyard_hort_Mitt_all.tif"
Arable_land_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Arable_land_Mitt_all.tif"
Meadows_farmpastures_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Meadows_farmpastures_Mitt_all.tif"
#settlements
Grey_settlements_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Grey_settlements_Mitt_all.tif"
Green_settlements_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Green_settlements_Mitt_all.tif"

Agriculture_density_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Agriculture_density_Mitt_all.tif"

Agriculture_density_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Agriculture_density_core.tif"
Orchard_vineyard_hort_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Orchard_vineyard_hort_core.tif"
Arable_land_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Arable_land_core.tif"
Meadows_farmpastures_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Meadows_farmpastures_core.tif"

Grey_settlements_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Grey_settlements_core.tif"
Green_settlements_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_landuse_cover\\Green_settlements_core.tif"

# Process: Point to Raster, follow old envs
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PointToRaster_conversion(Arealstatistik_points, "as09_72", Arealstatistik_raster_CH, "MOST_FREQUENT", "NONE", spp_records_raster)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Reclassifications to get all of the groups we need
 # Reclass1: All_agro
arcpy.gp.Reclassify_sa(Arealstatistik_raster_CH, "Value", "1 0;2 0;3 0;4 0;5 0;6 0;7 0;8 0;9 0;10 0;11 0;12 0;13 0;14 0;15 0;16 0;17 0;18 0;19 0;20 0;21 0;22 0;23 0;24 0;25 0;26 0;27 0;28 0;29 0;30 0;31 0;32 0;33 0;34 0;35 0;36 0;37 1;38 1;39 1;40 1;41 1;42 1;43 1;44 1;45 1;46 1;47 1;48 1;49 1;50 0;51 0;52 0;53 0;54 0;55 0;56 0;57 0;58 0;59 0;60 0;61 0;62 0;63 0;64 0;65 0;66 0;67 0;68 0;69 0;70 0;71 0;72 0;NODATA 0", Agriculture_only_Mitt_all, "DATA")
#meaning of terms: (inputraster,valuefield(to be reclasified), "originalvalueX newvalueX; originalvalueNoData newvalue0", outputRaster, "ChangemissingvaluestoNoData(Unchecked)")
    # Reclass2: Orchard, vineyard and horticulture areas
arcpy.gp.Reclassify_sa(Arealstatistik_raster_CH, "Value", "1 0;2 0;3 0;4 0;5 0;6 0;7 0;8 0;9 0;10 0;11 0;12 0;13 0;14 0;15 0;16 0;17 0;18 0;19 0;20 0;21 0;22 0;23 0;24 0;25 0;26 0;27 0;28 0;29 0;30 0;31 0;32 0;33 0;34 0;35 0;36 0;37 1;38 1;39 1;40 1;41 0;42 0;43 0;44 0;45 0;46 0;47 0;48 0;49 0;50 0;51 0;52 0;53 0;54 0;55 0;56 0;57 0;58 0;59 0;60 0;61 0;62 0;63 0;64 0;65 0;66 0;67 0;68 0;69 0;70 0;71 0;72 0;NODATA 0", Orchard_vineyard_hort_Mitt_all, "DATA")
    # Reclass3 - arable land
arcpy.gp.Reclassify_sa(Arealstatistik_raster_CH, "Value", "1 0;2 0;3 0;4 0;5 0;6 0;7 0;8 0;9 0;10 0;11 0;12 0;13 0;14 0;15 0;16 0;17 0;18 0;19 0;20 0;21 0;22 0;23 0;24 0;25 0;26 0;27 0;28 0;29 0;30 0;31 0;32 0;33 0;34 0;35 0;36 0;37 0;38 0;39 0;40 0;41 1;42 0;43 0;44 0;45 0;46 0;47 0;48 0;49 0;50 0;51 0;52 0;53 0;54 0;55 0;56 0;57 0;58 0;59 0;60 0;61 0;62 0;63 0;64 0;65 0;66 0;67 0;68 0;69 0;70 0;71 0;72 0;NODATA 0", Arable_land_Mitt_all, "DATA")
    # Reclass4 - Meadows_farmpastures
arcpy.gp.Reclassify_sa(Arealstatistik_raster_CH, "Value", "1 0;2 0;3 0;4 0;5 0;6 0;7 0;8 0;9 0;10 0;11 0;12 0;13 0;14 0;15 0;16 0;17 0;18 0;19 0;20 0;21 0;22 0;23 0;24 0;25 0;26 0;27 0;28 0;29 0;30 0;31 0;32 0;33 0;34 0;35 0;36 0;37 0;38 0;39 0;40 0;41 0;42 1;43 1;44 1;45 0;46 0;47 0;48 0;49 0;50 0;51 0;52 0;53 0;54 0;55 0;56 0;57 0;58 0;59 0;60 0;61 0;62 0;63 0;64 0;65 0;66 0;67 0;68 0;69 0;70 0;71 0;72 0;NODATA 0", Meadows_farmpastures_Mitt_all, "DATA")

 # Reclass5 - Grey_settlements (1-30,32)
arcpy.gp.Reclassify_sa(Arealstatistik_raster_CH, "Value", "1 1;2 1;3 1;4 1;5 1;6 1;7 1;8 1;9 1;10 1;11 1;12 1;13 1;14 1;15 1;16 1;17 1;18 1;19 1;20 1;21 1;22 1;23 1;24 1;25 1;26 1;27 1;28 1;29 1;30 1;31 0;32 1;33 0;34 0;35 0;36 0;37 0;38 0;39 0;40 0;41 0;42 0;43 0;44 0;45 0;46 0;47 0;48 0;49 0;50 0;51 0;52 0;53 0;54 0;55 0;56 0;57 0;58 0;59 0;60 0;61 0;62 0;63 0;64 0;65 0;66 0;67 0;68 0;69 0;70 0;71 0;72 0;NODATA 0", Grey_settlements_Mitt_all, "DATA")

 # Reclass6 - Green_settlements (31, 33-36)
arcpy.gp.Reclassify_sa(Arealstatistik_raster_CH, "Value", "1 0;2 0;3 0;4 0;5 0;6 0;7 0;8 0;9 0;10 0;11 0;12 0;13 0;14 0;15 0;16 0;17 0;18 0;19 0;20 0;21 0;22 0;23 0;24 0;25 0;26 0;27 0;28 0;29 0;30 0;31 1;32 0;33 1;34 1;35 1;36 1;37 0;38 0;39 0;40 0;41 0;42 0;43 0;44 0;45 0;46 0;47 0;48 0;49 0;50 0;51 0;52 0;53 0;54 0;55 0;56 0;57 0;58 0;59 0;60 0;61 0;62 0;63 0;64 0;65 0;66 0;67 0;68 0;69 0;70 0;71 0;72 0;NODATA 0", Green_settlements_Mitt_all, "DATA")

#Focal statistics - only for Agro_areas_all
arcpy.gp.FocalStatistics_sa(Agriculture_only_Mitt_all, Agriculture_density_Mitt_all, "Circle 2000 MAP", "MEAN", "DATA")

#Clip geometry
    #Agriculture_density
arcpy.Clip_management(Agriculture_density_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Agriculture_density_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Orchard, vineyard and horticulture areas
arcpy.Clip_management(Orchard_vineyard_hort_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Orchard_vineyard_hort_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #arable land
arcpy.Clip_management(Arable_land_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Arable_land_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Meadows_farmpastures
arcpy.Clip_management(Meadows_farmpastures_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Meadows_farmpastures_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Grey_settlements
arcpy.Clip_management(Grey_settlements_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Grey_settlements_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")
    #Green_settlements
arcpy.Clip_management(Green_settlements_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Green_settlements_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")


# ############################################################
# Presence of Rivers
#Local vars
Rivers_CH_shp = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\swissTLM3D_2016_LV03_LN02.gdb\\TLM_GEWAESSER\\TLM_FLIESSGEWAESSER"
#Rivers_CH_Noimaginarylines = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_water_bodies\\Rivers_CH_Noimaginarylines.shp"
Rivers_CH_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_water_bodies\\Rivers_CH_Mitt_all.tif"
Rec_Rivers_CH_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_water_bodies\\Rec_Rivers_CH_Mitt_all.tif"
Rivers_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_water_bodies\\Rivers_core.tif"

#Polyline to raster
            #To 100m (apply old_envs, field value= Objektart)
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PolylineToRaster_conversion(Rivers_CH_shp, "OBJEKTART", Rivers_CH_Mitt_all, "MAXIMUM_LENGTH", "NONE", "100")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Reclassify
arcpy.gp.Reclassify_sa(Rivers_CH_Mitt_all, "Value", "1 1;2 1;3 1;4 1;5 1;6 0;NODATA 0", Rec_Rivers_CH_Mitt_all, "DATA")

#clip geometry
arcpy.Clip_management(Rec_Rivers_CH_Mitt_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Rivers_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")

# ############################################################
#Slope
#vars
DEM_CH_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_elevation\\dhm25grid.tif"
Slope_CH_25 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_elevation\\Slope_CH_25.tif"
Slope_CH_100 = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_elevation\\Slope_CH_100.tif"
Slope_100_core = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_elevation\\Slope_CH_100_core.tif"

#Calculate slope
#outSlope = Slope(inRaster, outMeasurement, zFactor)
# Save the output
#outSlope.save("C:/sapyexamples/output/outslope02")
Slope_CH_25 = arcpy.sa.Slope(DEM_CH_25, "DEGREE", "")
Slope_CH_25.save("E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_elevation\\Slope_CH_25.tif")

#Resample
    #Use OldEnvs
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.Resample_management(Slope_CH_25, Slope_CH_100, "100", "BILINEAR")

#Clip geometry
arcpy.Clip_management(Slope_CH_100, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Slope_100_core, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")


# ######### Doing mask ################################################################
#All amphibian ocurrences + flachmoore + lakes + amphibian spawning sites + gravel pits

# ##All amphibian ocurrences
#When running complete script, use 1st argument (simple equality), when running by parts, excluding the species raster setup, define object
#outSHP = Presence_points_core_no_sal
Presence_points_core_no_sal = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Presence_points_core_no_sal.shp"
Presence_points_no_sal_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Presence_points_no_sal_100_Mitt_all.tif"
Binary_Presence_points_no_sal_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Binary_Presence_points_no_sal_100_Mitt_all.tif"
#Polygon to raster
    #old envs,100 cell size
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PointToRaster_conversion(Presence_points_core_no_sal, "ART", Presence_points_no_sal_100_Mitt_all, "MOST_FREQUENT", "NONE", "100")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1
#Reclassify, to binary map
arcpy.gp.Reclassify_sa(Presence_points_no_sal_100_Mitt_all, "ART", "' ' 0;'Bombina variegata' 1;'Bufo bufo' 1;'Epidalea  calamita' 1;'Hyla arborea' 1;'Ichthyosaura alpestris' 1;'Lissotriton helveticus' 1;'Pelophylax aggr.' 1;'Pelophylax ridibundus' 1;'Rana dalmatina' 1;'Rana temporaria' 1;'Triturus cristatus' 1;'Alytes obstetricans' 1;NODATA 0", Binary_Presence_points_no_sal_100_Mitt_all, "DATA")

# ### Lakes
Lakes_all_CH_shp = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_roads\\swissTLM3D_2016_LV03_LN02.gdb\\TLM_GEWAESSER\\TLM_STEHENDES_GEWAESSER"
Lakes_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Lakes_100_Mitt_all.tif"
Binary_Lakes_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Binary_Lakes_100_Mitt_all.tif"
#Feature to raster
    #OldEnvs
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.FeatureToRaster_conversion(Lakes_all_CH_shp, "OBJEKTART", Lakes_100_Mitt_all, "100")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1
#Reclassify to binary map
arcpy.gp.Reclassify_sa(Lakes_100_Mitt_all, "Value", "0 1;1 1;NODATA 0", Binary_Lakes_100_Mitt_all, "DATA")

# ### Swamps (Flachmoore)
Flachmoore_shp = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\fm.shp"
Flachmoore_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Flachmoore_100_Mitt_all.tif"
Binary_Flachmoore_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Binary_Flachmoore_100_Mitt_all.tif"
#Polygon to raster
    #OldEnvs
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PolygonToRaster_conversion(Flachmoore_shp, "FM_FL", Flachmoore_100_Mitt_all, "CELL_CENTER", "NONE", "100")
#Meaning of things inside parenthesis:(inputFeatures, valueField, outputRaster, assignmentType, priorityField, cellSize)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1
#Reclassify to binary map
arcpy.gp.Reclassify_sa(Flachmoore_100_Mitt_all, "VALUE", "0.066886 5.756350 1;5.756350 12.573606 1;12.573606 21.262525 1;21.262525 31.159669 1;31.159669 43.506616 1;43.506616 61.289550 1;61.289550 85.332813 1;85.332813 141.908488 1;141.908488 197.323000 1;NODATA 0", Binary_Flachmoore_100_Mitt_all, "DATA")

# ### Amphibian Spawning sites
AmphSpawnSites_shp = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\amphibien_l\\polygon"
AmphSpawnSites_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\AmphSpawnSites_raster100.tif"
Binary_AmphSpawnSites_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Binary_AmphSpawnSites_raster100.tif"
#Polygon to raster
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.PolygonToRaster_conversion(AmphSpawnSites_shp, "AM_L_BEREICH", AmphSpawnSites_100_Mitt_all, "CELL_CENTER", "NONE", "100")
#Meaning of things inside parenthesis:(inputFeatures, valueField, outputRaster, assignmentType, priorityField, cellSize)
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

#Reclassify to binary map
    #Use Old Envs
tempEnvironment0 = arcpy.env.snapRaster
arcpy.env.snapRaster = spp_records_raster
tempEnvironment1 = arcpy.env.extent
arcpy.env.extent = Mittelland_all
arcpy.gp.Reclassify_sa(AmphSpawnSites_100_Mitt_all, "Value", "1 1;2 1;3 1;NODATA 0", Binary_AmphSpawnSites_100_Mitt_all, "DATA")
arcpy.env.snapRaster = tempEnvironment0
arcpy.env.extent = tempEnvironment1

# ### # Join rasterized versions of Lakes, Flachmoore, and Amphibian Spawning Sites with tool 'Mosaic to New Raster'
    # Pixel depth (bit depth) has to be the same for all inputs
    #Copy rasters 1st to Make all layers the same depth (8-bit)
    #Binary_Flachmoore_100_Mitt_all already has this depth
#local vars
EightBit_Binary_Presence_points_no_sal_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\EightBit_Binary_Presence_points_no_sal_100_Mitt_all.tif"
EightBit_Binary_Lakes_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\EightBit_Binary_Lakes_100_Mitt_all.tif"
EightBit_Binary_AmphSpawnSites_100_Mitt_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\EightBit_Binary_AmphSpawnSites_100_Mitt_all.tif"

Mask_directory = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask"
Amph_habs_Mit_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Amph_habs_Mit_all.tif"
Binary_Amph_habs_Mit_all = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Binary_Amph_habs_Mit_all.tif"
Mask = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Mask_Amph_habs.tif"

# Copy Rasters
    #Binary_Presence_points_no_sal_100_Mitt_all
arcpy.CopyRaster_management(Binary_Presence_points_no_sal_100_Mitt_all, EightBit_Binary_Presence_points_no_sal_100_Mitt_all, "", "", "255", "NONE", "ColormapToRGB", "8_BIT_UNSIGNED", "NONE", "NONE", "TIFF", "NONE")
    #Binary_Lakes_100_Mitt_all
arcpy.CopyRaster_management(Binary_Lakes_100_Mitt_all, EightBit_Binary_Lakes_100_Mitt_all, "", "", "255", "NONE", "ColormapToRGB", "8_BIT_UNSIGNED", "NONE", "NONE", "TIFF", "NONE")
    #Binary_AmphSpawnSites_100_Mitt_all
arcpy.CopyRaster_management(Binary_AmphSpawnSites_100_Mitt_all, EightBit_Binary_AmphSpawnSites_100_Mitt_all, "", "", "255", "NONE", "ColormapToRGB", "8_BIT_UNSIGNED", "NONE", "NONE", "TIFF", "NONE")

# Process: Mosaic To New Raster
    #Very important: Mosaic Operator: Sum (Otherwise the nodata and 0 values overlap the actual '1's
#arcpy.MosaicToNewRaster_management(Binary_Presence_points_no_sal_100_Mitt_all;Binary_Lakes_100_Mitt_all;Binary_AmphSpawnSites_100_Mitt_all;Binary_Flachmoore_100_Mitt_all, Mask_directory, Amph_habs_Mit_all, "PROJCS['CH1903_LV03',GEOGCS['GCS_CH1903',DATUM['D_CH1903',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Hotine_Oblique_Mercator_Azimuth_Center'],PARAMETER['False_Easting',600000.0],PARAMETER['False_Northing',200000.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Azimuth',90.0],PARAMETER['Longitude_Of_Center',7.439583333333333],PARAMETER['Latitude_Of_Center',46.95240555555556],UNIT['Meter',1.0]]", "8_BIT_UNSIGNED", "100", "1", "SUM", "FIRST")
arcpy.MosaicToNewRaster_management([EightBit_Binary_Presence_points_no_sal_100_Mitt_all,EightBit_Binary_Lakes_100_Mitt_all,EightBit_Binary_AmphSpawnSites_100_Mitt_all,Binary_Flachmoore_100_Mitt_all], Mask_directory, "Amph_habs_Mit_all.tif", "", "8_BIT_UNSIGNED", "", "1", "SUM", "FIRST")

#Reclassify to binary map
arcpy.gp.Reclassify_sa(Amph_habs_Mit_all, "Value", "0 0;1 1;2 1;3 1;4 1;NODATA 0", Binary_Amph_habs_Mit_all, "DATA")

#clip geometry
arcpy.Clip_management(Binary_Amph_habs_Mit_all, "489594.481512878 112941.234701329 766238.240258206 289748.081452074", Mask, Study_area, "", "ClippingGeometry", "NO_MAINTAIN_EXTENT")


### ######### Defining projection of HSM outputs: In Script 'Mask application + Patch (node) delineation' ################################################################
##
###Define projection of HSM maps - In ArcGIS
##    #HSM
##arcpy.DefineProjection_management(Auto2ndHSMap, "PROJCS['CH1903_LV03',GEOGCS['GCS_CH1903',DATUM['D_CH1903',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Hotine_Oblique_Mercator_Azimuth_Center'],PARAMETER['False_Easting',600000.0],PARAMETER['False_Northing',200000.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Azimuth',90.0],PARAMETER['Longitude_Of_Center',7.439583333333333],PARAMETER['Latitude_Of_Center',46.95240555555556],UNIT['Meter',1.0]]")
##    #ROCBin
##arcpy.DefineProjection_management(ROCBin_2ndHSMap, "PROJCS['CH1903_LV03',GEOGCS['GCS_CH1903',DATUM['D_CH1903',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Hotine_Oblique_Mercator_Azimuth_Center'],PARAMETER['False_Easting',600000.0],PARAMETER['False_Northing',200000.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Azimuth',90.0],PARAMETER['Longitude_Of_Center',7.439583333333333],PARAMETER['Latitude_Of_Center',46.95240555555556],UNIT['Meter',1.0]]")
##    #Clampingmask
##arcpy.DefineProjection_management(ClampingMask_Auto2ndHSMap, "PROJCS['CH1903_LV03',GEOGCS['GCS_CH1903',DATUM['D_CH1903',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Hotine_Oblique_Mercator_Azimuth_Center'],PARAMETER['False_Easting',600000.0],PARAMETER['False_Northing',200000.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Azimuth',90.0],PARAMETER['Longitude_Of_Center',7.439583333333333],PARAMETER['Latitude_Of_Center',46.95240555555556],UNIT['Meter',1.0]]")




#Move this to the end of the whole set of processes (that require spatial analyst extension)
arcpy.CheckInExtension('Spatial')


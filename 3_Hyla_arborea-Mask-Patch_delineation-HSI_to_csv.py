#-------------------------------------------------------------------------------
# Script name:       Hyla_arborea-Mask-patch_delineation-HSI_to_csv
# Purpose:           Mask application + Patch (node) delineation + Calculation of Mean HSI per habitat patch
# Purpose:           For Hyla arborea
# From article:      "Sensitivity of habitat network models to changes in maximum dispersal distance"
# Originally developed for article "Predicting species occurrences with habitat network models"
# Author:      Damian O. Ortiz-Rodriguez, Antoine Guisan, Rolf Holderegger, Maarten J. van Strien
# 1st version created:     29/11/2017

#-------------------------------------------------------------------------------

def main():
    pass

if __name__ == '__main__':
    main()

# importing functions
import arcpy, time
from arcpy import env
from arcpy.sa import *
import igraph as ig
import numpy as np

#Import functions
import arcgisscripting
import os
import enum
import dbf

# Setting the initial environment

arcpy.SetProduct('ArcInfo')
arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True

# ######### Applying Mask (& Defining projection of HSM outputs ################################################################
Mask = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Mask\\Mask_Amph_habs.tif"
ContHSMap = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\Auto3rd_Hyarb_ensemble_1.tif"
binaryHSMap = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\ROCBin_Auto3rd_Hyarb_ensemble_1.tif"
Masked_binaryHSMap = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\Masked_Hy_arb_NoPseudorep_ROCbin.tif"

#Define projection of HSM maps
    #Continuous
arcpy.DefineProjection_management(ContHSMap, "PROJCS['CH1903_LV03',GEOGCS['GCS_CH1903',DATUM['D_CH1903',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Hotine_Oblique_Mercator_Azimuth_Center'],PARAMETER['False_Easting',600000.0],PARAMETER['False_Northing',200000.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Azimuth',90.0],PARAMETER['Longitude_Of_Center',7.439583333333333],PARAMETER['Latitude_Of_Center',46.95240555555556],UNIT['Meter',1.0]]")
    #ROCBin
arcpy.DefineProjection_management(binaryHSMap, "PROJCS['CH1903_LV03',GEOGCS['GCS_CH1903',DATUM['D_CH1903',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Hotine_Oblique_Mercator_Azimuth_Center'],PARAMETER['False_Easting',600000.0],PARAMETER['False_Northing',200000.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Azimuth',90.0],PARAMETER['Longitude_Of_Center',7.439583333333333],PARAMETER['Latitude_Of_Center',46.95240555555556],UNIT['Meter',1.0]]")

#Reclassify mask to have only '1' values
Mask_Amph_habs_only1 = Reclassify(Mask, "Value",RemapRange([[0,"NODATA"],[1,1]]))

# Execute ExtractByMask
Masked_binaryHSMap = ExtractByMask(binaryHSMap, Mask_Amph_habs_only1)
# Save the output
Masked_binaryHSMap.save("E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\Masked_Hy_arb_NoPseudorep_ROCbin.tif")


# ################################ Get coded patches ################################################
# From the masked ROCBinarized HSM output, this performs the RegionGroup Algorithm

# Input settings
wrkspc = 'E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\'
env.workspace = wrkspc
##Masked_binaryHSMap = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\Masked_Hy_arb_NoPseudorep_ROCbin.tif" #Define here if you want to run by parts

#Change the binary habitat suitability (HS) file to an integer as this will save huge amounts of space.
binaryHSMapInt = Int(Masked_binaryHSMap)

#Reclassify the HS map to have only 1 values for the suitable patches. All the rest becomes NoData
OnlysuitableHSMap = Reclassify(binaryHSMapInt, "Value",RemapRange([[0,"NODATA"],[1,1]]))

#Code patches of continuous habitat with a unique number
habPatchCode = RegionGroup(OnlysuitableHSMap, "EIGHT", "WITHIN", "NO_LINK", "")

#Save the patch code raster
habPatchCode.save("habPatchCode_NoPseudorep.tif")



# ################################ Rasterize all spp records ################################################
#### This only needs to be done once, so this process is only present in this script, and not in the scripts of other species ### 

# ###Do Shp's for ind. spp. by year #######################
#Set variables
inDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Ind_spp\\By_year\\"
outDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Ind_spp\\By_year\\Raster\\" #<-- CHANGE

#Set the working directory
os.chdir(inDir)

# do the same for all the shapefiles in the directory - https://community.esri.com/thread/50708
#Set folder from which files are selected and processed
for filename in os.listdir(inDir):
    if filename.endswith(".shp"):
        fName = arcpy.Describe(filename).basename
        outTIF = outDir + "\\" + fName + u".tif"
        arcpy.PointToRaster_conversion(filename, "ART", outTIF, "MOST_FREQUENT", "NONE", "100")

del filename, outTIF
#END

# ################################ Get raster of patches with presence/absences of all species ###########################
#Get zonal statistics as table to get count rasters of amphibian observed in a certain year in the different habitat patches, binarize in R
#Set variables
##habPatchCode = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif"
inDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Spp\\Ind_spp\\By_year\\Raster\\"
outDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Presence_absence\\Zonal_stats_tables\\"


#Set the working directory
os.chdir(inDir)

for filename in os.listdir(inDir):
    if filename.endswith(".tif"):
        fName = arcpy.Describe(filename).basename
        outDBF = outDir + "\\" + fName + u"_zstat.dbf"

        ZonalStatisticsAsTable(habPatchCode, "VALUE", filename, outDBF, "DATA", "ALL")

del filename, outDBF
#END

#  ################################ Append to attribute table of each year ################################

inDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Presence_absence\\Zonal_stats_tables\\"
outDir = u"E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Presence_absence\\Zonal_stats_tables\\"
#Do this with the files involved not open in ArcGIS, specially HabPatchCode

#Set the working directory
os.chdir(inDir)

for filename in os.listdir(inDir):
    if filename.endswith(".dbf"):
        fName = arcpy.Describe(filename).basename
        print(fName)
        fieldName = fName.replace("_zstat", "_c").replace("Hyla_arborea_20", "Hy_arb").replace("Alytes_obstetricans_20", "Al_obs").replace("Bombina_variegata_20", "Bo_var").replace("Bufo_bufo_20", "Bu_buf").replace("Epidalea__calamita_20", "Ep_cal").replace("Ichthyosaura_alpestris_20", "Ic_alp").replace("Lissotriton_helveticus_20", "Li_hel").replace("Pelophylax_aggr_20", "Pe_agg").replace("Pelophylax_ridibundus_20", "Pe_rid").replace("Rana_dalmatina_20", "Ra_dal").replace("Rana_temporaria_20", "Ra_tem").replace("Triturus_carnifex_20", "Tr_car").replace("Triturus_cristatus_20", "Tr_cri")
        print(fieldName)


        arcpy.JoinField_management (habPatchCode, "VALUE", filename, "VALUE", ["SUM"])
        #arcpy.AlterField_management ("E:\\PhD\\Automated_Run_BioCHECNET\\Network_setup\\habPatchCode2.tif", ["SUM"], [fieldName])

        arcpy.AddField_management (habPatchCode, fieldName, "LONG", 10)
        arcpy.CalculateField_management (habPatchCode, fieldName, "!SUM!", "PYTHON")
        arcpy.DeleteField_management (habPatchCode, "SUM")

del filename, fieldName,
#END


# #################### Export attribute table as .csv (for input into R)
# Local variables:
Presence_absence_folder = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Presence_absence\\"
Macrotable_csv = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Presence_absence\\Macrotable_NoPseudorep.csv"
# Process: Table to Table (to csv)
arcpy.TableToTable_conversion(habPatchCode, Presence_absence_folder, "Macrotable_NoPseudorep.csv")




# ################### HSI to CSV #############################################
# Input settings
wrkspc = 'E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup'
env.workspace = wrkspc

# #Turn HSI raster into integer, so as to get attribute table
ContHSMap = "E:\\PhD\\Automated_Run_BioCHECNET\\\\NoPseudorep\\HSM\\Auto3rd_Hyarb_ensemble_1.tif"
filename = Int(ContHSMap)
filename.save("E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\Auto3rd_Hyarb_NoPseudorep_int.tif")

# #Calculate HSI zonal stats
#vars
#habPatchCode = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif"
#filename = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\Auto3rd_Hyarb_NoPseudorep_int.tif"
HSI_DBF = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\HSI_habPatchCode_NoPseudorep.dbf"


# # Zonal stats as table
ZonalStatisticsAsTable(habPatchCode, "VALUE", filename, HSI_DBF, "DATA", "ALL")

#Write csv to import into R
arcpy.TableToTable_conversion(HSI_DBF, wrkspc, "HSI_habPatchCode_NoPseudorep.csv")
#Attention: HSI relevant Field is the "MEAN"


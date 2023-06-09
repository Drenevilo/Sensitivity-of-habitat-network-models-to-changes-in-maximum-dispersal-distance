#-------------------------------------------------------------------------------
# Name:        Pelophylax_agg-Mask-patch_delineation-HSI_to_csv
# Purpose:     Mask application + Patches_then_Networks + Processing of HSI
# Purpose:     For Pelophylax agg.

# From article: "Sensitivity of habitat network models to changes in maximum dispersal distance"

# Author:      Damian O. Ortiz-Rodriguez, Antoine Guisan, Rolf Holderegger, Maarten J. van Strien
# 1st version created:  15-04-2019
# Copyright:   (c) damiano 2019
# Licence:     <your licence>
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
#arcpy.env.parallelProcessingFactor = "75%" # Use 75% of the cores on the machine.

# ######### NOTE: Definition of projection of HSM outputs (continuous and binarized) was done in the ArcMap GUI for this species.
# #########        For the scripted process in Python, see the equivalent script for Hyla arborea

# ################################ Apply mask #######################################################
wrkspc = 'C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\'
env.workspace = wrkspc

# Local variables
binaryHSMap = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\ROCBin_Auto1st_Peagg_ensemble_1.tif"
Mask_Amph_habs = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Mask_Amph_habs.tif"
Masked_binaryHSMap = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Masked_Pe_agg_ensemble_1_ROCbin.tif"

#Reclassify mask to have only '1' values
Mask_Amph_habs_only1 = Reclassify(Mask_Amph_habs, "Value",RemapRange([[0,"NODATA"],[1,1]]))

# Execute ExtractByMask
#Change coord. system 1st #Transformations
tempEnvironment0 = arcpy.env.outputCoordinateSystem
arcpy.env.outputCoordinateSystem = "PROJCS['CH1903_LV03',GEOGCS['GCS_CH1903',DATUM['D_CH1903',SPHEROID['Bessel_1841',6377397.155,299.1528128]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],PROJECTION['Hotine_Oblique_Mercator_Azimuth_Center'],PARAMETER['False_Easting',600000.0],PARAMETER['False_Northing',200000.0],PARAMETER['Scale_Factor',1.0],PARAMETER['Azimuth',90.0],PARAMETER['Longitude_Of_Center',7.439583333333333],PARAMETER['Latitude_Of_Center',46.95240555555556],UNIT['Meter',1.0]]"
tempEnvironment1 = arcpy.env.geographicTransformations
arcpy.env.geographicTransformations = ""
Masked_binaryHSMap = ExtractByMask(binaryHSMap, Mask_Amph_habs_only1)
# Save the output
Masked_binaryHSMap.save("C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Masked_Pe_agg_ensemble_1_ROCbin.tif")



# ################################ Get coded patches ################################################
# From the masked ROCBinarized HSM output, this performs the RegionGroup Algorithm

# Input settings
wrkspc = 'C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\'
env.workspace = wrkspc

Masked_binaryHSMap = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Masked_Pe_agg_ensemble_1_ROCbin.tif"


#Change the binary habitat suitability (HS) file to an integer as this will save huge amounts of space.
binaryHSMapInt = Int(Masked_binaryHSMap)

#Reclassify the HS map to have only 1 values for the suitable patches. All the rest becomes NoData
OnlysuitableHSMap = Reclassify(binaryHSMapInt, "Value",RemapRange([[0,"NODATA"],[1,1]]))

#Code patches of continuous habitat with a unique number
habPatchCode = RegionGroup(OnlysuitableHSMap, "EIGHT", "WITHIN", "NO_LINK", "")

#Save the patch code raster
habPatchCode.save("habPatchCode_Peagg1.tif")



# ################################ Rasterize all spp records ################################################
#### This process only needs to be done once, so see it in the equivalent script for Hyla arborea ####


# ################################ Get raster of patches with presence/absences of all species ###########################
#Get zonal statistics as table to get count rasters of amphibian observed in a certain year in the different habitat patches, binarize in R

#Set variables
##habPatchCode = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\habPatchCode_Peagg1.tif"
inDir = u"C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\GeoData\\Spp\\Ind_spp\\By_year\\Raster\\"
outDir = u"C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Presence_absence\\Pe_agg\\Zonal_stats_tables\\"


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

inDir = u"C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Presence_absence\\Pe_agg\\Zonal_stats_tables\\"
outDir = u"C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Presence_absence\\Pe_agg\\Zonal_stats_tables\\"
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
Presence_absence_folder = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Presence_absence\\"
Macrotable_csv = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Presence_absence\\Macrotable_Peagg.csv"

# Process: Table to Table (to csv)
arcpy.TableToTable_conversion(habPatchCode, Presence_absence_folder, "Macrotable_Peagg.csv")



# ################### HSI to CSV #############################################
# Input settings
wrkspc = 'C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup'
env.workspace = wrkspc

# #Turn HSI raster into integer, so as to get attribute table
HSI_double = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Auto1st_Peagg_ensemble_1.tif"

filename = Int(HSI_double)
filename.save("C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Auto1st_Peagg_ensemble_int.tif")


# #Calculate HSI zonal stats
#vars
habPatchCode = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\habPatchCode_Peagg1.tif"
filename = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Auto1st_Peagg_ensemble_int.tif"
outDBF = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\HSI_habPatchCode_Peagg1.dbf"

# # Zonal stats as table
ZonalStatisticsAsTable(habPatchCode, "VALUE", filename, outDBF, "DATA", "ALL")

#Write csv to import into R
arcpy.TableToTable_conversion(outDBF, wrkspc, "HSI_habPatchCode_Peagg1.csv")
#Attention: HSI relevant Field is the "MEAN"
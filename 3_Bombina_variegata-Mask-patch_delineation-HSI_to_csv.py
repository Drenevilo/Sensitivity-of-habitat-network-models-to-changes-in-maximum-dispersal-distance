#-------------------------------------------------------------------------------
# Name:        Bombina_variegata-Mask-patch_delineation-HSI_to_csv
# Purpose:     Mask application + Patches_then_Networks + Processing of HSI
# Purpose:     For Bombina_variegata

# From article: "Sensitivity of habitat network models to changes in maximum dispersal distance"

# Author:      Damian O. Ortiz-Rodriguez, Antoine Guisan, Rolf Holderegger, Maarten J. van Strien
# 1st version created:  29-03-2019
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

# Local variables:
binaryHSMap = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\ROCBin_Auto1st_Bovar_ensemble_1.tif"
Mask_Amph_habs = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Mask_Amph_habs.tif"
Masked_binaryHSMap = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Masked_Bo_var_ensemble_1_ROCbin.tif"


#Reclassify mask to have only '1' values
Mask_Amph_habs_only1 = Reclassify(Mask_Amph_habs, "Value",RemapRange([[0,"NODATA"],[1,1]]))

# Execute ExtractByMask
Masked_binaryHSMap = ExtractByMask(binaryHSMap, Mask_Amph_habs_only1)
# Save the output
Masked_binaryHSMap.save("C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Masked_Bo_var_ensemble_1_ROCbin.tif")



# ################################ Get coded patches ################################################
# From the masked ROCBinarized HSM output, this performs the RegionGroup Algorithm

# Input settings
wrkspc = 'C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\'
env.workspace = wrkspc

Masked_binaryHSMap = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Masked_Bo_var_ensemble_1_ROCbin.tif"


#Change the binary habitat suitability (HS) file to an integer as this will save huge amounts of space.
binaryHSMapInt = Int(Masked_binaryHSMap)

#Reclassify the HS map to have only 1 values for the suitable patches. All the rest becomes NoData
OnlysuitableHSMap = Reclassify(binaryHSMapInt, "Value",RemapRange([[0,"NODATA"],[1,1]]))

#Code patches of continuous habitat with a unique number
habPatchCode = RegionGroup(OnlysuitableHSMap, "EIGHT", "WITHIN", "NO_LINK", "")

#Save the patch code raster and save it as a polygon and point shapefile (centroids)

habPatchCode.save("habPatchCode_Bovar1.tif")



# ################################ Rasterize all spp records ################################################
#### This process only needs to be done once, so see it in the equivalent script for Hyla arborea ####


# ################################ Get raster of patches with presence/absences of all species ###########################
#Get zonal statistics as table to get count rasters of amphibian observed in a certain year in the different habitat patches, binarize in R

#Set variables
##habPatchCode = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\habPatchCode_Bovar1.tif"

inDir = u"C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\GeoData\\Spp\\Ind_spp\\By_year\\Raster\\"
outDir = u"C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Presence_absence\\Bo_var\\Zonal_stats_tables\\"


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

inDir = u"C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Presence_absence\\Bo_var\\Zonal_stats_tables\\"
outDir = u"C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Presence_absence\\Bo_var\\Zonal_stats_tables\\"
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
arcpy.TableToTable_conversion(habPatchCode, Presence_absence_folder, "Macrotable_NoPseudorep.csv", "", "VALUE \"VALUE\" false true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,VALUE,-1,-1;COUNT \"COUNT\" false true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,COUNT,-1,-1;Pe_agg06_c \"Pe_agg06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg06_c,-1,-1;Bu_buf12_c \"Bu_buf12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf12_c,-1,-1;Pe_rid13_c \"Pe_rid13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid13_c,-1,-1;Tr_car13_c \"Tr_car13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_car13_c,-1,-1;Pe_agg07_c \"Pe_agg07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg07_c,-1,-1;Hy_arb12_c \"Hy_arb12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb12_c,-1,-1;Tr_cri09_c \"Tr_cri09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri09_c,-1,-1;Bu_buf13_c \"Bu_buf13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf13_c,-1,-1;Pe_rid06_c \"Pe_rid06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid06_c,-1,-1;Ra_tem06_c \"Ra_tem06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem06_c,-1,-1;Pe_agg08_c \"Pe_agg08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg08_c,-1,-1;Hy_arb13_c \"Hy_arb13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb13_c,-1,-1;Tr_cri10_c \"Tr_cri10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri10_c,-1,-1;Pe_agg09_c \"Pe_agg09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg09_c,-1,-1;Ra_tem07_c \"Ra_tem07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem07_c,-1,-1;Pe_rid14_c \"Pe_rid14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid14_c,-1,-1;Bu_buf14_c \"Bu_buf14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf14_c,-1,-1;Pe_rid07_c \"Pe_rid07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid07_c,-1,-1;Hy_arb14_c \"Hy_arb14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb14_c,-1,-1;Tr_cri11_c \"Tr_cri11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri11_c,-1,-1;Ra_tem08_c \"Ra_tem08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem08_c,-1,-1;Pe_agg10_c \"Pe_agg10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg10_c,-1,-1;Bu_buf15_c \"Bu_buf15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf15_c,-1,-1;Tr_cri12_c \"Tr_cri12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri12_c,-1,-1;Ra_tem09_c \"Ra_tem09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem09_c,-1,-1;Pe_agg11_c \"Pe_agg11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg11_c,-1,-1;Hy_arb15_c \"Hy_arb15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb15_c,-1,-1;Pe_rid15_c \"Pe_rid15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid15_c,-1,-1;Pe_rid08_c \"Pe_rid08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid08_c,-1,-1;Ep_cal06_c \"Ep_cal06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal06_c,-1,-1;Tr_cri13_c \"Tr_cri13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri13_c,-1,-1;Ra_tem10_c \"Ra_tem10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem10_c,-1,-1;Bu_buf06_c \"Bu_buf06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf06_c,-1,-1;Ep_cal07_c \"Ep_cal07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal07_c,-1,-1;Pe_agg12_c \"Pe_agg12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg12_c,-1,-1;Ra_tem11_c \"Ra_tem11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem11_c,-1,-1;Tr_cri14_c \"Tr_cri14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri14_c,-1,-1;Ep_cal08_c \"Ep_cal08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal08_c,-1,-1;Pe_agg13_c \"Pe_agg13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg13_c,-1,-1;Tr_cri15_c \"Tr_cri15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri15_c,-1,-1;Ra_tem12_c \"Ra_tem12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem12_c,-1,-1;Pe_rid09_c \"Pe_rid09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid09_c,-1,-1;Ep_cal09_c \"Ep_cal09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal09_c,-1,-1;Bu_buf07_c \"Bu_buf07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf07_c,-1,-1;Ra_tem13_c \"Ra_tem13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem13_c,-1,-1;Pe_agg14_c \"Pe_agg14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg14_c,-1,-1;Ep_cal10_c \"Ep_cal10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal10_c,-1,-1;Hy_arb06_c \"Hy_arb06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb06_c,-1,-1;Pe_rid10_c \"Pe_rid10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid10_c,-1,-1;Ra_tem14_c \"Ra_tem14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem14_c,-1,-1;Bu_buf08_c \"Bu_buf08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf08_c,-1,-1;Ep_cal11_c \"Ep_cal11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal11_c,-1,-1;Pe_agg15_c \"Pe_agg15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_agg15_c,-1,-1;Ra_tem15_c \"Ra_tem15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_tem15_c,-1,-1;Hy_arb07_c \"Hy_arb07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb07_c,-1,-1;Ep_cal12_c \"Ep_cal12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal12_c,-1,-1;Pe_rid11_c \"Pe_rid11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid11_c,-1,-1;Bu_buf09_c \"Bu_buf09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf09_c,-1,-1;Ep_cal13_c \"Ep_cal13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal13_c,-1,-1;Tr_cri06_c \"Tr_cri06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri06_c,-1,-1;Tr_car06_c \"Tr_car06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_car06_c,-1,-1;Hy_arb08_c \"Hy_arb08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb08_c,-1,-1;Tr_cri07_c \"Tr_cri07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri07_c,-1,-1;Ep_cal14_c \"Ep_cal14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal14_c,-1,-1;Bu_buf10_c \"Bu_buf10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf10_c,-1,-1;Li_hel06_c \"Li_hel06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel06_c,-1,-1;Pe_rid12_c \"Pe_rid12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Pe_rid12_c,-1,-1;Ep_cal15_c \"Ep_cal15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ep_cal15_c,-1,-1;Hy_arb09_c \"Hy_arb09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb09_c,-1,-1;Tr_cri08_c \"Tr_cri08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_cri08_c,-1,-1;Tr_car07_c \"Tr_car07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_car07_c,-1,-1;Tr_car10_c \"Tr_car10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_car10_c,-1,-1;Ic_alp07_c \"Ic_alp07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp07_c,-1,-1;Li_hel07_c \"Li_hel07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel07_c,-1,-1;Hy_arb10_c \"Hy_arb10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb10_c,-1,-1;Tr_car11_c \"Tr_car11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_car11_c,-1,-1;Ra_dal13_c \"Ra_dal13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal13_c,-1,-1;Bu_buf11_c \"Bu_buf11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bu_buf11_c,-1,-1;Tr_car08_c \"Tr_car08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_car08_c,-1,-1;Hy_arb11_c \"Hy_arb11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Hy_arb11_c,-1,-1;Tr_car12_c \"Tr_car12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_car12_c,-1,-1;Al_obs11_c \"Al_obs11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs11_c,-1,-1;Ic_alp08_c \"Ic_alp08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp08_c,-1,-1;Al_obs12_c \"Al_obs12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs12_c,-1,-1;Tr_car09_c \"Tr_car09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Tr_car09_c,-1,-1;Al_obs13_c \"Al_obs13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs13_c,-1,-1;Ra_dal14_c \"Ra_dal14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal14_c,-1,-1;Al_obs14_c \"Al_obs14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs14_c,-1,-1;Ic_alp15_c \"Ic_alp15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp15_c,-1,-1;Al_obs15_c \"Al_obs15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs15_c,-1,-1;Li_hel08_c \"Li_hel08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel08_c,-1,-1;Ic_alp09_c \"Ic_alp09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp09_c,-1,-1;Li_hel09_c \"Li_hel09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel09_c,-1,-1;Ra_dal15_c \"Ra_dal15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal15_c,-1,-1;Li_hel10_c \"Li_hel10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel10_c,-1,-1;Ic_alp10_c \"Ic_alp10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp10_c,-1,-1;Ra_dal06_c \"Ra_dal06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal06_c,-1,-1;Bo_var06_c \"Bo_var06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var06_c,-1,-1;Li_hel11_c \"Li_hel11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel11_c,-1,-1;Ic_alp11_c \"Ic_alp11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp11_c,-1,-1;Ra_dal07_c \"Ra_dal07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal07_c,-1,-1;Bo_var07_c \"Bo_var07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var07_c,-1,-1;Li_hel12_c \"Li_hel12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel12_c,-1,-1;Ic_alp12_c \"Ic_alp12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp12_c,-1,-1;Ra_dal08_c \"Ra_dal08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal08_c,-1,-1;Bo_var08_c \"Bo_var08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var08_c,-1,-1;Li_hel13_c \"Li_hel13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel13_c,-1,-1;Ic_alp13_c \"Ic_alp13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp13_c,-1,-1;Ra_dal09_c \"Ra_dal09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal09_c,-1,-1;Bo_var09_c \"Bo_var09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var09_c,-1,-1;Li_hel14_c \"Li_hel14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel14_c,-1,-1;Ic_alp14_c \"Ic_alp14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp14_c,-1,-1;Li_hel15_c \"Li_hel15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Li_hel15_c,-1,-1;Ra_dal10_c \"Ra_dal10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal10_c,-1,-1;Bo_var10_c \"Bo_var10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var10_c,-1,-1;Ra_dal11_c \"Ra_dal11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal11_c,-1,-1;Bo_var11_c \"Bo_var11_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var11_c,-1,-1;Ic_alp06_c \"Ic_alp06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ic_alp06_c,-1,-1;Ra_dal12_c \"Ra_dal12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Ra_dal12_c,-1,-1;Bo_var12_c \"Bo_var12_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var12_c,-1,-1;Bo_var13_c \"Bo_var13_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var13_c,-1,-1;Bo_var14_c \"Bo_var14_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var14_c,-1,-1;Bo_var15_c \"Bo_var15_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Bo_var15_c,-1,-1;Al_obs06_c \"Al_obs06_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs06_c,-1,-1;Al_obs07_c \"Al_obs07_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs07_c,-1,-1;Al_obs08_c \"Al_obs08_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs08_c,-1,-1;Al_obs09_c \"Al_obs09_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs09_c,-1,-1;Al_obs10_c \"Al_obs10_c\" true true false 10 Long 0 10 ,First,#,E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\habPatchCode_NoPseudorep.tif\\Band_1,Al_obs10_c,-1,-1", "")



# ################### HSI to CSV #############################################
# Input settings
wrkspc = 'C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup'
#'E:\\PhD\\Automated_Run_BioCHECNET\\Network_setup'
env.workspace = wrkspc

# #Turn HSI raster into integer, so as to get attribute table
HSI_double = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Auto1st_Bovar_ensemble_1.tif"

filename = Int(HSI_double)
filename.save("C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Auto1st_Bovar_ensemble_int.tif")


# #Calculate HSI zonal stats
#vars
habPatchCode = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\habPatchCode_Bovar1.tif"
filename = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\HSM\\Auto1st_Bovar_ensemble_int.tif"
outDBF = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\HSI_habPatchCode_Bovar1.dbf"


# # Zonal stats as table
ZonalStatisticsAsTable(habPatchCode, "VALUE", filename, outDBF, "DATA", "ALL")

#Write csv to import into R
arcpy.TableToTable_conversion(outDBF, wrkspc, "HSI_habPatchCode_Bovar1.csv")
#Attention: HSI relevant Field is the "MEAN"

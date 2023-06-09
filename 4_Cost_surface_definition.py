#-------------------------------------------------------------------------------
# Script name:      Cost_surface_definition

# Purpose:      Cost surface definition

# Author:      Damian O. Ortiz-Rodriguez, Antoine Guisan, Rolf Holderegger, Maarten J. van Strien
# Created:     07/12/2017

# Used for article "Sensitivity of habitat network models to changes in maximum dispersal distance"
# Originally developed for article "Predicting species occurrences with habitat network models"
#-------------------------------------------------------------------------------

#### NOTE: Commented-out lines (#) and elements with ###*** were not used in the article "Sensitivity of habitat network models to changes in maximum dispersal distance" ####


def main():
    pass

if __name__ == '__main__':
    main()


# Input settings
binaryHSMap = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\ROCBin_Auto3rd_Hyarb_ensemble_1.tif"
#unmasked binary map is used for snap raster and extent
trafShape = "E:\\PhD\\Automated_Run_BioCHECNET\\GeoData\\Processes_w_traffic\\Traffic_CH_bothdirs_notunnels.shp" ###***
hsRast = "E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\HSM\\Auto3rd_Hyarb_ensemble_1.tif" #continuous hab.suit. map (###***)
wrkspc = 'E:\\PhD\\Automated_Run_BioCHECNET\\NoPseudorep\\Network_setup\\'


# Values below are for the tree frog Hyla arborea: (###***)
d0 = 200 #Median dispersal distance
wc = 0.55 #Width of car (m)
wa = 0.02 #Width of animal (m)
lc = 5. #length of car (m)
la = 0.1 #length of animal (m)
vc = 22.2 # = 80 km/h. Speed of car(m/s) Can only cross outside of built area
va = 0.033 #Speed of animal (m/s)

# importing functions
import arcpy, time, os
from arcpy import env
from arcpy.sa import *
import numpy as np

# Set the python working directory
os.chdir(wrkspc)

# Setting the initial environment
env.workspace = wrkspc
arcpy.SetProduct('ArcInfo')
arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True

#### Extent, cell size, and coordinate system are the same for all the species, so it does not matter which species' suitability map is taken as reference to generate a uniform cost surface ####

# Set the extent of all other raster to this buffered shapefile.
arcpy.env.extent = arcpy.Describe(binaryHSMap).extent
#Snap all other rasters to this habitat map
arcpy.env.snapRaster = Raster(binaryHSMap)
# Get the cell size of the raster
cellSize = int(arcpy.GetRasterProperties_management(binaryHSMap, "CELLSIZEX").getOutput(0))
# Get the coordinate system
coordSys = arcpy.Describe(binaryHSMap).spatialReference

# CREATE UNIFORM RESISTANCE SURFACE
costSurface = CreateConstantRaster(1, "INTEGER", "", "")
arcpy.DefineProjection_management(costSurface, coordSys)
costSurface.save("UniformCS.tif")

#### From here onwards ###***

## CREATE TRAFFIC RESISTANCE SURFACE
## The default cost value of the raster is 1.
## Calculate the capacity of each edge in the forest network following formula's from: van Langevelde F, Jaarsma CF (2009).
## The probability that an animal successfully crosses the road Pa = exp(-volij*(((wc+la)/va)+((lc+wa)/vc))). The variables are described below.
## The probability of dispersal can be translated to a cost value with this equation: cost = ln(prob)/log(.5)*d0
## We divide daily traffic by (24*60*60) to convert to cars per seconds.
#arcpy.PolylineToRaster_conversion(trafShape, "TOT_TRAFF", "traffic.tif", "MAXIMUM_LENGTH", "", 10) #Read lineshapefile of NPVM data into python.
#costSurface = Raster("traffic.tif")
#costSurface = Float(costSurface)
#costSurface = Aggregate(costSurface, cellSize/10, "MAXIMUM") #From the fine resolution traffic.tif raster we make a course resolution raster.
#costSurface = costSurface/(24*60*60) #Cars per day are transformed to vehicles per second.
#costSurface.save("TrafficPSec.tif")
#costSurface = Exp(-costSurface*(((wc+la)/va)+((lc+wa)/vc))) #Transform traffic volumes to probabilities.
#costSurface = (Ln(costSurface)/np.log(.5)*d0)/cellSize #Probabilities are transformed to costs. Since the costs in least-cost paths will be multiplied by the raster cell size, the costs are devided by the cell size.
#costSurface = costSurface+1 # 1 is added to the cost values to represent the "background cost value"
#costSurface = Con(IsNull(costSurface),1,costSurface)
#costSurface.save("TrafficCS.tif")

#del costSurface


## CREATE HABITAT SUITABILITY RESISTANCE SURFACE
## The Habitat Suitability raster ranges between 0 (unsuitable) and 1000 (suitable).
## We simply assume that the probability of dispersing through the most suitable terrain is 1 and throught the most unsuitable terrain is 0
## Therefore, the HS raster is divided by the maximum suitability (this way the most suitable terrain becomes 1).
## Then we transform probabilities in the same way as the raster above.
#hsSurface = Raster(hsRast)
#maxVal = int(arcpy.GetRasterProperties_management(hsSurface, "MAXIMUM").getOutput(0))
#hsSurface = Float(hsSurface)
#hsSurface = hsSurface/maxVal
#hsSurface = (Ln(hsSurface)/np.log(.5)*d0)/cellSize
#hsSurface = hsSurface+1 # 1 is added to the cost values to represent the "background cost value"
#hsSurface.save("habsuitCS.tif")

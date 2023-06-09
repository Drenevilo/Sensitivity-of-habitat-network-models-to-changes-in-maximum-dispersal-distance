#-------------------------------------------------------------------------------
# Name:        Create_networks_plus_HabAv_Bovar_d0_alternatives
# Purpose:     Implement different maximum dispersal distances
#              Create different networks changing the d0 parameter
#              Edge definition (Generation of networks) +
#              calculation of Habitat availability predictor
# Purpose:     for Bombina_variegata
# Author:      Damian O. Ortiz-Rodriguez, Antoine Guisan, Maarten J. van Strien
#
# 1st version created:     04-11-2021
# Copyright:   (c) damiano 2021
# Licence:     <your licence>
#-------------------------------------------------------------------------------

def main():
    pass

if __name__ == '__main__':
    main()


# ##########################################
# ### Define THE FUNCTIONS ####################
# ##########################################

# Function to check the coordinate system
def check_ref(feature,coord_sys):
    if arcpy.Describe(feature).spatialReference.name != coord_sys:
        print "WARNING!"+str(feature)+" is in a different coordinate system as "+coord_sys

# Function to save the output of the functions "make_delaunay_edges" and "make_gabriel_edges".
def save_edge_array_as_shapefile(edges, shapefile, attrib_dat = -1):
    """
    Function to save the output of the functions "make_delaunay_edges" and "make_gabriel_edges".
    Requires:
        edges: output array from "make_delaunay_edges" and "make_gabriel_edges"
        shapefile: name of the shapefile (ending with ".shp")
        attrib_dat: (optional) array of equal length as the number of edges with values that have to be added to the attribute table of the shapefile.
    From: http://resources.arcgis.com/en/help/main/10.1/index.html#/Array/018z0000006n000000/
    """
    import arcpy
    arcpy.env.overwriteOutput = True
    n_row = edges.shape[0]
    edge_coords = edges.reshape([n_row,2,2]).tolist()
    # A list of features and coordinate pairs
    feature_info = edge_coords
    # A list that will hold each of the Polyline objects
    features = []
    for feature in feature_info:
        features.append(
            arcpy.Polyline(
                arcpy.Array([arcpy.Point(*coords) for coords in feature])))
    try:
        # Persist a copy of the Polyline objects using CopyFeatures
        arcpy.CopyFeatures_management(features, shapefile)
        # Repair the geometry of the shapefile
        arcpy.RepairGeometry_management(shapefile)
        arcpy.RepairGeometry_management(shapefile)
    except:
        print arcpy.GetMessages()
    # Loop to add the values in attrib_dat to the column "attrib" in the attribute table of the shapefile.
    if attrib_dat != -1:
        arcpy.AddField_management(shapefile,"attrib","DOUBLE")
        cur = arcpy.UpdateCursor(shapefile)
        counter = 0
        for row in cur:
            row.setValue('attrib', attrib_dat[counter])
            counter += 1
            cur.updateRow(row)
        del cur, row


# Function to change the array returned from the function "spatial_join_edge_attributes" to which the capacity of the edges has been added into a igraph graph.
def edge_table_to_igraph(patch_table, edge_table, xy_names, patch_dict = -1, capacity_col = -1, weight_col = -1, probability_col = -1):
    """
    Function to change the array returned from the function "spatial_join_edge_attributes" to which the capacity of the edges has been added into a igraph graph.
    Requires:
        patch_table: python array with the names of all of the nodes
        edge_table: numpy array with a row for each link
        xy_names: two columns in "edge_table" containing the names of the nodes.
        patch_dict: (optional) A dictionary with values that should be assigned to the nodes (NEEDS TO BE UPDATES)
        capacity_col: (optional) Column in "edge_table" containing the capacity of each edge.
        weight_col: (optional) Column in the "edge_table" containing the weights of each edge.
    """
    import igraph as ig
    names_x1y1 = edge_table[:,xy_names[0]].astype(int).astype(str)
    names_x2y2 = edge_table[:,xy_names[1]].astype(int).astype(str)
    G = ig.Graph(directed = False)
    G.add_vertices(np.array(patch_table).astype(str))
    G.add_edges(np.column_stack((names_x1y1,names_x2y2)).tolist())
    if patch_dict != -1:
        G.vs['patch'] = [patch_dict.get(int(enn), [])[0] for enn in G.vs['name']] # The number of the patch to with a node belongs
        G.vs['edge_bool'] = [patch_dict.get(int(enn), [])[1] for enn in G.vs['name']] # Boolean indicator whether a batch is on the edge of a patch or on the edge of a gap
        G.vs['patch_edge'] = np.multiply(G.vs['patch'], G.vs['edge_bool']) #Only the patches on an edge are assigned their patch numbers. All other patches have 0.
    if weight_col != -1:
        G.es["weight"] = edge_table[:,weight_col].tolist()
    if capacity_col != -1:
        G.es["capacity"] = edge_table[:,capacity_col].tolist()
    if probability_col != -1:
        G.es["probability"] = edge_table[:,probability_col].tolist()
    return G

# ##########################################


# Input settings
wrkspc = "C:\\Users\\damiano\\Documents\\PhD\\Sensitivity_DispDist\\Networks\\" #Where the network files are created
inputFiles = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\" #Where the cost surfaces are sourced from
costSurface = inputFiles+"UniformCS.tif" #UniformCostSurface
run_name = "Bovar_10KmDispDist" #The name that will be given to the output files
# #### For each run of the different dispersal distance just change the distance to name the run
# #### i.e. 300m, 1Km, 2Km, 4Km, 6Km, 8Km, 10Km

# With the p2p function in the R-package PopGenReport cost distances are converted to probabilities:
#           The function is: prob = exp(((x/d0) * log(p)))
#           x = cost distance
#           d0 = dispersal distance of proportion p individuals
#           p = proportion of individuals dispersing over distance d0. For example, d0=100, p = 0.5 ->
#           50% procent of all migrating individuals go up to 100 m.
# We fix p at 0.5, so then d0 is equal to the median dispersal distance. Costs are transformed to probabilities with
# the above function
d0 = 753 # Median dispersal distance, change for each species
# In order to speed up the calculations a maximum cost-distance can be indicated in the arcpy cost-distance function.
#  This maximum cost-distance is the distance beyond which the probability of dispersal drops below
# a certain probability min_prob and is calculated as follows: max_cost = log(min_prob)/log(p)*d0 with p = 0.5 (see above)

#Different dispersal distance settings applied - d0 chosen was closest possible to each maximum dispersal distance in non-fractionary meters
# [Done]        300m:   d0 = 23 -> Max.Disp.dist. = 305.6173847
# [Done]        1Km:    d0 = 75 -> Max.Disp.dist. = 996.5784285
# [Done]        2Km :   d0 = 151 -> Max.Disp.dist. = 2006.44
# [Done]        4Km:    d0 = 301 -> Max.Disp.dist. = 3999.601426
# [Done]		6Km:    d0 = 452 -> Max.Disp.dist. = 6006.045996
# [Done]        8Km:    d0 = 602 -> Max.Disp.dist. = 7999.202852
# [Done]        10Km:   d0 = 753 -> Max.Disp.dist. = 10005.64742    #The last one to be processed, which is why this distance is in the script

min_prob = 0.0001 #The minimum probability of dispersal between patches with which max_cost is calculated.


# importing functions
import arcpy, time, os
from arcpy import env
from arcpy.sa import *
import igraph as ig
import numpy as np

# Set the python working directory
os.chdir(wrkspc)

# Setting the initial environment
env.workspace = wrkspc
arcpy.SetProduct('ArcInfo')
arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True


# ##########################################
#Work with saved patch code raster, and save it as a polygon and point shapefile (centroids)

patchCodeRast = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\habPatchCode_Bovar1.tif"
habPatchCode = Raster(patchCodeRast)

arcpy.RasterToPolygon_conversion(habPatchCode, "habPatchCode_Bovar_10KmDispDist.shp", "NO_SIMPLIFY", "VALUE") #Change raster to polygon
arcpy.Dissolve_management("habPatchCode_Bovar_10KmDispDist.shp", "habPatchCode_Bovar_10KmDispDist_Diss.shp", "GRIDCODE", "", "MULTI_PART", "") #Disolve the polygons based on the patch code
arcpy.FeatureToPoint_management("habPatchCode_Bovar_10KmDispDist_Diss.shp", "habPatchCode_Bovar_10KmDispDist_Cntrd.shp", "INSIDE") #Find the centroid of all the patches.
arcpy.AddXY_management ("habPatchCode_Bovar_10KmDispDist_Cntrd.shp") #Add the X and Y coordinates of each centroid to the attribute table
NAhabPatchCode = int(habPatchCode.noDataValue) #Determine the NoData value for later

#Convert the habPatchCode to numpy array
PC = arcpy.RasterToNumPyArray(habPatchCode) #PC includes the ArcGIS patch numbers

#List all the unique patch numbers
patch_array = np.unique(PC)
patch_array = patch_array[patch_array!=NAhabPatchCode].tolist() #Patch array is a 1D vector of all the patch numbers

# ##########################################

# Calculate the max_cost
#log = natural logarithm
max_cost = np.log(min_prob)/np.log(.5)*d0

# Make dictionaries
CD_dict={}  #Dictionary of cost-distance values
COUNT_tab=np.column_stack((patch_array, np.repeat(0,len(patch_array)))) #Dictionary of neighboring patches to each patch
Coord_dict={} #Dictionary containing the coordinates of each patch
edge_tab = [] #An empty array containing an edge in each row

# Add the coordinates of the centroid of each patch to the Coord_dict
fields = ["GRIDCODE","POINT_X","POINT_Y"]
with arcpy.da.SearchCursor("C:\\Users\\damiano\\Documents\\PhD\\Sensitivity_DispDist\\Networks\\habPatchCode_Bovar_10KmDispDist_Cntrd.shp", fields) as cursor:
    for row in cursor:
        Coord_dict["p{0}".format(row[0])] = np.array([row[1],row[2]])

#Build a timer to check progress of the modelling.
counter = 0
perc = 0
interval = patch_array.__len__()/10 #Progress with 10 % increments

# Start the time clock for duration of the calculations
t0 = time.clock()

seen = set() #To prevent double counting of edges (i.e. if A-B is present, B-A will not be included).

for patch in patch_array: #Loop over all the patches
    #For each calculate the cost distance raster (as numpy array) if not in CD_dict
    if "p{0}".format(patch) not in CD_dict:
        outRas = ExtractByAttributes(habPatchCode, "Value = "+str(patch)) #Make a raster layer for the patch in this loop
        CD_dict["p{0}".format(patch)] = arcpy.RasterToNumPyArray(CostDistance(outRas, costSurface, max_cost), nodata_to_value = np.nan) #Calculate the cost-distance raster
        del outRas #Cleanup

    #For each patch place it's reachable neighbors in the dictionary.
    PtPC = np.multiply((np.multiply(CD_dict["p"+str(patch)],0)+1),PC) #Set all cells in a cost-distance array to 1 (excluding all the NoData cells) and multiply it with the PC array. This will result in a subset of the PC array that only included the patches that can be reached by patch.
    reachPatch = np.unique(PtPC[~np.isnan(PtPC)]).astype(int) #Extract the unique patch numers for the patches in PtPC
    reachPatch = reachPatch[reachPatch!=NAhabPatchCode] #Exclude the NoData value of the patch codes.
    for rPatch in reachPatch: #Loop over all the patches in reachPatch
        if patch != rPatch: #Exclude the option patch = patch
            if (rPatch,patch) not in seen: #If the combination in the oposite direction has already been calculated (e.g. A to B), do not calculate it again (e.g. B to A).
                if "p{0}".format(rPatch) not in CD_dict:
                    outRas = ExtractByAttributes(habPatchCode, "Value = "+str(rPatch)) #Make a raster layer for the patch in this loop
                    CD_dict["p{0}".format(rPatch)] = arcpy.RasterToNumPyArray(CostDistance(outRas, costSurface, max_cost), nodata_to_value = np.nan) #Calculate the cost-distance raster
                    del outRas
                c1 = Coord_dict["p"+str(patch)].tolist() #Extract the coordinates patch
                c2 = Coord_dict["p"+str(rPatch)].tolist() #Extract the coordinates of rPatch
                cost = np.nanmin(CD_dict["p"+str(patch)] + CD_dict["p"+str(rPatch)]).tolist()
                prob = np.exp(((cost/d0) * np.log(0.5)))
                edge_tab.append([patch,rPatch,cost,prob,c1[0],c1[1],c2[0],c2[1]]) #Array with columns: codeP1, codeP2, cost, probability, P1_x, P1_y, P2_x, P2_y
                seen.add((patch,rPatch))

    # Update the patch counter that is used to delete cost distance rasters from CD_dict that are no longer used.
    COUNT_tab[:,1] = COUNT_tab[:,1] - 1
    for p in reachPatch:
        if COUNT_tab[np.where(COUNT_tab[:,0] == p),1] < 1:
            COUNT_tab[np.where(COUNT_tab[:,0] == p),1] = 20
    for row in COUNT_tab:
        if row[1] == 0:
            if "p{0}".format(row[0]) in CD_dict:
                del CD_dict["p{0}".format(row[0])]

    # Update timer and print progress
    counter = counter + 1
    if counter == interval:
        counter = 0
        perc = perc + 10
        print("Calculating Cost distances: "+str(perc)+" % done")

# Stop the clock
print("Total time: "+str(int((time.clock() - t0)/60))+" minutes process time") #End the time

edge_tab = np.array(edge_tab) #Numpy array with columns: (0)codeP1, (1)codeP2, (2)cost, (3)probability, (4)P1_x, (5)P1_y, (6)P2_x, (7)P2_y

# Set the python working directory
os.chdir(wrkspc)


# Save the network as a shapefile
save_edge_array_as_shapefile(edge_tab[:,4:8], run_name+".shp", attrib_dat=edge_tab[:,3].tolist())

# Generate an igraph object from the network and write it to a file
CDgraph = edge_table_to_igraph(patch_array, edge_tab, xy_names = (0,1), weight_col = 3)
if CDgraph.ecount() != len(edge_tab) or CDgraph.vcount() != len(patch_array):
    print("Something went wrong with building the igraph graph")
CDgraph.write_graphml(run_name+".graphml")



# ##########################################
# ########### Habitat Availability calculation ##### #
# ##########################################

resRast = 100 #Spatial resolution (m) of the resistance surface and landscape raster.

# importing functions
import arcpy, time, os
from arcpy import env
from arcpy.sa import *
import igraph as ig
import numpy as np
from itertools import combinations

# Set the python working directory
wrkspc = "C:\\Users\\damiano\\Documents\\PhD\\Sensitivity_DispDist\\Networks\\"
os.chdir(wrkspc)
# Setting the initial environment
env.workspace = wrkspc
arcpy.SetProduct('ArcInfo')
arcpy.CheckOutExtension('Spatial')
arcpy.env.overwriteOutput = True


H_CDgraph = ig.Graph.Read_GraphML("Bovar_10KmDispDist.graphml")
habPatchCode = "C:\\Users\\damiano\\Documents\\PhD\\Additional_species_runs\\Network_setup\\habPatchCode_Bovar1.tif" #change for species names or codes

# Make a table of the suitabilities of each habitat patch --> Size ('COUNT'), in the case of the 1st networks
suit_tab = arcpy.da.TableToNumPyArray(habPatchCode,['VALUE', 'COUNT']) #Read the raster attribute table of the habitat patches. This table contains the patch code (VALUE) and the patch size (COUNT)
patch_ID = suit_tab['VALUE'] #Vector of patch codes/IDs
patch_size = (suit_tab['COUNT']*resRast*resRast) #Vector of patch sizes (m2): multiplied by resRast*resRast

patch_suit = patch_size
P_suit_dict={}
for p in range(len(patch_ID)):
    P_suit_dict["p{0}".format(patch_ID[p])] = patch_suit[p] #size

node_suit_dict = P_suit_dict


neighbors = H_CDgraph.neighborhood(vertices=None, order=2, mode="all") #Identify all the neighbors that can be reached in order number of steps.

graph_vnames = H_CDgraph.vs['name']

#Build a timer to check progress of the modelling.
counter = 0
perc = 0
interval = len(neighbors)/10 #Progress with 10 % increments

habAv = np.array([])
for node in neighbors:
    areas = np.array([])
    for neigh in node:
        i = node[0]
        j = neigh
        i_name = graph_vnames[i]
        j_name = graph_vnames[j]
        edges = H_CDgraph.get_shortest_paths(i, to=j, weights=-np.log(H_CDgraph.es['weight']), mode = 'OUT', output="epath")
        if len(edges[0]) == 0:
            areas = np.append(areas, node_suit_dict["p"+str(i_name)])
        else:
            calc = np.prod(H_CDgraph.es[edges[0]]['weight'])*node_suit_dict["p"+str(j_name)]
            areas = np.append(areas, calc)
    habAv = np.append(habAv, sum(areas))
    # Update timer and print progress
    counter = counter + 1
    if counter == interval:
        counter = 0
        perc = perc + 10
        print("Calculating habAv: "+str(perc)+" % done")

np.savetxt('Bovar_10KmDispDist_habAv.txt', habAv, fmt='%.10f')



#### Script for: ##################################################################################################
#### Comparison of the BRT fit of the different network models with different d0s (which give different maximum dispersal distances) #### 
###################################################################################################################
#### Author: Damian O. Ortiz-Rodríguez, Antoine Guisan, Maarten J. van Strien #####################################
#### Article: "Sensitivity of habitat network models to changes in maximum dispersal distance" ####################
###################################################################################################################


library(ggplot2)
library(sp)
library(raster)
library(rgdal)
library(tools)
library(usdm)
library(igraph)
library(plyr)
library(dplyr)
library(data.table)
library(nnet)
library(forcats)
library(dismo)
library(ROCR)
library(rnetcarto)

### Function definition
#moveme : https://stackoverflow.com/a/18540144

moveMe <- function(data, tomove, where = "last", ba = NULL) {
  temp <- setdiff(names(data), tomove)
  x <- switch(
    where,
    first = data[c(tomove, temp)],
    last = data[c(temp, tomove)],
    before = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)-1))]
    },
    after = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)))]
    })
  x
}


#### Species abbreviations ####
# Alobs = Alytes obstetricans
# Bovar = Bombina variegata
# Epcal = Epidalea calamita
# Hyarb = Hyla arborea
# Peagg = Pelophylax agg (P. esculentus + P. lesonae)
# Perid = Pelophylax ridibundus


######################################################################################################################################
#### Preparation for having the network objects and model variables for all the different tried d0's (and hence max. disp. dists) ####
#######################################################################################################################################

#### Load default (species-specific) d0 df's for each species with all BRT predictors and response ####
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/BRTs')

Hyarb_stattest_dispdist_default <- read.csv("Hyarb_stattest_t6.csv") # Also uniform_stattest_6Vt.csv
Bovar_stattest_dispdist_default <- read.csv("Bovar_stattest_t4.csv")
Alobs_stattest_dispdist_default <- read.csv("Alobs_stattest_t4.csv")
Epcal_stattest_dispdist_default <- read.csv("Epcal_stattest_t4.csv")
Peagg_stattest_dispdist_default <- read.csv("Peagg_stattest_t3.csv")
Perid_stattest_dispdist_default <- read.csv("Perid_stattest_t4.csv")


#### Load habpatchcode of each species to make non-default d0 BRT's ####
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Presence_absence')

Pres_abs_allpatches_Hy_arb <- read.csv("Pres_abs_Vs_Vt_Hy_arb_nonpositionaldefinition.csv")
Pres_abs_allpatches_Bo_var <- read.csv("Pres_abs_Vs_Vt_Bo_var.csv")
Pres_abs_allpatches_Al_obs <- read.csv("Pres_abs_Vs_Vt_Al_obs.csv")
Pres_abs_allpatches_Ep_cal <- read.csv("Pres_abs_Vs_Vt_Ep_cal.csv")
Pres_abs_allpatches_Pe_agg <- read.csv("Pres_abs_Vs_Vt_Pe_agg.csv")
Pres_abs_allpatches_Pe_rid <- read.csv("Pres_abs_Vs_Vt_Pe_rid.csv")


### Add HSI to habpatchcode of each species
HSI_Hyarb = read.csv("C:/Users/damiano/Documents/PhD/Automated_Run_BioCHECNET/Corrected/HSI_habPatchCode_NoPseudorep.csv")

setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Network_setup')
HSI_Alobs = read.csv("HSI_habPatchCode_Alobs1.csv")
HSI_Bovar = read.csv("HSI_habPatchCode_Bovar1.csv")
HSI_Epcal = read.csv("HSI_habPatchCode_Epcal1.csv")
HSI_Peagg = read.csv("HSI_habPatchCode_Peagg1.csv")
HSI_Perid = read.csv("HSI_habPatchCode_Perid1.csv")

Pres_abs_allpatches_Hy_arb$HSI <- HSI_Hyarb$"MEAN"
Pres_abs_allpatches_Al_obs$HSI <- HSI_Alobs$"MEAN"
Pres_abs_allpatches_Bo_var$HSI <- HSI_Bovar$"MEAN"
Pres_abs_allpatches_Ep_cal$HSI <- HSI_Epcal$"MEAN"
Pres_abs_allpatches_Pe_agg$HSI <- HSI_Peagg$"MEAN"
Pres_abs_allpatches_Pe_rid$HSI <- HSI_Perid$"MEAN"


#### Define occurrence_state (pres_abs) for each species####
### Classify the dataframe in subsets that classify all the patches in presences, absences or questionmarks/NotDefined (NA) 

#subset of questionmarks, #All the records, except the ones that comply with the following commands
Pres_abs_allpatches_Hy_arb[,"pres_abs"] = NA 
Pres_abs_allpatches_Al_obs[,"pres_abs"] = NA 
Pres_abs_allpatches_Bo_var[,"pres_abs"] = NA 
Pres_abs_allpatches_Ep_cal[,"pres_abs"] = NA 
Pres_abs_allpatches_Pe_agg[,"pres_abs"] = NA 
Pres_abs_allpatches_Pe_rid[,"pres_abs"] = NA 

#Subset of likely absences (With threshold defined by pres_abs decider plot ))
Pres_abs_allpatches_Hy_arb[which(Pres_abs_allpatches_Hy_arb$V_t >= 6 & Pres_abs_allpatches_Hy_arb$V_Hy_arb ==0),"pres_abs"] = 0
Pres_abs_allpatches_Bo_var[which(Pres_abs_allpatches_Bo_var$V_t >= 4 & Pres_abs_allpatches_Bo_var$V_Bo_var ==0),"pres_abs"] = 0
Pres_abs_allpatches_Al_obs[which(Pres_abs_allpatches_Al_obs$V_t >= 4 & Pres_abs_allpatches_Al_obs$V_Al_obs ==0),"pres_abs"] = 0
Pres_abs_allpatches_Ep_cal[which(Pres_abs_allpatches_Ep_cal$V_t >= 4 & Pres_abs_allpatches_Ep_cal$V_Ep_cal ==0),"pres_abs"] = 0
Pres_abs_allpatches_Pe_agg[which(Pres_abs_allpatches_Pe_agg$V_t >= 3 & Pres_abs_allpatches_Pe_agg$V_Pe_agg ==0),"pres_abs"] = 0
Pres_abs_allpatches_Pe_rid[which(Pres_abs_allpatches_Pe_rid$V_t >= 4 & Pres_abs_allpatches_Pe_rid$V_Pe_rid ==0),"pres_abs"] = 0

#Subset of confirmed presences
Pres_abs_allpatches_Hy_arb[which(Pres_abs_allpatches_Hy_arb$V_Hy_arb > 0),"pres_abs"] = 1
Pres_abs_allpatches_Bo_var[which(Pres_abs_allpatches_Bo_var$V_Bo_var > 0),"pres_abs"] = 1
Pres_abs_allpatches_Al_obs[which(Pres_abs_allpatches_Al_obs$V_Al_obs > 0),"pres_abs"] = 1
Pres_abs_allpatches_Ep_cal[which(Pres_abs_allpatches_Ep_cal$V_Ep_cal > 0),"pres_abs"] = 1
Pres_abs_allpatches_Pe_agg[which(Pres_abs_allpatches_Pe_agg$V_Pe_agg > 0),"pres_abs"] = 1
Pres_abs_allpatches_Pe_rid[which(Pres_abs_allpatches_Pe_rid$V_Pe_rid > 0),"pres_abs"] = 1

#### Display the 'occurrence state' (0/1/NA) of every patch as defined above (pres_abs) ####
#Get total amount of values for presence and for absence
table(Pres_abs_allpatches_Hy_arb$pres_abs)
table(Pres_abs_allpatches_Bo_var$pres_abs)
table(Pres_abs_allpatches_Al_obs$pres_abs)
table(Pres_abs_allpatches_Ep_cal$pres_abs)
table(Pres_abs_allpatches_Pe_agg$pres_abs)
table(Pres_abs_allpatches_Pe_rid$pres_abs)


# Save habpatchcode df w/HSI and pres_abs
setwd('C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/Network_setup')

write.csv(Pres_abs_allpatches_Hy_arb, file = "Pres_abs_allpatches_Hy_arb_w_HSI_presabs.csv")
write.csv(Pres_abs_allpatches_Al_obs, file = "Pres_abs_allpatches_Al_obs_w_HSI_presabs.csv")
write.csv(Pres_abs_allpatches_Bo_var, file = "Pres_abs_allpatches_Bo_var_w_HSI_presabs.csv")
write.csv(Pres_abs_allpatches_Ep_cal, file = "Pres_abs_allpatches_Ep_cal_w_HSI_presabs.csv")
write.csv(Pres_abs_allpatches_Pe_agg, file = "Pres_abs_allpatches_Pe_agg_w_HSI_presabs.csv")
write.csv(Pres_abs_allpatches_Pe_rid, file = "Pres_abs_allpatches_Pe_rid_w_HSI_presabs.csv")


#### Rename ID & area columns### Change also V_Ge_spe for each species to V_s for all,
#to have it more generic
Pres_abs_allpatches_Hy_arb <- setnames(Pres_abs_allpatches_Hy_arb, "Value", "PatchID")
Pres_abs_allpatches_Hy_arb <- setnames(Pres_abs_allpatches_Hy_arb, "Count", "Patch_Area")
Pres_abs_allpatches_Hy_arb <- setnames(Pres_abs_allpatches_Hy_arb, "V_Hy_arb", "V_s")

Pres_abs_allpatches_Al_obs <- setnames(Pres_abs_allpatches_Al_obs, "Value", "PatchID")
Pres_abs_allpatches_Al_obs <- setnames(Pres_abs_allpatches_Al_obs, "Count", "Patch_Area")
Pres_abs_allpatches_Al_obs <- setnames(Pres_abs_allpatches_Al_obs, "V_Al_obs", "V_s")

Pres_abs_allpatches_Bo_var <- setnames(Pres_abs_allpatches_Bo_var, "Value", "PatchID")
Pres_abs_allpatches_Bo_var <- setnames(Pres_abs_allpatches_Bo_var, "Count", "Patch_Area")
Pres_abs_allpatches_Bo_var <- setnames(Pres_abs_allpatches_Bo_var, "V_Bo_var", "V_s")

Pres_abs_allpatches_Ep_cal <- setnames(Pres_abs_allpatches_Ep_cal, "Value", "PatchID")
Pres_abs_allpatches_Ep_cal <- setnames(Pres_abs_allpatches_Ep_cal, "Count", "Patch_Area")
Pres_abs_allpatches_Ep_cal <- setnames(Pres_abs_allpatches_Ep_cal, "V_Ep_cal", "V_s")

Pres_abs_allpatches_Pe_agg <- setnames(Pres_abs_allpatches_Pe_agg, "Value", "PatchID")
Pres_abs_allpatches_Pe_agg <- setnames(Pres_abs_allpatches_Pe_agg, "Count", "Patch_Area")
Pres_abs_allpatches_Pe_agg <- setnames(Pres_abs_allpatches_Pe_agg, "V_Pe_agg", "V_s")

Pres_abs_allpatches_Pe_rid <- setnames(Pres_abs_allpatches_Pe_rid, "Value", "PatchID")
Pres_abs_allpatches_Pe_rid <- setnames(Pres_abs_allpatches_Pe_rid, "Count", "Patch_Area")
Pres_abs_allpatches_Pe_rid <- setnames(Pres_abs_allpatches_Pe_rid, "V_Pe_rid", "V_s")

Pres_abs_allpatches_Pe_rid <- moveMe(Pres_abs_allpatches_Pe_rid, "PatchID", "first")
names(Pres_abs_allpatches_Pe_rid)


#### Load networks #########################################################################

#Original species-specific d0 (default) networks
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Networks')

Alobs_graph_default_DispDist <- read_graph("Alobs.graphml", format = "graphml") 
Bovar_graph_default_DispDist <- read_graph("Bovar.graphml", format = "graphml") 
Hyarb_graph_default_DispDist <- read_graph("uniform.graphml", format = "graphml")
Epcal_graph_default_DispDist <- read_graph("Epcal.graphml", format = "graphml")
Peagg_graph_default_DispDist <- read_graph("Peagg.graphml", format = "graphml")
Perid_graph_default_DispDist <- read_graph("Perid.graphml", format = "graphml")


#New networks with d0 variations
setwd('C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/Networks')

#300 m disp dist networks
Hyarb_graph_300m_DispDist <- read_graph("Hyarb_300mDispDist.graphml", format = "graphml")
Bovar_graph_300m_DispDist <- read_graph("Bovar_300mDispDist.graphml", format = "graphml")
Alobs_graph_300m_DispDist <- read_graph("Alobs_300mDispDist.graphml", format = "graphml")
Epcal_graph_300m_DispDist <- read_graph("Epcal_300mDispDist.graphml", format = "graphml")
Peagg_graph_300m_DispDist <- read_graph("Peagg_300mDispDist.graphml", format = "graphml")
Perid_graph_300m_DispDist <- read_graph("Perid_300mDispDist.graphml", format = "graphml")

#1 Km disp dist networks
Hyarb_graph_1Km_DispDist <- read_graph("Hyarb_1KmDispDist.graphml", format = "graphml")
Bovar_graph_1Km_DispDist <- read_graph("Bovar_1KmDispDist.graphml", format = "graphml")
Alobs_graph_1Km_DispDist <- read_graph("Alobs_1KmDispDist.graphml", format = "graphml")
Epcal_graph_1Km_DispDist <- read_graph("Epcal_1KmDispDist.graphml", format = "graphml")
Peagg_graph_1Km_DispDist <- read_graph("Peagg_1KmDispDist.graphml", format = "graphml")
Perid_graph_1Km_DispDist <- read_graph("Perid_1KmDispDist.graphml", format = "graphml")

#2 Km disp dist networks
Hyarb_graph_2Km_DispDist <- read_graph("Hyarb_2KmDispDist.graphml", format = "graphml")
Bovar_graph_2Km_DispDist <- read_graph("Bovar_2KmDispDist.graphml", format = "graphml")
Alobs_graph_2Km_DispDist <- read_graph("Alobs_2KmDispDist.graphml", format = "graphml")
Epcal_graph_2Km_DispDist <- read_graph("Epcal_2KmDispDist.graphml", format = "graphml")
Peagg_graph_2Km_DispDist <- read_graph("Peagg_2KmDispDist.graphml", format = "graphml")
Perid_graph_2Km_DispDist <- read_graph("Perid_2KmDispDist.graphml", format = "graphml")

#4 Km disp dist networks
Hyarb_graph_4Km_DispDist <- read_graph("Hyarb_4KmDispDist.graphml", format = "graphml")
# Bovar_graph_4Km_DispDist <- read_graph("Bovar_4KmDispDist.graphml", format = "graphml") #4 km is the default disp dist of this species, so no need to do it twice
Alobs_graph_4Km_DispDist <- read_graph("Alobs_4KmDispDist.graphml", format = "graphml")
Epcal_graph_4Km_DispDist <- read_graph("Epcal_4KmDispDist.graphml", format = "graphml")
Peagg_graph_4Km_DispDist <- read_graph("Peagg_4KmDispDist.graphml", format = "graphml")
Perid_graph_4Km_DispDist <- read_graph("Perid_4KmDispDist.graphml", format = "graphml")

#6 Km disp dist networks
Hyarb_graph_6Km_DispDist <- read_graph("Hyarb_6KmDispDist.graphml", format = "graphml")
Bovar_graph_6Km_DispDist <- read_graph("Bovar_6KmDispDist.graphml", format = "graphml")
Alobs_graph_6Km_DispDist <- read_graph("Alobs_6KmDispDist.graphml", format = "graphml")
Epcal_graph_6Km_DispDist <- read_graph("Epcal_6KmDispDist.graphml", format = "graphml")
Peagg_graph_6Km_DispDist <- read_graph("Peagg_6KmDispDist.graphml", format = "graphml")
Perid_graph_6Km_DispDist <- read_graph("Perid_6KmDispDist.graphml", format = "graphml")

#8 Km disp dist networks
Hyarb_graph_8Km_DispDist <- read_graph("Hyarb_8KmDispDist.graphml", format = "graphml")
Bovar_graph_8Km_DispDist <- read_graph("Bovar_8KmDispDist.graphml", format = "graphml") 
Alobs_graph_8Km_DispDist <- read_graph("Alobs_8KmDispDist.graphml", format = "graphml")
Epcal_graph_8Km_DispDist <- read_graph("Epcal_8KmDispDist.graphml", format = "graphml")
Peagg_graph_8Km_DispDist <- read_graph("Peagg_8KmDispDist.graphml", format = "graphml")
Perid_graph_8Km_DispDist <- read_graph("Perid_8KmDispDist.graphml", format = "graphml")

#10 Km disp dist networks
Hyarb_graph_10Km_DispDist <- read_graph("Hyarb_10KmDispDist.graphml", format = "graphml")
Bovar_graph_10Km_DispDist <- read_graph("Bovar_10KmDispDist.graphml", format = "graphml")
Alobs_graph_10Km_DispDist <- read_graph("Alobs_10KmDispDist.graphml", format = "graphml")
Epcal_graph_10Km_DispDist <- read_graph("Epcal_10KmDispDist.graphml", format = "graphml")
Peagg_graph_10Km_DispDist <- read_graph("Peagg_10KmDispDist.graphml", format = "graphml")
Perid_graph_10Km_DispDist <- read_graph("Perid_10KmDispDist.graphml", format = "graphml")


#Look at vertices and edges of graph
#300m disp.dist. networks
E(Hyarb_graph_300m_DispDist)
V(Hyarb_graph_300m_DispDist)
E(Bovar_graph_300m_DispDist)
V(Bovar_graph_300m_DispDist)
E(Alobs_graph_300m_DispDist)
V(Alobs_graph_300m_DispDist)
E(Epcal_graph_300m_DispDist)
V(Epcal_graph_300m_DispDist)
E(Peagg_graph_300m_DispDist)
V(Peagg_graph_300m_DispDist)
E(Perid_graph_300m_DispDist)
V(Perid_graph_300m_DispDist)

#1 Km disp dist networks
E(Hyarb_graph_1Km_DispDist)
V(Hyarb_graph_1Km_DispDist)
E(Bovar_graph_1Km_DispDist)
V(Bovar_graph_1Km_DispDist)
E(Alobs_graph_1Km_DispDist)
V(Alobs_graph_1Km_DispDist)
E(Epcal_graph_1Km_DispDist)
V(Epcal_graph_1Km_DispDist)
E(Peagg_graph_1Km_DispDist)
V(Peagg_graph_1Km_DispDist)
E(Perid_graph_1Km_DispDist)
V(Perid_graph_1Km_DispDist)

#2Km disp.dist. networks
E(Hyarb_graph_2Km_DispDist)
V(Hyarb_graph_2Km_DispDist)
E(Bovar_graph_2Km_DispDist)
V(Bovar_graph_2Km_DispDist)
E(Alobs_graph_2Km_DispDist)
V(Alobs_graph_2Km_DispDist)
E(Epcal_graph_2Km_DispDist)
V(Epcal_graph_2Km_DispDist)
E(Peagg_graph_2Km_DispDist)
V(Peagg_graph_2Km_DispDist)
E(Perid_graph_2Km_DispDist)
V(Perid_graph_2Km_DispDist)

#4 Km disp dist networks
E(Hyarb_graph_4Km_DispDist)
V(Hyarb_graph_4Km_DispDist)
# E(Bovar_graph_4Km_DispDist)
# V(Bovar_graph_4Km_DispDist)
E(Alobs_graph_4Km_DispDist)
V(Alobs_graph_4Km_DispDist)
E(Epcal_graph_4Km_DispDist)
V(Epcal_graph_4Km_DispDist)
E(Peagg_graph_4Km_DispDist)
V(Peagg_graph_4Km_DispDist)
E(Perid_graph_4Km_DispDist)
V(Perid_graph_4Km_DispDist)

#6 Km disp dist networks
E(Hyarb_graph_6Km_DispDist)
V(Hyarb_graph_6Km_DispDist)
E(Bovar_graph_6Km_DispDist)
V(Bovar_graph_6Km_DispDist)
E(Alobs_graph_6Km_DispDist)
V(Alobs_graph_6Km_DispDist)
E(Epcal_graph_6Km_DispDist)
V(Epcal_graph_6Km_DispDist)
E(Peagg_graph_6Km_DispDist)
V(Peagg_graph_6Km_DispDist)
E(Perid_graph_6Km_DispDist)
V(Perid_graph_6Km_DispDist)

#8 Km disp dist networks
E(Hyarb_graph_8Km_DispDist)
V(Hyarb_graph_8Km_DispDist)
E(Bovar_graph_8Km_DispDist)
V(Bovar_graph_8Km_DispDist)
E(Alobs_graph_8Km_DispDist)
V(Alobs_graph_8Km_DispDist)
E(Epcal_graph_8Km_DispDist)
V(Epcal_graph_8Km_DispDist)
E(Peagg_graph_8Km_DispDist)
V(Peagg_graph_8Km_DispDist)
E(Perid_graph_8Km_DispDist)
V(Perid_graph_8Km_DispDist)

#10Km disp.dist. networks
E(Hyarb_graph_10Km_DispDist)
V(Hyarb_graph_10Km_DispDist)
E(Bovar_graph_10Km_DispDist)
V(Bovar_graph_10Km_DispDist)
E(Alobs_graph_10Km_DispDist)
V(Alobs_graph_10Km_DispDist)
E(Epcal_graph_10Km_DispDist)
V(Epcal_graph_10Km_DispDist)
E(Peagg_graph_10Km_DispDist) 
V(Peagg_graph_10Km_DispDist) 
E(Perid_graph_10Km_DispDist) 
V(Perid_graph_10Km_DispDist) 



#### Get number of components per network #################################################
#default = species-specific maximum dispersal distance from literature

count_components(Alobs_graph_default_DispDist)
count_components(Alobs_graph_300m_DispDist)
count_components(Alobs_graph_1Km_DispDist)
count_components(Alobs_graph_2Km_DispDist)
count_components(Alobs_graph_4Km_DispDist)
count_components(Alobs_graph_6Km_DispDist)
count_components(Alobs_graph_8Km_DispDist)
count_components(Alobs_graph_10Km_DispDist)

count_components(Bovar_graph_default_DispDist)
count_components(Bovar_graph_300m_DispDist)
count_components(Bovar_graph_1Km_DispDist)
count_components(Bovar_graph_2Km_DispDist)
# count_components(Bovar_graph_4Km_DispDist)
count_components(Bovar_graph_6Km_DispDist)
count_components(Bovar_graph_8Km_DispDist)
count_components(Bovar_graph_10Km_DispDist)

count_components(Epcal_graph_default_DispDist)
count_components(Epcal_graph_300m_DispDist)
count_components(Epcal_graph_1Km_DispDist)
count_components(Epcal_graph_2Km_DispDist)
count_components(Epcal_graph_4Km_DispDist)
count_components(Epcal_graph_6Km_DispDist)
count_components(Epcal_graph_8Km_DispDist)
count_components(Epcal_graph_10Km_DispDist)

count_components(Hyarb_graph_default_DispDist)
count_components(Hyarb_graph_300m_DispDist)
count_components(Hyarb_graph_1Km_DispDist)
count_components(Hyarb_graph_2Km_DispDist)
count_components(Hyarb_graph_4Km_DispDist)
count_components(Hyarb_graph_6Km_DispDist)
count_components(Hyarb_graph_8Km_DispDist)
count_components(Hyarb_graph_10Km_DispDist)

count_components(Peagg_graph_default_DispDist)
count_components(Peagg_graph_300m_DispDist)
count_components(Peagg_graph_1Km_DispDist)
count_components(Peagg_graph_2Km_DispDist)
count_components(Peagg_graph_4Km_DispDist)
count_components(Peagg_graph_6Km_DispDist)
count_components(Peagg_graph_8Km_DispDist)
count_components(Peagg_graph_10Km_DispDist)

count_components(Perid_graph_default_DispDist)
count_components(Perid_graph_300m_DispDist)
count_components(Perid_graph_1Km_DispDist)
count_components(Perid_graph_2Km_DispDist)
count_components(Perid_graph_4Km_DispDist)
count_components(Perid_graph_6Km_DispDist)
count_components(Perid_graph_8Km_DispDist)
count_components(Perid_graph_10Km_DispDist)

count_components(Bovar_graph_default_DispDist)
count_components(Hyarb_graph_default_DispDist)
count_components(Epcal_graph_default_DispDist)
count_components(Peagg_graph_default_DispDist)
count_components(Perid_graph_default_DispDist)


###################################################################################################
#### Calculate igraph network properties for newly developed networks (non-default disp_dists) ####
#### get metrics as attributes and into a data.frame ####

## Hyarb
#2Km
Hyarb_graph_2Km_DispDist_metrics <- data.frame(
  PatchID = V(Hyarb_graph_2Km_DispDist)$name,
  deg=degree(Hyarb_graph_2Km_DispDist),
  strength = strength(Hyarb_graph_2Km_DispDist, vids = V(Hyarb_graph_2Km_DispDist), mode = "all", loops = FALSE, weights = E(Hyarb_graph_2Km_DispDist)$weight),
  EgoSize = ego_size(Hyarb_graph_2Km_DispDist, order = 3, nodes = V(Hyarb_graph_2Km_DispDist))
)
#300m
Hyarb_graph_300m_DispDist_metrics <- data.frame(
  PatchID = V(Hyarb_graph_300m_DispDist)$name,
  deg=degree(Hyarb_graph_300m_DispDist),
  strength = strength(Hyarb_graph_300m_DispDist, vids = V(Hyarb_graph_300m_DispDist), mode = "all", loops = FALSE, weights = E(Hyarb_graph_300m_DispDist)$weight),
  EgoSize = ego_size(Hyarb_graph_300m_DispDist, order = 3, nodes = V(Hyarb_graph_300m_DispDist))
)
#10Km
Hyarb_graph_10Km_DispDist_metrics <- data.frame(
  PatchID = V(Hyarb_graph_10Km_DispDist)$name,
  deg=degree(Hyarb_graph_10Km_DispDist),
  strength = strength(Hyarb_graph_10Km_DispDist, vids = V(Hyarb_graph_10Km_DispDist), mode = "all", loops = FALSE, weights = E(Hyarb_graph_10Km_DispDist)$weight),
  EgoSize = ego_size(Hyarb_graph_10Km_DispDist, order = 3, nodes = V(Hyarb_graph_10Km_DispDist))
)
#1Km
Hyarb_graph_1Km_DispDist_metrics <- data.frame(
  PatchID = V(Hyarb_graph_1Km_DispDist)$name,
  deg=degree(Hyarb_graph_1Km_DispDist),
  strength = strength(Hyarb_graph_1Km_DispDist, vids = V(Hyarb_graph_1Km_DispDist), mode = "all", loops = FALSE, weights = E(Hyarb_graph_1Km_DispDist)$weight),
  EgoSize = ego_size(Hyarb_graph_1Km_DispDist, order = 3, nodes = V(Hyarb_graph_1Km_DispDist))
)
#6Km
Hyarb_graph_6Km_DispDist_metrics <- data.frame(
  PatchID = V(Hyarb_graph_6Km_DispDist)$name,
  deg=degree(Hyarb_graph_6Km_DispDist),
  strength = strength(Hyarb_graph_6Km_DispDist, vids = V(Hyarb_graph_6Km_DispDist), mode = "all", loops = FALSE, weights = E(Hyarb_graph_6Km_DispDist)$weight),
  EgoSize = ego_size(Hyarb_graph_6Km_DispDist, order = 3, nodes = V(Hyarb_graph_6Km_DispDist))
)
#4Km
Hyarb_graph_4Km_DispDist_metrics <- data.frame(
  PatchID = V(Hyarb_graph_4Km_DispDist)$name,
  deg=degree(Hyarb_graph_4Km_DispDist),
  strength = strength(Hyarb_graph_4Km_DispDist, vids = V(Hyarb_graph_4Km_DispDist), mode = "all", loops = FALSE, weights = E(Hyarb_graph_4Km_DispDist)$weight),
  EgoSize = ego_size(Hyarb_graph_4Km_DispDist, order = 3, nodes = V(Hyarb_graph_4Km_DispDist))
)
#8Km
Hyarb_graph_8Km_DispDist_metrics <- data.frame(
  PatchID = V(Hyarb_graph_8Km_DispDist)$name,
  deg=degree(Hyarb_graph_8Km_DispDist),
  strength = strength(Hyarb_graph_8Km_DispDist, vids = V(Hyarb_graph_8Km_DispDist), mode = "all", loops = FALSE, weights = E(Hyarb_graph_8Km_DispDist)$weight),
  EgoSize = ego_size(Hyarb_graph_8Km_DispDist, order = 3, nodes = V(Hyarb_graph_8Km_DispDist))
)

head(Hyarb_graph_2Km_DispDist_metrics)
head(Hyarb_graph_300m_DispDist_metrics)
head(Hyarb_graph_10Km_DispDist_metrics)

head(Hyarb_graph_1Km_DispDist_metrics)
head(Hyarb_graph_6Km_DispDist_metrics)
head(Hyarb_graph_4Km_DispDist_metrics)
head(Hyarb_graph_8Km_DispDist_metrics)


## Bovar
#2Km
Bovar_graph_2Km_DispDist_metrics <- data.frame(
  PatchID = V(Bovar_graph_2Km_DispDist)$name,
  deg=degree(Bovar_graph_2Km_DispDist),
  strength = strength(Bovar_graph_2Km_DispDist, vids = V(Bovar_graph_2Km_DispDist), mode = "all", loops = FALSE, weights = E(Bovar_graph_2Km_DispDist)$weight),
  EgoSize = ego_size(Bovar_graph_2Km_DispDist, order = 3, nodes = V(Bovar_graph_2Km_DispDist))
)
#300m
Bovar_graph_300m_DispDist_metrics <- data.frame(
  PatchID = V(Bovar_graph_300m_DispDist)$name,
  deg=degree(Bovar_graph_300m_DispDist),
  strength = strength(Bovar_graph_300m_DispDist, vids = V(Bovar_graph_300m_DispDist), mode = "all", loops = FALSE, weights = E(Bovar_graph_300m_DispDist)$weight),
  EgoSize = ego_size(Bovar_graph_300m_DispDist, order = 3, nodes = V(Bovar_graph_300m_DispDist))
)
#10Km
Bovar_graph_10Km_DispDist_metrics <- data.frame(
  PatchID = V(Bovar_graph_10Km_DispDist)$name,
  deg=degree(Bovar_graph_10Km_DispDist),
  strength = strength(Bovar_graph_10Km_DispDist, vids = V(Bovar_graph_10Km_DispDist), mode = "all", loops = FALSE, weights = E(Bovar_graph_10Km_DispDist)$weight),
  EgoSize = ego_size(Bovar_graph_10Km_DispDist, order = 3, nodes = V(Bovar_graph_10Km_DispDist))
)
#1Km
Bovar_graph_1Km_DispDist_metrics <- data.frame(
  PatchID = V(Bovar_graph_1Km_DispDist)$name,
  deg=degree(Bovar_graph_1Km_DispDist),
  strength = strength(Bovar_graph_1Km_DispDist, vids = V(Bovar_graph_1Km_DispDist), mode = "all", loops = FALSE, weights = E(Bovar_graph_1Km_DispDist)$weight),
  EgoSize = ego_size(Bovar_graph_1Km_DispDist, order = 3, nodes = V(Bovar_graph_1Km_DispDist))
)
#6Km
Bovar_graph_6Km_DispDist_metrics <- data.frame(
  PatchID = V(Bovar_graph_6Km_DispDist)$name,
  deg=degree(Bovar_graph_6Km_DispDist),
  strength = strength(Bovar_graph_6Km_DispDist, vids = V(Bovar_graph_6Km_DispDist), mode = "all", loops = FALSE, weights = E(Bovar_graph_6Km_DispDist)$weight),
  EgoSize = ego_size(Bovar_graph_6Km_DispDist, order = 3, nodes = V(Bovar_graph_6Km_DispDist))
)
# #4Km #Default disp dist for this sp. - not required
# Bovar_graph_4Km_DispDist_metrics <- data.frame(
#   PatchID = V(Bovar_graph_4Km_DispDist)$name,
#   deg=degree(Bovar_graph_4Km_DispDist),
#   strength = strength(Bovar_graph_4Km_DispDist, vids = V(Bovar_graph_4Km_DispDist), mode = "all", loops = FALSE, weights = E(Bovar_graph_4Km_DispDist)$weight),
#   EgoSize = ego_size(Bovar_graph_4Km_DispDist, order = 3, nodes = V(Bovar_graph_4Km_DispDist))
# )
# 
#8Km
Bovar_graph_8Km_DispDist_metrics <- data.frame(
  PatchID = V(Bovar_graph_8Km_DispDist)$name,
  deg=degree(Bovar_graph_8Km_DispDist),
  strength = strength(Bovar_graph_8Km_DispDist, vids = V(Bovar_graph_8Km_DispDist), mode = "all", loops = FALSE, weights = E(Bovar_graph_8Km_DispDist)$weight),
  EgoSize = ego_size(Bovar_graph_8Km_DispDist, order = 3, nodes = V(Bovar_graph_8Km_DispDist))
)


##Alobs
#2Km
Alobs_graph_2Km_DispDist_metrics <- data.frame(
  PatchID = V(Alobs_graph_2Km_DispDist)$name,
  deg=degree(Alobs_graph_2Km_DispDist),
  strength = strength(Alobs_graph_2Km_DispDist, vids = V(Alobs_graph_2Km_DispDist), mode = "all", loops = FALSE, weights = E(Alobs_graph_2Km_DispDist)$weight),
  EgoSize = ego_size(Alobs_graph_2Km_DispDist, order = 3, nodes = V(Alobs_graph_2Km_DispDist))
)
#300m
Alobs_graph_300m_DispDist_metrics <- data.frame(
  PatchID = V(Alobs_graph_300m_DispDist)$name,
  deg=degree(Alobs_graph_300m_DispDist),
  strength = strength(Alobs_graph_300m_DispDist, vids = V(Alobs_graph_300m_DispDist), mode = "all", loops = FALSE, weights = E(Alobs_graph_300m_DispDist)$weight),
  EgoSize = ego_size(Alobs_graph_300m_DispDist, order = 3, nodes = V(Alobs_graph_300m_DispDist))
)
#10Km
Alobs_graph_10Km_DispDist_metrics <- data.frame(
  PatchID = V(Alobs_graph_10Km_DispDist)$name,
  deg=degree(Alobs_graph_10Km_DispDist),
  strength = strength(Alobs_graph_10Km_DispDist, vids = V(Alobs_graph_10Km_DispDist), mode = "all", loops = FALSE, weights = E(Alobs_graph_10Km_DispDist)$weight),
  EgoSize = ego_size(Alobs_graph_10Km_DispDist, order = 3, nodes = V(Alobs_graph_10Km_DispDist))
)
#1Km
Alobs_graph_1Km_DispDist_metrics <- data.frame(
  PatchID = V(Alobs_graph_1Km_DispDist)$name,
  deg=degree(Alobs_graph_1Km_DispDist),
  strength = strength(Alobs_graph_1Km_DispDist, vids = V(Alobs_graph_1Km_DispDist), mode = "all", loops = FALSE, weights = E(Alobs_graph_1Km_DispDist)$weight),
  EgoSize = ego_size(Alobs_graph_1Km_DispDist, order = 3, nodes = V(Alobs_graph_1Km_DispDist))
)
#6Km
Alobs_graph_6Km_DispDist_metrics <- data.frame(
  PatchID = V(Alobs_graph_6Km_DispDist)$name,
  deg=degree(Alobs_graph_6Km_DispDist),
  strength = strength(Alobs_graph_6Km_DispDist, vids = V(Alobs_graph_6Km_DispDist), mode = "all", loops = FALSE, weights = E(Alobs_graph_6Km_DispDist)$weight),
  EgoSize = ego_size(Alobs_graph_6Km_DispDist, order = 3, nodes = V(Alobs_graph_6Km_DispDist))
)
#4Km
Alobs_graph_4Km_DispDist_metrics <- data.frame(
  PatchID = V(Alobs_graph_4Km_DispDist)$name,
  deg=degree(Alobs_graph_4Km_DispDist),
  strength = strength(Alobs_graph_4Km_DispDist, vids = V(Alobs_graph_4Km_DispDist), mode = "all", loops = FALSE, weights = E(Alobs_graph_4Km_DispDist)$weight),
  EgoSize = ego_size(Alobs_graph_4Km_DispDist, order = 3, nodes = V(Alobs_graph_4Km_DispDist))
)
#8Km
Alobs_graph_8Km_DispDist_metrics <- data.frame(
  PatchID = V(Alobs_graph_8Km_DispDist)$name,
  deg=degree(Alobs_graph_8Km_DispDist),
  strength = strength(Alobs_graph_8Km_DispDist, vids = V(Alobs_graph_8Km_DispDist), mode = "all", loops = FALSE, weights = E(Alobs_graph_8Km_DispDist)$weight),
  EgoSize = ego_size(Alobs_graph_8Km_DispDist, order = 3, nodes = V(Alobs_graph_8Km_DispDist))
)


## Epcal
#2Km
Epcal_graph_2Km_DispDist_metrics <- data.frame(
  PatchID = V(Epcal_graph_2Km_DispDist)$name,
  deg=degree(Epcal_graph_2Km_DispDist),
  strength = strength(Epcal_graph_2Km_DispDist, vids = V(Epcal_graph_2Km_DispDist), mode = "all", loops = FALSE, weights = E(Epcal_graph_2Km_DispDist)$weight),
  EgoSize = ego_size(Epcal_graph_2Km_DispDist, order = 3, nodes = V(Epcal_graph_2Km_DispDist))
)
#300m
Epcal_graph_300m_DispDist_metrics <- data.frame(
  PatchID = V(Epcal_graph_300m_DispDist)$name,
  deg=degree(Epcal_graph_300m_DispDist),
  strength = strength(Epcal_graph_300m_DispDist, vids = V(Epcal_graph_300m_DispDist), mode = "all", loops = FALSE, weights = E(Epcal_graph_300m_DispDist)$weight),
  EgoSize = ego_size(Epcal_graph_300m_DispDist, order = 3, nodes = V(Epcal_graph_300m_DispDist))
)
#10Km
Epcal_graph_10Km_DispDist_metrics <- data.frame(
  PatchID = V(Epcal_graph_10Km_DispDist)$name,
  deg=degree(Epcal_graph_10Km_DispDist),
  strength = strength(Epcal_graph_10Km_DispDist, vids = V(Epcal_graph_10Km_DispDist), mode = "all", loops = FALSE, weights = E(Epcal_graph_10Km_DispDist)$weight),
  EgoSize = ego_size(Epcal_graph_10Km_DispDist, order = 3, nodes = V(Epcal_graph_10Km_DispDist))
)

#1Km
Epcal_graph_1Km_DispDist_metrics <- data.frame(
  PatchID = V(Epcal_graph_1Km_DispDist)$name,
  deg=degree(Epcal_graph_1Km_DispDist),
  strength = strength(Epcal_graph_1Km_DispDist, vids = V(Epcal_graph_1Km_DispDist), mode = "all", loops = FALSE, weights = E(Epcal_graph_1Km_DispDist)$weight),
  EgoSize = ego_size(Epcal_graph_1Km_DispDist, order = 3, nodes = V(Epcal_graph_1Km_DispDist))
)
#6Km
Epcal_graph_6Km_DispDist_metrics <- data.frame(
  PatchID = V(Epcal_graph_6Km_DispDist)$name,
  deg=degree(Epcal_graph_6Km_DispDist),
  strength = strength(Epcal_graph_6Km_DispDist, vids = V(Epcal_graph_6Km_DispDist), mode = "all", loops = FALSE, weights = E(Epcal_graph_6Km_DispDist)$weight),
  EgoSize = ego_size(Epcal_graph_6Km_DispDist, order = 3, nodes = V(Epcal_graph_6Km_DispDist))
)
#4Km
Epcal_graph_4Km_DispDist_metrics <- data.frame(
  PatchID = V(Epcal_graph_4Km_DispDist)$name,
  deg=degree(Epcal_graph_4Km_DispDist),
  strength = strength(Epcal_graph_4Km_DispDist, vids = V(Epcal_graph_4Km_DispDist), mode = "all", loops = FALSE, weights = E(Epcal_graph_4Km_DispDist)$weight),
  EgoSize = ego_size(Epcal_graph_4Km_DispDist, order = 3, nodes = V(Epcal_graph_4Km_DispDist))
)
#8Km
Epcal_graph_8Km_DispDist_metrics <- data.frame(
  PatchID = V(Epcal_graph_8Km_DispDist)$name,
  deg=degree(Epcal_graph_8Km_DispDist),
  strength = strength(Epcal_graph_8Km_DispDist, vids = V(Epcal_graph_8Km_DispDist), mode = "all", loops = FALSE, weights = E(Epcal_graph_8Km_DispDist)$weight),
  EgoSize = ego_size(Epcal_graph_8Km_DispDist, order = 3, nodes = V(Epcal_graph_8Km_DispDist))
)


## Peagg
#2Km
Peagg_graph_2Km_DispDist_metrics <- data.frame(
  PatchID = V(Peagg_graph_2Km_DispDist)$name,
  deg=degree(Peagg_graph_2Km_DispDist),
  strength = strength(Peagg_graph_2Km_DispDist, vids = V(Peagg_graph_2Km_DispDist), mode = "all", loops = FALSE, weights = E(Peagg_graph_2Km_DispDist)$weight),
  EgoSize = ego_size(Peagg_graph_2Km_DispDist, order = 3, nodes = V(Peagg_graph_2Km_DispDist))
)
#300m
Peagg_graph_300m_DispDist_metrics <- data.frame(
  PatchID = V(Peagg_graph_300m_DispDist)$name,
  deg=degree(Peagg_graph_300m_DispDist),
  strength = strength(Peagg_graph_300m_DispDist, vids = V(Peagg_graph_300m_DispDist), mode = "all", loops = FALSE, weights = E(Peagg_graph_300m_DispDist)$weight),
  EgoSize = ego_size(Peagg_graph_300m_DispDist, order = 3, nodes = V(Peagg_graph_300m_DispDist))
)
#10Km  
Peagg_graph_10Km_DispDist_metrics <- data.frame(
  PatchID = V(Peagg_graph_10Km_DispDist)$name,
  deg=degree(Peagg_graph_10Km_DispDist),
  strength = strength(Peagg_graph_10Km_DispDist, vids = V(Peagg_graph_10Km_DispDist), mode = "all", loops = FALSE, weights = E(Peagg_graph_10Km_DispDist)$weight),
  EgoSize = ego_size(Peagg_graph_10Km_DispDist, order = 3, nodes = V(Peagg_graph_10Km_DispDist))
)
#1Km
Peagg_graph_1Km_DispDist_metrics <- data.frame(
  PatchID = V(Peagg_graph_1Km_DispDist)$name,
  deg=degree(Peagg_graph_1Km_DispDist),
  strength = strength(Peagg_graph_1Km_DispDist, vids = V(Peagg_graph_1Km_DispDist), mode = "all", loops = FALSE, weights = E(Peagg_graph_1Km_DispDist)$weight),
  EgoSize = ego_size(Peagg_graph_1Km_DispDist, order = 3, nodes = V(Peagg_graph_1Km_DispDist))
)
#6Km
Peagg_graph_6Km_DispDist_metrics <- data.frame(
  PatchID = V(Peagg_graph_6Km_DispDist)$name,
  deg=degree(Peagg_graph_6Km_DispDist),
  strength = strength(Peagg_graph_6Km_DispDist, vids = V(Peagg_graph_6Km_DispDist), mode = "all", loops = FALSE, weights = E(Peagg_graph_6Km_DispDist)$weight),
  EgoSize = ego_size(Peagg_graph_6Km_DispDist, order = 3, nodes = V(Peagg_graph_6Km_DispDist))
)
#4Km
Peagg_graph_4Km_DispDist_metrics <- data.frame(
  PatchID = V(Peagg_graph_4Km_DispDist)$name,
  deg=degree(Peagg_graph_4Km_DispDist),
  strength = strength(Peagg_graph_4Km_DispDist, vids = V(Peagg_graph_4Km_DispDist), mode = "all", loops = FALSE, weights = E(Peagg_graph_4Km_DispDist)$weight),
  EgoSize = ego_size(Peagg_graph_4Km_DispDist, order = 3, nodes = V(Peagg_graph_4Km_DispDist))
)
#8Km
Peagg_graph_8Km_DispDist_metrics <- data.frame(
  PatchID = V(Peagg_graph_8Km_DispDist)$name,
  deg=degree(Peagg_graph_8Km_DispDist),
  strength = strength(Peagg_graph_8Km_DispDist, vids = V(Peagg_graph_8Km_DispDist), mode = "all", loops = FALSE, weights = E(Peagg_graph_8Km_DispDist)$weight),
  EgoSize = ego_size(Peagg_graph_8Km_DispDist, order = 3, nodes = V(Peagg_graph_8Km_DispDist))
)


## Perid
#2Km
Perid_graph_2Km_DispDist_metrics <- data.frame(
  PatchID = V(Perid_graph_2Km_DispDist)$name,
  deg=degree(Perid_graph_2Km_DispDist),
  strength = strength(Perid_graph_2Km_DispDist, vids = V(Perid_graph_2Km_DispDist), mode = "all", loops = FALSE, weights = E(Perid_graph_2Km_DispDist)$weight),
  EgoSize = ego_size(Perid_graph_2Km_DispDist, order = 3, nodes = V(Perid_graph_2Km_DispDist))
)
#300m
Perid_graph_300m_DispDist_metrics <- data.frame(
  PatchID = V(Perid_graph_300m_DispDist)$name,
  deg=degree(Perid_graph_300m_DispDist),
  strength = strength(Perid_graph_300m_DispDist, vids = V(Perid_graph_300m_DispDist), mode = "all", loops = FALSE, weights = E(Perid_graph_300m_DispDist)$weight),
  EgoSize = ego_size(Perid_graph_300m_DispDist, order = 3, nodes = V(Perid_graph_300m_DispDist))
)
#10Km  
Perid_graph_10Km_DispDist_metrics <- data.frame(
  PatchID = V(Perid_graph_10Km_DispDist)$name,
  deg=degree(Perid_graph_10Km_DispDist),
  strength = strength(Perid_graph_10Km_DispDist, vids = V(Perid_graph_10Km_DispDist), mode = "all", loops = FALSE, weights = E(Perid_graph_10Km_DispDist)$weight),
  EgoSize = ego_size(Perid_graph_10Km_DispDist, order = 3, nodes = V(Perid_graph_10Km_DispDist))
)
#1Km
Perid_graph_1Km_DispDist_metrics <- data.frame(
  PatchID = V(Perid_graph_1Km_DispDist)$name,
  deg=degree(Perid_graph_1Km_DispDist),
  strength = strength(Perid_graph_1Km_DispDist, vids = V(Perid_graph_1Km_DispDist), mode = "all", loops = FALSE, weights = E(Perid_graph_1Km_DispDist)$weight),
  EgoSize = ego_size(Perid_graph_1Km_DispDist, order = 3, nodes = V(Perid_graph_1Km_DispDist))
)
#6Km
Perid_graph_6Km_DispDist_metrics <- data.frame(
  PatchID = V(Perid_graph_6Km_DispDist)$name,
  deg=degree(Perid_graph_6Km_DispDist),
  strength = strength(Perid_graph_6Km_DispDist, vids = V(Perid_graph_6Km_DispDist), mode = "all", loops = FALSE, weights = E(Perid_graph_6Km_DispDist)$weight),
  EgoSize = ego_size(Perid_graph_6Km_DispDist, order = 3, nodes = V(Perid_graph_6Km_DispDist))
)
# #4Km
Perid_graph_4Km_DispDist_metrics <- data.frame(
  PatchID = V(Perid_graph_4Km_DispDist)$name,
  deg=degree(Perid_graph_4Km_DispDist),
  strength = strength(Perid_graph_4Km_DispDist, vids = V(Perid_graph_4Km_DispDist), mode = "all", loops = FALSE, weights = E(Perid_graph_4Km_DispDist)$weight),
  EgoSize = ego_size(Perid_graph_4Km_DispDist, order = 3, nodes = V(Perid_graph_4Km_DispDist))
)
# #8Km
Perid_graph_8Km_DispDist_metrics <- data.frame(
  PatchID = V(Perid_graph_8Km_DispDist)$name,
  deg=degree(Perid_graph_8Km_DispDist),
  strength = strength(Perid_graph_8Km_DispDist, vids = V(Perid_graph_8Km_DispDist), mode = "all", loops = FALSE, weights = E(Perid_graph_8Km_DispDist)$weight),
  EgoSize = ego_size(Perid_graph_8Km_DispDist, order = 3, nodes = V(Perid_graph_8Km_DispDist))
)



#### join the tables of patches and network metrics by attribute ####

Bovar_2Km_DispDist_topo_attributes <- Pres_abs_allpatches_Bo_var
Bovar_300m_DispDist_topo_attributes <- Pres_abs_allpatches_Bo_var
Bovar_10Km_DispDist_topo_attributes <- Pres_abs_allpatches_Bo_var
Bovar_1Km_DispDist_topo_attributes <- Pres_abs_allpatches_Bo_var
Bovar_6Km_DispDist_topo_attributes <- Pres_abs_allpatches_Bo_var
# Bovar_4Km_DispDist_topo_attributes <- Pres_abs_allpatches_Bo_var
Bovar_8Km_DispDist_topo_attributes <- Pres_abs_allpatches_Bo_var

Hyarb_2Km_DispDist_topo_attributes <- Pres_abs_allpatches_Hy_arb
Hyarb_300m_DispDist_topo_attributes <- Pres_abs_allpatches_Hy_arb
Hyarb_10Km_DispDist_topo_attributes <- Pres_abs_allpatches_Hy_arb
Hyarb_1Km_DispDist_topo_attributes <- Pres_abs_allpatches_Hy_arb
Hyarb_6Km_DispDist_topo_attributes <- Pres_abs_allpatches_Hy_arb
Hyarb_4Km_DispDist_topo_attributes <- Pres_abs_allpatches_Hy_arb
Hyarb_8Km_DispDist_topo_attributes <- Pres_abs_allpatches_Hy_arb

Alobs_2Km_DispDist_topo_attributes <- Pres_abs_allpatches_Al_obs
Alobs_300m_DispDist_topo_attributes <- Pres_abs_allpatches_Al_obs
Alobs_10Km_DispDist_topo_attributes <- Pres_abs_allpatches_Al_obs
Alobs_1Km_DispDist_topo_attributes <- Pres_abs_allpatches_Al_obs
Alobs_6Km_DispDist_topo_attributes <- Pres_abs_allpatches_Al_obs
Alobs_4Km_DispDist_topo_attributes <- Pres_abs_allpatches_Al_obs
Alobs_8Km_DispDist_topo_attributes <- Pres_abs_allpatches_Al_obs

Epcal_2Km_DispDist_topo_attributes <- Pres_abs_allpatches_Ep_cal
Epcal_300m_DispDist_topo_attributes <- Pres_abs_allpatches_Ep_cal
Epcal_10Km_DispDist_topo_attributes <- Pres_abs_allpatches_Ep_cal
Epcal_1Km_DispDist_topo_attributes <- Pres_abs_allpatches_Ep_cal
Epcal_6Km_DispDist_topo_attributes <- Pres_abs_allpatches_Ep_cal
Epcal_4Km_DispDist_topo_attributes <- Pres_abs_allpatches_Ep_cal
Epcal_8Km_DispDist_topo_attributes <- Pres_abs_allpatches_Ep_cal

Peagg_2Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_agg
Peagg_300m_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_agg
Peagg_10Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_agg
Peagg_1Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_agg
Peagg_6Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_agg
Peagg_4Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_agg
Peagg_8Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_agg

Perid_2Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_rid
Perid_300m_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_rid
Perid_10Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_rid
Perid_1Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_rid
Perid_6Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_rid
Perid_4Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_rid
Perid_8Km_DispDist_topo_attributes <- Pres_abs_allpatches_Pe_rid


### inner join= join by common variable names
# Bovar
Bovar_2Km_DispDist_topo_attributes <- merge(Bovar_2Km_DispDist_topo_attributes, Bovar_graph_2Km_DispDist_metrics, by = "PatchID")
Bovar_300m_DispDist_topo_attributes <- merge(Bovar_300m_DispDist_topo_attributes, Bovar_graph_300m_DispDist_metrics, by = "PatchID")
Bovar_10Km_DispDist_topo_attributes <- merge(Bovar_10Km_DispDist_topo_attributes, Bovar_graph_10Km_DispDist_metrics, by = "PatchID")
Bovar_1Km_DispDist_topo_attributes <- merge(Bovar_1Km_DispDist_topo_attributes, Bovar_graph_1Km_DispDist_metrics, by = "PatchID")
Bovar_6Km_DispDist_topo_attributes <- merge(Bovar_6Km_DispDist_topo_attributes, Bovar_graph_6Km_DispDist_metrics, by = "PatchID")
# Bovar_4Km_DispDist_topo_attributes <- merge(Bovar_4Km_DispDist_topo_attributes, Bovar_graph_4Km_DispDist_metrics, by = "PatchID") #default disp dist
Bovar_8Km_DispDist_topo_attributes <- merge(Bovar_8Km_DispDist_topo_attributes, Bovar_graph_8Km_DispDist_metrics, by = "PatchID")

#Hyarb
Hyarb_2Km_DispDist_topo_attributes <- merge(Hyarb_2Km_DispDist_topo_attributes, Hyarb_graph_2Km_DispDist_metrics, by = "PatchID")
Hyarb_300m_DispDist_topo_attributes <- merge(Hyarb_300m_DispDist_topo_attributes, Hyarb_graph_300m_DispDist_metrics, by = "PatchID")
Hyarb_10Km_DispDist_topo_attributes <- merge(Hyarb_10Km_DispDist_topo_attributes, Hyarb_graph_10Km_DispDist_metrics, by = "PatchID")
Hyarb_1Km_DispDist_topo_attributes <- merge(Hyarb_1Km_DispDist_topo_attributes, Hyarb_graph_1Km_DispDist_metrics, by = "PatchID")
Hyarb_6Km_DispDist_topo_attributes <- merge(Hyarb_6Km_DispDist_topo_attributes, Hyarb_graph_6Km_DispDist_metrics, by = "PatchID")
Hyarb_4Km_DispDist_topo_attributes <- merge(Hyarb_4Km_DispDist_topo_attributes, Hyarb_graph_4Km_DispDist_metrics, by = "PatchID") 
Hyarb_8Km_DispDist_topo_attributes <- merge(Hyarb_8Km_DispDist_topo_attributes, Hyarb_graph_8Km_DispDist_metrics, by = "PatchID")

#Alobs
Alobs_2Km_DispDist_topo_attributes <- merge(Alobs_2Km_DispDist_topo_attributes, Alobs_graph_2Km_DispDist_metrics, by = "PatchID")
Alobs_300m_DispDist_topo_attributes <- merge(Alobs_300m_DispDist_topo_attributes, Alobs_graph_300m_DispDist_metrics, by = "PatchID")
Alobs_10Km_DispDist_topo_attributes <- merge(Alobs_10Km_DispDist_topo_attributes, Alobs_graph_10Km_DispDist_metrics, by = "PatchID")
Alobs_1Km_DispDist_topo_attributes <- merge(Alobs_1Km_DispDist_topo_attributes, Alobs_graph_1Km_DispDist_metrics, by = "PatchID")
Alobs_6Km_DispDist_topo_attributes <- merge(Alobs_6Km_DispDist_topo_attributes, Alobs_graph_6Km_DispDist_metrics, by = "PatchID")
Alobs_4Km_DispDist_topo_attributes <- merge(Alobs_4Km_DispDist_topo_attributes, Alobs_graph_4Km_DispDist_metrics, by = "PatchID") 
Alobs_8Km_DispDist_topo_attributes <- merge(Alobs_8Km_DispDist_topo_attributes, Alobs_graph_8Km_DispDist_metrics, by = "PatchID")

#Epcal
Epcal_2Km_DispDist_topo_attributes <- merge(Epcal_2Km_DispDist_topo_attributes, Epcal_graph_2Km_DispDist_metrics, by = "PatchID")
Epcal_300m_DispDist_topo_attributes <- merge(Epcal_300m_DispDist_topo_attributes, Epcal_graph_300m_DispDist_metrics, by = "PatchID")
Epcal_10Km_DispDist_topo_attributes <- merge(Epcal_10Km_DispDist_topo_attributes, Epcal_graph_10Km_DispDist_metrics, by = "PatchID")
Epcal_1Km_DispDist_topo_attributes <- merge(Epcal_1Km_DispDist_topo_attributes, Epcal_graph_1Km_DispDist_metrics, by = "PatchID")
Epcal_6Km_DispDist_topo_attributes <- merge(Epcal_6Km_DispDist_topo_attributes, Epcal_graph_6Km_DispDist_metrics, by = "PatchID")
Epcal_4Km_DispDist_topo_attributes <- merge(Epcal_4Km_DispDist_topo_attributes, Epcal_graph_4Km_DispDist_metrics, by = "PatchID") 
Epcal_8Km_DispDist_topo_attributes <- merge(Epcal_8Km_DispDist_topo_attributes, Epcal_graph_8Km_DispDist_metrics, by = "PatchID")

#Peagg
Peagg_2Km_DispDist_topo_attributes <- merge(Peagg_2Km_DispDist_topo_attributes, Peagg_graph_2Km_DispDist_metrics, by = "PatchID")
Peagg_300m_DispDist_topo_attributes <- merge(Peagg_300m_DispDist_topo_attributes, Peagg_graph_300m_DispDist_metrics, by = "PatchID")
Peagg_10Km_DispDist_topo_attributes <- merge(Peagg_10Km_DispDist_topo_attributes, Peagg_graph_10Km_DispDist_metrics, by = "PatchID")
Peagg_1Km_DispDist_topo_attributes <- merge(Peagg_1Km_DispDist_topo_attributes, Peagg_graph_1Km_DispDist_metrics, by = "PatchID")
Peagg_6Km_DispDist_topo_attributes <- merge(Peagg_6Km_DispDist_topo_attributes, Peagg_graph_6Km_DispDist_metrics, by = "PatchID")
Peagg_4Km_DispDist_topo_attributes <- merge(Peagg_4Km_DispDist_topo_attributes, Peagg_graph_4Km_DispDist_metrics, by = "PatchID") 
Peagg_8Km_DispDist_topo_attributes <- merge(Peagg_8Km_DispDist_topo_attributes, Peagg_graph_8Km_DispDist_metrics, by = "PatchID")

#Perid
Perid_2Km_DispDist_topo_attributes <- merge(Perid_2Km_DispDist_topo_attributes, Perid_graph_2Km_DispDist_metrics, by = "PatchID")
Perid_300m_DispDist_topo_attributes <- merge(Perid_300m_DispDist_topo_attributes, Perid_graph_300m_DispDist_metrics, by = "PatchID")
Perid_10Km_DispDist_topo_attributes <- merge(Perid_10Km_DispDist_topo_attributes, Perid_graph_10Km_DispDist_metrics, by = "PatchID")
Perid_1Km_DispDist_topo_attributes <- merge(Perid_1Km_DispDist_topo_attributes, Perid_graph_1Km_DispDist_metrics, by = "PatchID")
Perid_6Km_DispDist_topo_attributes <- merge(Perid_6Km_DispDist_topo_attributes, Perid_graph_6Km_DispDist_metrics, by = "PatchID")
Perid_4Km_DispDist_topo_attributes <- merge(Perid_4Km_DispDist_topo_attributes, Perid_graph_4Km_DispDist_metrics, by = "PatchID") 
Perid_8Km_DispDist_topo_attributes <- merge(Perid_8Km_DispDist_topo_attributes, Perid_graph_8Km_DispDist_metrics, by = "PatchID")



### Add unweighted b_c to all new network df's
#### Measure betweenness centrality without weight included ####
#Re-import graphs, change name of 'weight' attribute
setwd('C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/Networks')
#Bovar
Bovar_graph_2Km_DispDist_2 <- read_graph("Bovar_2KmDispDist.graphml", format = "graphml")
Bovar_graph_300m_DispDist_2 <- read_graph("Bovar_300mDispDist.graphml", format = "graphml")
Bovar_graph_10Km_DispDist_2 <- read_graph("Bovar_10KmDispDist.graphml", format = "graphml")
Bovar_graph_1Km_DispDist_2 <- read_graph("Bovar_1KmDispDist.graphml", format = "graphml")
Bovar_graph_6Km_DispDist_2 <- read_graph("Bovar_6KmDispDist.graphml", format = "graphml")
# Bovar_graph_4Km_DispDist_2 <- read_graph("Bovar_4KmDispDist.graphml", format = "graphml")
Bovar_graph_8Km_DispDist_2 <- read_graph("Bovar_8KmDispDist.graphml", format = "graphml")

#Hyarb
Hyarb_graph_2Km_DispDist_2 <- read_graph("Hyarb_2KmDispDist.graphml", format = "graphml")
Hyarb_graph_300m_DispDist_2 <- read_graph("Hyarb_300mDispDist.graphml", format = "graphml")
Hyarb_graph_10Km_DispDist_2 <- read_graph("Hyarb_10KmDispDist.graphml", format = "graphml")
Hyarb_graph_1Km_DispDist_2 <- read_graph("Hyarb_1KmDispDist.graphml", format = "graphml")
Hyarb_graph_6Km_DispDist_2 <- read_graph("Hyarb_6KmDispDist.graphml", format = "graphml")
Hyarb_graph_4Km_DispDist_2 <- read_graph("Hyarb_4KmDispDist.graphml", format = "graphml")
Hyarb_graph_8Km_DispDist_2 <- read_graph("Hyarb_8KmDispDist.graphml", format = "graphml")

#Alobs
Alobs_graph_2Km_DispDist_2 <- read_graph("Alobs_2KmDispDist.graphml", format = "graphml")
Alobs_graph_300m_DispDist_2 <- read_graph("Alobs_300mDispDist.graphml", format = "graphml")
Alobs_graph_10Km_DispDist_2 <- read_graph("Alobs_10KmDispDist.graphml", format = "graphml")
Alobs_graph_1Km_DispDist_2 <- read_graph("Alobs_1KmDispDist.graphml", format = "graphml")
Alobs_graph_6Km_DispDist_2 <- read_graph("Alobs_6KmDispDist.graphml", format = "graphml")
Alobs_graph_4Km_DispDist_2 <- read_graph("Alobs_4KmDispDist.graphml", format = "graphml")
Alobs_graph_8Km_DispDist_2 <- read_graph("Alobs_8KmDispDist.graphml", format = "graphml")

#Epcal
Epcal_graph_2Km_DispDist_2 <- read_graph("Epcal_2KmDispDist.graphml", format = "graphml")
Epcal_graph_300m_DispDist_2 <- read_graph("Epcal_300mDispDist.graphml", format = "graphml")
Epcal_graph_10Km_DispDist_2 <- read_graph("Epcal_10KmDispDist.graphml", format = "graphml")
Epcal_graph_1Km_DispDist_2 <- read_graph("Epcal_1KmDispDist.graphml", format = "graphml")
Epcal_graph_6Km_DispDist_2 <- read_graph("Epcal_6KmDispDist.graphml", format = "graphml")
Epcal_graph_4Km_DispDist_2 <- read_graph("Epcal_4KmDispDist.graphml", format = "graphml")
Epcal_graph_8Km_DispDist_2 <- read_graph("Epcal_8KmDispDist.graphml", format = "graphml")

#Peagg
Peagg_graph_2Km_DispDist_2 <- read_graph("Peagg_2KmDispDist.graphml", format = "graphml")
Peagg_graph_300m_DispDist_2 <- read_graph("Peagg_300mDispDist.graphml", format = "graphml")
Peagg_graph_10Km_DispDist_2 <- read_graph("Peagg_10KmDispDist.graphml", format = "graphml")
Peagg_graph_1Km_DispDist_2 <- read_graph("Peagg_1KmDispDist.graphml", format = "graphml")
Peagg_graph_6Km_DispDist_2 <- read_graph("Peagg_6KmDispDist.graphml", format = "graphml")
Peagg_graph_4Km_DispDist_2 <- read_graph("Peagg_4KmDispDist.graphml", format = "graphml")
Peagg_graph_8Km_DispDist_2 <- read_graph("Peagg_8KmDispDist.graphml", format = "graphml")

#Perid
Perid_graph_2Km_DispDist_2 <- read_graph("Perid_2KmDispDist.graphml", format = "graphml")
Perid_graph_300m_DispDist_2 <- read_graph("Perid_300mDispDist.graphml", format = "graphml")
Perid_graph_10Km_DispDist_2 <- read_graph("Perid_10KmDispDist.graphml", format = "graphml")
Perid_graph_1Km_DispDist_2 <- read_graph("Perid_1KmDispDist.graphml", format = "graphml")
Perid_graph_6Km_DispDist_2 <- read_graph("Perid_6KmDispDist.graphml", format = "graphml")
Perid_graph_4Km_DispDist_2 <- read_graph("Perid_4KmDispDist.graphml", format = "graphml")
Perid_graph_8Km_DispDist_2 <- read_graph("Perid_8KmDispDist.graphml", format = "graphml")


#rename weight attribute
#Actually copying into a new attr. and deleting original attribute
#Bovar
E(Bovar_graph_2Km_DispDist_2)$Alter_wght <- E(Bovar_graph_2Km_DispDist_2)$weight
Bovar_graph_2Km_DispDist_2 <- remove.edge.attribute(Bovar_graph_2Km_DispDist_2, "weight")
E(Bovar_graph_300m_DispDist_2)$Alter_wght <- E(Bovar_graph_300m_DispDist_2)$weight
Bovar_graph_300m_DispDist_2 <- remove.edge.attribute(Bovar_graph_300m_DispDist_2, "weight")
E(Bovar_graph_10Km_DispDist_2)$Alter_wght <- E(Bovar_graph_10Km_DispDist_2)$weight
Bovar_graph_10Km_DispDist_2 <- remove.edge.attribute(Bovar_graph_10Km_DispDist_2, "weight")
E(Bovar_graph_1Km_DispDist_2)$Alter_wght <- E(Bovar_graph_1Km_DispDist_2)$weight
Bovar_graph_1Km_DispDist_2 <- remove.edge.attribute(Bovar_graph_1Km_DispDist_2, "weight")
E(Bovar_graph_6Km_DispDist_2)$Alter_wght <- E(Bovar_graph_6Km_DispDist_2)$weight
Bovar_graph_6Km_DispDist_2 <- remove.edge.attribute(Bovar_graph_6Km_DispDist_2, "weight")
E(Bovar_graph_8Km_DispDist_2)$Alter_wght <- E(Bovar_graph_8Km_DispDist_2)$weight
Bovar_graph_8Km_DispDist_2 <- remove.edge.attribute(Bovar_graph_8Km_DispDist_2, "weight")
#Hyarb
E(Hyarb_graph_2Km_DispDist_2)$Alter_wght <- E(Hyarb_graph_2Km_DispDist_2)$weight
Hyarb_graph_2Km_DispDist_2 <- remove.edge.attribute(Hyarb_graph_2Km_DispDist_2, "weight")
E(Hyarb_graph_300m_DispDist_2)$Alter_wght <- E(Hyarb_graph_300m_DispDist_2)$weight
Hyarb_graph_300m_DispDist_2 <- remove.edge.attribute(Hyarb_graph_300m_DispDist_2, "weight")
E(Hyarb_graph_10Km_DispDist_2)$Alter_wght <- E(Hyarb_graph_10Km_DispDist_2)$weight
Hyarb_graph_10Km_DispDist_2 <- remove.edge.attribute(Hyarb_graph_10Km_DispDist_2, "weight")
E(Hyarb_graph_1Km_DispDist_2)$Alter_wght <- E(Hyarb_graph_1Km_DispDist_2)$weight
Hyarb_graph_1Km_DispDist_2 <- remove.edge.attribute(Hyarb_graph_1Km_DispDist_2, "weight")
E(Hyarb_graph_4Km_DispDist_2)$Alter_wght <- E(Hyarb_graph_4Km_DispDist_2)$weight
Hyarb_graph_4Km_DispDist_2 <- remove.edge.attribute(Hyarb_graph_4Km_DispDist_2, "weight")
E(Hyarb_graph_6Km_DispDist_2)$Alter_wght <- E(Hyarb_graph_6Km_DispDist_2)$weight
Hyarb_graph_6Km_DispDist_2 <- remove.edge.attribute(Hyarb_graph_6Km_DispDist_2, "weight")
E(Hyarb_graph_8Km_DispDist_2)$Alter_wght <- E(Hyarb_graph_8Km_DispDist_2)$weight
Hyarb_graph_8Km_DispDist_2 <- remove.edge.attribute(Hyarb_graph_8Km_DispDist_2, "weight")
# Alobs
E(Alobs_graph_2Km_DispDist_2)$Alter_wght <- E(Alobs_graph_2Km_DispDist_2)$weight
Alobs_graph_2Km_DispDist_2 <- remove.edge.attribute(Alobs_graph_2Km_DispDist_2, "weight")
E(Alobs_graph_300m_DispDist_2)$Alter_wght <- E(Alobs_graph_300m_DispDist_2)$weight
Alobs_graph_300m_DispDist_2 <- remove.edge.attribute(Alobs_graph_300m_DispDist_2, "weight")
E(Alobs_graph_10Km_DispDist_2)$Alter_wght <- E(Alobs_graph_10Km_DispDist_2)$weight
Alobs_graph_10Km_DispDist_2 <- remove.edge.attribute(Alobs_graph_10Km_DispDist_2, "weight")
E(Alobs_graph_1Km_DispDist_2)$Alter_wght <- E(Alobs_graph_1Km_DispDist_2)$weight
Alobs_graph_1Km_DispDist_2 <- remove.edge.attribute(Alobs_graph_1Km_DispDist_2, "weight")
E(Alobs_graph_4Km_DispDist_2)$Alter_wght <- E(Alobs_graph_4Km_DispDist_2)$weight
Alobs_graph_4Km_DispDist_2 <- remove.edge.attribute(Alobs_graph_4Km_DispDist_2, "weight")
E(Alobs_graph_6Km_DispDist_2)$Alter_wght <- E(Alobs_graph_6Km_DispDist_2)$weight
Alobs_graph_6Km_DispDist_2 <- remove.edge.attribute(Alobs_graph_6Km_DispDist_2, "weight")
E(Alobs_graph_8Km_DispDist_2)$Alter_wght <- E(Alobs_graph_8Km_DispDist_2)$weight
Alobs_graph_8Km_DispDist_2 <- remove.edge.attribute(Alobs_graph_8Km_DispDist_2, "weight")
# Epcal
E(Epcal_graph_2Km_DispDist_2)$Alter_wght <- E(Epcal_graph_2Km_DispDist_2)$weight
Epcal_graph_2Km_DispDist_2 <- remove.edge.attribute(Epcal_graph_2Km_DispDist_2, "weight")
E(Epcal_graph_300m_DispDist_2)$Alter_wght <- E(Epcal_graph_300m_DispDist_2)$weight
Epcal_graph_300m_DispDist_2 <- remove.edge.attribute(Epcal_graph_300m_DispDist_2, "weight")
E(Epcal_graph_10Km_DispDist_2)$Alter_wght <- E(Epcal_graph_10Km_DispDist_2)$weight
Epcal_graph_10Km_DispDist_2 <- remove.edge.attribute(Epcal_graph_10Km_DispDist_2, "weight")
E(Epcal_graph_1Km_DispDist_2)$Alter_wght <- E(Epcal_graph_1Km_DispDist_2)$weight
Epcal_graph_1Km_DispDist_2 <- remove.edge.attribute(Epcal_graph_1Km_DispDist_2, "weight")
E(Epcal_graph_4Km_DispDist_2)$Alter_wght <- E(Epcal_graph_4Km_DispDist_2)$weight
Epcal_graph_4Km_DispDist_2 <- remove.edge.attribute(Epcal_graph_4Km_DispDist_2, "weight")
E(Epcal_graph_6Km_DispDist_2)$Alter_wght <- E(Epcal_graph_6Km_DispDist_2)$weight
Epcal_graph_6Km_DispDist_2 <- remove.edge.attribute(Epcal_graph_6Km_DispDist_2, "weight")
E(Epcal_graph_8Km_DispDist_2)$Alter_wght <- E(Epcal_graph_8Km_DispDist_2)$weight
Epcal_graph_8Km_DispDist_2 <- remove.edge.attribute(Epcal_graph_8Km_DispDist_2, "weight")
# Peagg
E(Peagg_graph_2Km_DispDist_2)$Alter_wght <- E(Peagg_graph_2Km_DispDist_2)$weight
Peagg_graph_2Km_DispDist_2 <- remove.edge.attribute(Peagg_graph_2Km_DispDist_2, "weight")
E(Peagg_graph_300m_DispDist_2)$Alter_wght <- E(Peagg_graph_300m_DispDist_2)$weight
Peagg_graph_300m_DispDist_2 <- remove.edge.attribute(Peagg_graph_300m_DispDist_2, "weight")
E(Peagg_graph_10Km_DispDist_2)$Alter_wght <- E(Peagg_graph_10Km_DispDist_2)$weight
Peagg_graph_10Km_DispDist_2 <- remove.edge.attribute(Peagg_graph_10Km_DispDist_2, "weight")
E(Peagg_graph_1Km_DispDist_2)$Alter_wght <- E(Peagg_graph_1Km_DispDist_2)$weight
Peagg_graph_1Km_DispDist_2 <- remove.edge.attribute(Peagg_graph_1Km_DispDist_2, "weight")
E(Peagg_graph_4Km_DispDist_2)$Alter_wght <- E(Peagg_graph_4Km_DispDist_2)$weight
Peagg_graph_4Km_DispDist_2 <- remove.edge.attribute(Peagg_graph_4Km_DispDist_2, "weight")
E(Peagg_graph_6Km_DispDist_2)$Alter_wght <- E(Peagg_graph_6Km_DispDist_2)$weight
Peagg_graph_6Km_DispDist_2 <- remove.edge.attribute(Peagg_graph_6Km_DispDist_2, "weight")
E(Peagg_graph_8Km_DispDist_2)$Alter_wght <- E(Peagg_graph_8Km_DispDist_2)$weight
Peagg_graph_8Km_DispDist_2 <- remove.edge.attribute(Peagg_graph_8Km_DispDist_2, "weight")
# Perid
E(Perid_graph_2Km_DispDist_2)$Alter_wght <- E(Perid_graph_2Km_DispDist_2)$weight
Perid_graph_2Km_DispDist_2 <- remove.edge.attribute(Perid_graph_2Km_DispDist_2, "weight")
E(Perid_graph_300m_DispDist_2)$Alter_wght <- E(Perid_graph_300m_DispDist_2)$weight
Perid_graph_300m_DispDist_2 <- remove.edge.attribute(Perid_graph_300m_DispDist_2, "weight")
E(Perid_graph_10Km_DispDist_2)$Alter_wght <- E(Perid_graph_10Km_DispDist_2)$weight
Perid_graph_10Km_DispDist_2 <- remove.edge.attribute(Perid_graph_10Km_DispDist_2, "weight")
E(Perid_graph_1Km_DispDist_2)$Alter_wght <- E(Perid_graph_1Km_DispDist_2)$weight
Perid_graph_1Km_DispDist_2 <- remove.edge.attribute(Perid_graph_1Km_DispDist_2, "weight")
E(Perid_graph_4Km_DispDist_2)$Alter_wght <- E(Perid_graph_4Km_DispDist_2)$weight
Perid_graph_4Km_DispDist_2 <- remove.edge.attribute(Perid_graph_4Km_DispDist_2, "weight")
E(Perid_graph_6Km_DispDist_2)$Alter_wght <- E(Perid_graph_6Km_DispDist_2)$weight
Perid_graph_6Km_DispDist_2 <- remove.edge.attribute(Perid_graph_6Km_DispDist_2, "weight")
E(Perid_graph_8Km_DispDist_2)$Alter_wght <- E(Perid_graph_8Km_DispDist_2)$weight
Perid_graph_8Km_DispDist_2 <- remove.edge.attribute(Perid_graph_8Km_DispDist_2, "weight")


#Calculate explicitely unweighted b_c
#Set dataframe w/only PatchID & unw_b_c

# Bovar
Bovar_2Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Bovar_graph_2Km_DispDist_2)$name,
  b_c=betweenness(Bovar_graph_2Km_DispDist_2, weights = NULL))
Bovar_300m_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Bovar_graph_300m_DispDist_2)$name,
  b_c=betweenness(Bovar_graph_300m_DispDist_2, weights = NULL))
Bovar_10Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Bovar_graph_10Km_DispDist_2)$name,
  b_c=betweenness(Bovar_graph_10Km_DispDist_2, weights = NULL))
Bovar_1Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Bovar_graph_1Km_DispDist_2)$name,
  b_c=betweenness(Bovar_graph_1Km_DispDist_2, weights = NULL))
Bovar_6Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Bovar_graph_6Km_DispDist_2)$name,
  b_c=betweenness(Bovar_graph_6Km_DispDist_2, weights = NULL))
Bovar_8Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Bovar_graph_8Km_DispDist_2)$name,
  b_c=betweenness(Bovar_graph_8Km_DispDist_2, weights = NULL))

# Hyarb
Hyarb_2Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Hyarb_graph_2Km_DispDist_2)$name,
  b_c=betweenness(Hyarb_graph_2Km_DispDist_2, weights = NULL))
Hyarb_300m_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Hyarb_graph_300m_DispDist_2)$name,
  b_c=betweenness(Hyarb_graph_300m_DispDist_2, weights = NULL))
Hyarb_10Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Hyarb_graph_10Km_DispDist_2)$name,
  b_c=betweenness(Hyarb_graph_10Km_DispDist_2, weights = NULL))
Hyarb_1Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Hyarb_graph_1Km_DispDist_2)$name,
  b_c=betweenness(Hyarb_graph_1Km_DispDist_2, weights = NULL))
Hyarb_6Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Hyarb_graph_6Km_DispDist_2)$name,
  b_c=betweenness(Hyarb_graph_6Km_DispDist_2, weights = NULL))
Hyarb_4Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Hyarb_graph_4Km_DispDist_2)$name,
  b_c=betweenness(Hyarb_graph_4Km_DispDist_2, weights = NULL))
Hyarb_8Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Hyarb_graph_8Km_DispDist_2)$name,
  b_c=betweenness(Hyarb_graph_8Km_DispDist_2, weights = NULL))

# Alobs
Alobs_2Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Alobs_graph_2Km_DispDist_2)$name,
  b_c=betweenness(Alobs_graph_2Km_DispDist_2, weights = NULL))
Alobs_300m_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Alobs_graph_300m_DispDist_2)$name,
  b_c=betweenness(Alobs_graph_300m_DispDist_2, weights = NULL))
Alobs_10Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Alobs_graph_10Km_DispDist_2)$name,
  b_c=betweenness(Alobs_graph_10Km_DispDist_2, weights = NULL))
Alobs_1Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Alobs_graph_1Km_DispDist_2)$name,
  b_c=betweenness(Alobs_graph_1Km_DispDist_2, weights = NULL))
Alobs_6Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Alobs_graph_6Km_DispDist_2)$name,
  b_c=betweenness(Alobs_graph_6Km_DispDist_2, weights = NULL))
Alobs_4Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Alobs_graph_4Km_DispDist_2)$name,
  b_c=betweenness(Alobs_graph_4Km_DispDist_2, weights = NULL))
Alobs_8Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Alobs_graph_8Km_DispDist_2)$name,
  b_c=betweenness(Alobs_graph_8Km_DispDist_2, weights = NULL))

# Epcal
Epcal_2Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Epcal_graph_2Km_DispDist_2)$name,
  b_c=betweenness(Epcal_graph_2Km_DispDist_2, weights = NULL))
Epcal_300m_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Epcal_graph_300m_DispDist_2)$name,
  b_c=betweenness(Epcal_graph_300m_DispDist_2, weights = NULL))
Epcal_10Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Epcal_graph_10Km_DispDist_2)$name,
  b_c=betweenness(Epcal_graph_10Km_DispDist_2, weights = NULL))
Epcal_1Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Epcal_graph_1Km_DispDist_2)$name,
  b_c=betweenness(Epcal_graph_1Km_DispDist_2, weights = NULL))
Epcal_6Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Epcal_graph_6Km_DispDist_2)$name,
  b_c=betweenness(Epcal_graph_6Km_DispDist_2, weights = NULL))
Epcal_4Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Epcal_graph_4Km_DispDist_2)$name,
  b_c=betweenness(Epcal_graph_4Km_DispDist_2, weights = NULL))
Epcal_8Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Epcal_graph_8Km_DispDist_2)$name,
  b_c=betweenness(Epcal_graph_8Km_DispDist_2, weights = NULL))

# Peagg
Peagg_2Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Peagg_graph_2Km_DispDist_2)$name,
  b_c=betweenness(Peagg_graph_2Km_DispDist_2, weights = NULL))
Peagg_300m_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Peagg_graph_300m_DispDist_2)$name,
  b_c=betweenness(Peagg_graph_300m_DispDist_2, weights = NULL))
Peagg_10Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Peagg_graph_10Km_DispDist_2)$name,
  b_c=betweenness(Peagg_graph_10Km_DispDist_2, weights = NULL))
Peagg_1Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Peagg_graph_1Km_DispDist_2)$name,
  b_c=betweenness(Peagg_graph_1Km_DispDist_2, weights = NULL))
Peagg_6Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Peagg_graph_6Km_DispDist_2)$name,
  b_c=betweenness(Peagg_graph_6Km_DispDist_2, weights = NULL))
Peagg_4Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Peagg_graph_4Km_DispDist_2)$name,
  b_c=betweenness(Peagg_graph_4Km_DispDist_2, weights = NULL))
Peagg_8Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Peagg_graph_8Km_DispDist_2)$name,
  b_c=betweenness(Peagg_graph_8Km_DispDist_2, weights = NULL))

# Perid
Perid_2Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Perid_graph_2Km_DispDist_2)$name,
  b_c=betweenness(Perid_graph_2Km_DispDist_2, weights = NULL))
Perid_300m_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Perid_graph_300m_DispDist_2)$name,
  b_c=betweenness(Perid_graph_300m_DispDist_2, weights = NULL))
Perid_10Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Perid_graph_10Km_DispDist_2)$name,
  b_c=betweenness(Perid_graph_10Km_DispDist_2, weights = NULL))
Perid_1Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Perid_graph_1Km_DispDist_2)$name,
  b_c=betweenness(Perid_graph_1Km_DispDist_2, weights = NULL))
Perid_6Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Perid_graph_6Km_DispDist_2)$name,
  b_c=betweenness(Perid_graph_6Km_DispDist_2, weights = NULL))
Perid_4Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Perid_graph_4Km_DispDist_2)$name,
  b_c=betweenness(Perid_graph_4Km_DispDist_2, weights = NULL))
Perid_8Km_DispDist_unweighted_b_c <- data.frame(
  PatchID = V(Perid_graph_8Km_DispDist_2)$name,
  b_c=betweenness(Perid_graph_8Km_DispDist_2, weights = NULL))


### rename unw_b_c
# Bovar
Bovar_2Km_DispDist_unweighted_b_c <- rename(Bovar_2Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Bovar_300m_DispDist_unweighted_b_c <- rename(Bovar_300m_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Bovar_10Km_DispDist_unweighted_b_c <- rename(Bovar_10Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Bovar_1Km_DispDist_unweighted_b_c <- rename(Bovar_1Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Bovar_6Km_DispDist_unweighted_b_c <- rename(Bovar_6Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
# Bovar_4Km_DispDist_unweighted_b_c <- rename(Bovar_4Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Bovar_8Km_DispDist_unweighted_b_c <- rename(Bovar_8Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")

# Hyarb
Hyarb_2Km_DispDist_unweighted_b_c <- rename(Hyarb_2Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Hyarb_300m_DispDist_unweighted_b_c <- rename(Hyarb_300m_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Hyarb_10Km_DispDist_unweighted_b_c <- rename(Hyarb_10Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Hyarb_1Km_DispDist_unweighted_b_c <- rename(Hyarb_1Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Hyarb_6Km_DispDist_unweighted_b_c <- rename(Hyarb_6Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Hyarb_4Km_DispDist_unweighted_b_c <- rename(Hyarb_4Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Hyarb_8Km_DispDist_unweighted_b_c <- rename(Hyarb_8Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")

# Alobs
Alobs_2Km_DispDist_unweighted_b_c <- rename(Alobs_2Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Alobs_300m_DispDist_unweighted_b_c <- rename(Alobs_300m_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Alobs_10Km_DispDist_unweighted_b_c <- rename(Alobs_10Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Alobs_1Km_DispDist_unweighted_b_c <- rename(Alobs_1Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Alobs_6Km_DispDist_unweighted_b_c <- rename(Alobs_6Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Alobs_4Km_DispDist_unweighted_b_c <- rename(Alobs_4Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Alobs_8Km_DispDist_unweighted_b_c <- rename(Alobs_8Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")

# Epcal
Epcal_2Km_DispDist_unweighted_b_c <- rename(Epcal_2Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Epcal_300m_DispDist_unweighted_b_c <- rename(Epcal_300m_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Epcal_10Km_DispDist_unweighted_b_c <- rename(Epcal_10Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Epcal_1Km_DispDist_unweighted_b_c <- rename(Epcal_1Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Epcal_6Km_DispDist_unweighted_b_c <- rename(Epcal_6Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Epcal_4Km_DispDist_unweighted_b_c <- rename(Epcal_4Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Epcal_8Km_DispDist_unweighted_b_c <- rename(Epcal_8Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")

# Peagg
Peagg_2Km_DispDist_unweighted_b_c <- rename(Peagg_2Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Peagg_300m_DispDist_unweighted_b_c <- rename(Peagg_300m_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Peagg_10Km_DispDist_unweighted_b_c <- rename(Peagg_10Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Peagg_1Km_DispDist_unweighted_b_c <- rename(Peagg_1Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Peagg_6Km_DispDist_unweighted_b_c <- rename(Peagg_6Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Peagg_4Km_DispDist_unweighted_b_c <- rename(Peagg_4Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Peagg_8Km_DispDist_unweighted_b_c <- rename(Peagg_8Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")

# Perid
Perid_2Km_DispDist_unweighted_b_c <- rename(Perid_2Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Perid_300m_DispDist_unweighted_b_c <- rename(Perid_300m_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Perid_10Km_DispDist_unweighted_b_c <- rename(Perid_10Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Perid_1Km_DispDist_unweighted_b_c <- rename(Perid_1Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Perid_6Km_DispDist_unweighted_b_c <- rename(Perid_6Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Perid_4Km_DispDist_unweighted_b_c <- rename(Perid_4Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")
Perid_8Km_DispDist_unweighted_b_c <- rename(Perid_8Km_DispDist_unweighted_b_c,  "unw_b_c" = "b_c")


### join the tables by attribute
# Bovar
Bovar_2Km_DispDist_topo_attributes <- merge(x = Bovar_2Km_DispDist_topo_attributes, y = Bovar_2Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Bovar_2Km_DispDist_topo_attributes)
Bovar_300m_DispDist_topo_attributes <- merge(x = Bovar_300m_DispDist_topo_attributes, y = Bovar_300m_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Bovar_300m_DispDist_topo_attributes)
Bovar_10Km_DispDist_topo_attributes <- merge(x = Bovar_10Km_DispDist_topo_attributes, y = Bovar_10Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Bovar_10Km_DispDist_topo_attributes)
Bovar_1Km_DispDist_topo_attributes <- merge(x = Bovar_1Km_DispDist_topo_attributes, y = Bovar_1Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Bovar_1Km_DispDist_topo_attributes)
Bovar_6Km_DispDist_topo_attributes <- merge(x = Bovar_6Km_DispDist_topo_attributes, y = Bovar_6Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Bovar_6Km_DispDist_topo_attributes)
# Bovar_4Km_DispDist_topo_attributes <- merge(x = Bovar_4Km_DispDist_topo_attributes, y = Bovar_4Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
# names(Bovar_4Km_DispDist_topo_attributes)
Bovar_8Km_DispDist_topo_attributes <- merge(x = Bovar_8Km_DispDist_topo_attributes, y = Bovar_8Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Bovar_8Km_DispDist_topo_attributes)

# Hyarb
Hyarb_2Km_DispDist_topo_attributes <- merge(x = Hyarb_2Km_DispDist_topo_attributes, y = Hyarb_2Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Hyarb_2Km_DispDist_topo_attributes)
Hyarb_300m_DispDist_topo_attributes <- merge(x = Hyarb_300m_DispDist_topo_attributes, y = Hyarb_300m_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Hyarb_300m_DispDist_topo_attributes)
Hyarb_10Km_DispDist_topo_attributes <- merge(x = Hyarb_10Km_DispDist_topo_attributes, y = Hyarb_10Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Hyarb_10Km_DispDist_topo_attributes)
Hyarb_1Km_DispDist_topo_attributes <- merge(x = Hyarb_1Km_DispDist_topo_attributes, y = Hyarb_1Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Hyarb_6Km_DispDist_topo_attributes <- merge(x = Hyarb_6Km_DispDist_topo_attributes, y = Hyarb_6Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Hyarb_4Km_DispDist_topo_attributes <- merge(x = Hyarb_4Km_DispDist_topo_attributes, y = Hyarb_4Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Hyarb_8Km_DispDist_topo_attributes <- merge(x = Hyarb_8Km_DispDist_topo_attributes, y = Hyarb_8Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)

# Alobs
Alobs_2Km_DispDist_topo_attributes <- merge(x = Alobs_2Km_DispDist_topo_attributes, y = Alobs_2Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Alobs_2Km_DispDist_topo_attributes)
Alobs_300m_DispDist_topo_attributes <- merge(x = Alobs_300m_DispDist_topo_attributes, y = Alobs_300m_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Alobs_300m_DispDist_topo_attributes)
Alobs_10Km_DispDist_topo_attributes <- merge(x = Alobs_10Km_DispDist_topo_attributes, y = Alobs_10Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Alobs_10Km_DispDist_topo_attributes)
Alobs_1Km_DispDist_topo_attributes <- merge(x = Alobs_1Km_DispDist_topo_attributes, y = Alobs_1Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Alobs_6Km_DispDist_topo_attributes <- merge(x = Alobs_6Km_DispDist_topo_attributes, y = Alobs_6Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Alobs_4Km_DispDist_topo_attributes <- merge(x = Alobs_4Km_DispDist_topo_attributes, y = Alobs_4Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Alobs_8Km_DispDist_topo_attributes <- merge(x = Alobs_8Km_DispDist_topo_attributes, y = Alobs_8Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)

# Epcal
Epcal_2Km_DispDist_topo_attributes <- merge(x = Epcal_2Km_DispDist_topo_attributes, y = Epcal_2Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Epcal_2Km_DispDist_topo_attributes)
Epcal_300m_DispDist_topo_attributes <- merge(x = Epcal_300m_DispDist_topo_attributes, y = Epcal_300m_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Epcal_300m_DispDist_topo_attributes)
Epcal_10Km_DispDist_topo_attributes <- merge(x = Epcal_10Km_DispDist_topo_attributes, y = Epcal_10Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Epcal_10Km_DispDist_topo_attributes)
Epcal_1Km_DispDist_topo_attributes <- merge(x = Epcal_1Km_DispDist_topo_attributes, y = Epcal_1Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Epcal_6Km_DispDist_topo_attributes <- merge(x = Epcal_6Km_DispDist_topo_attributes, y = Epcal_6Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Epcal_4Km_DispDist_topo_attributes <- merge(x = Epcal_4Km_DispDist_topo_attributes, y = Epcal_4Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Epcal_8Km_DispDist_topo_attributes <- merge(x = Epcal_8Km_DispDist_topo_attributes, y = Epcal_8Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)

# Peagg
Peagg_2Km_DispDist_topo_attributes <- merge(x = Peagg_2Km_DispDist_topo_attributes, y = Peagg_2Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Peagg_2Km_DispDist_topo_attributes)
Peagg_300m_DispDist_topo_attributes <- merge(x = Peagg_300m_DispDist_topo_attributes, y = Peagg_300m_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Peagg_300m_DispDist_topo_attributes)
Peagg_10Km_DispDist_topo_attributes <- merge(x = Peagg_10Km_DispDist_topo_attributes, y = Peagg_10Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Peagg_10Km_DispDist_topo_attributes)
Peagg_1Km_DispDist_topo_attributes <- merge(x = Peagg_1Km_DispDist_topo_attributes, y = Peagg_1Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Peagg_6Km_DispDist_topo_attributes <- merge(x = Peagg_6Km_DispDist_topo_attributes, y = Peagg_6Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Peagg_4Km_DispDist_topo_attributes <- merge(x = Peagg_4Km_DispDist_topo_attributes, y = Peagg_4Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Peagg_8Km_DispDist_topo_attributes <- merge(x = Peagg_8Km_DispDist_topo_attributes, y = Peagg_8Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)

# Perid
Perid_2Km_DispDist_topo_attributes <- merge(x = Perid_2Km_DispDist_topo_attributes, y = Perid_2Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Perid_2Km_DispDist_topo_attributes)
Perid_300m_DispDist_topo_attributes <- merge(x = Perid_300m_DispDist_topo_attributes, y = Perid_300m_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Perid_300m_DispDist_topo_attributes)
Perid_10Km_DispDist_topo_attributes <- merge(x = Perid_10Km_DispDist_topo_attributes, y = Perid_10Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(Perid_10Km_DispDist_topo_attributes)
Perid_1Km_DispDist_topo_attributes <- merge(x = Perid_1Km_DispDist_topo_attributes, y = Perid_1Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Perid_6Km_DispDist_topo_attributes <- merge(x = Perid_6Km_DispDist_topo_attributes, y = Perid_6Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Perid_4Km_DispDist_topo_attributes <- merge(x = Perid_4Km_DispDist_topo_attributes, y = Perid_4Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)
Perid_8Km_DispDist_topo_attributes <- merge(x = Perid_8Km_DispDist_topo_attributes, y = Perid_8Km_DispDist_unweighted_b_c, by = "PatchID", all.x = TRUE)



#### Add Habitat Availability predictor (HabAv)for newly developed networks ####
# Bovar
Bovar_2KmDispDist_habAv <- scan("Bovar_2KmDispDist_habAv.txt")
Bovar_2Km_DispDist_topo_attributes$habAv <- Bovar_2KmDispDist_habAv
names(Bovar_2Km_DispDist_topo_attributes)
Bovar_300mDispDist_habAv <- scan("Bovar_300mDispDist_habAv.txt")
Bovar_300m_DispDist_topo_attributes$habAv <- Bovar_300mDispDist_habAv
names(Bovar_300m_DispDist_topo_attributes)
Bovar_10KmDispDist_habAv <- scan("Bovar_10KmDispDist_habAv.txt")
Bovar_10Km_DispDist_topo_attributes$habAv <- Bovar_10KmDispDist_habAv
names(Bovar_10Km_DispDist_topo_attributes)
Bovar_1KmDispDist_habAv <- scan("Bovar_1KmDispDist_habAv.txt")
Bovar_1Km_DispDist_topo_attributes$habAv <- Bovar_1KmDispDist_habAv
Bovar_6KmDispDist_habAv <- scan("Bovar_6KmDispDist_habAv.txt")
Bovar_6Km_DispDist_topo_attributes$habAv <- Bovar_6KmDispDist_habAv
# Bovar_4KmDispDist_habAv <- scan("Bovar_4KmDispDist_habAv.txt")
# Bovar_4Km_DispDist_topo_attributes$habAv <- Bovar_4KmDispDist_habAv
Bovar_8KmDispDist_habAv <- scan("Bovar_8KmDispDist_habAv.txt")
Bovar_8Km_DispDist_topo_attributes$habAv <- Bovar_8KmDispDist_habAv

# Hyarb
Hyarb_2KmDispDist_habAv <- scan("Hyarb_2KmDispDist_habAv.txt")
Hyarb_2Km_DispDist_topo_attributes$habAv <- Hyarb_2KmDispDist_habAv
names(Hyarb_2Km_DispDist_topo_attributes)
Hyarb_300mDispDist_habAv <- scan("Hyarb_300mDispDist_habAv.txt")
Hyarb_300m_DispDist_topo_attributes$habAv <- Hyarb_300mDispDist_habAv
names(Hyarb_300m_DispDist_topo_attributes)
Hyarb_10KmDispDist_habAv <- scan("Hyarb_10KmDispDist_habAv.txt")
Hyarb_10Km_DispDist_topo_attributes$habAv <- Hyarb_10KmDispDist_habAv
names(Hyarb_10Km_DispDist_topo_attributes)
Hyarb_1KmDispDist_habAv <- scan("Hyarb_1KmDispDist_habAv.txt")
Hyarb_1Km_DispDist_topo_attributes$habAv <- Hyarb_1KmDispDist_habAv
Hyarb_6KmDispDist_habAv <- scan("Hyarb_6KmDispDist_habAv.txt")
Hyarb_6Km_DispDist_topo_attributes$habAv <- Hyarb_6KmDispDist_habAv
Hyarb_4KmDispDist_habAv <- scan("Hyarb_4KmDispDist_habAv.txt")
Hyarb_4Km_DispDist_topo_attributes$habAv <- Hyarb_4KmDispDist_habAv
Hyarb_8KmDispDist_habAv <- scan("Hyarb_8KmDispDist_habAv.txt")
Hyarb_8Km_DispDist_topo_attributes$habAv <- Hyarb_8KmDispDist_habAv

# Alobs
Alobs_2KmDispDist_habAv <- scan("Alobs_2KmDispDist_habAv.txt")
Alobs_2Km_DispDist_topo_attributes$habAv <- Alobs_2KmDispDist_habAv
names(Alobs_2Km_DispDist_topo_attributes)
Alobs_300mDispDist_habAv <- scan("Alobs_300mDispDist_habAv.txt")
Alobs_300m_DispDist_topo_attributes$habAv <- Alobs_300mDispDist_habAv
names(Alobs_300m_DispDist_topo_attributes)
Alobs_10KmDispDist_habAv <- scan("Alobs_10KmDispDist_habAv.txt")
Alobs_10Km_DispDist_topo_attributes$habAv <- Alobs_10KmDispDist_habAv
names(Alobs_10Km_DispDist_topo_attributes)
Alobs_1KmDispDist_habAv <- scan("Alobs_1KmDispDist_habAv.txt")
Alobs_1Km_DispDist_topo_attributes$habAv <- Alobs_1KmDispDist_habAv
Alobs_6KmDispDist_habAv <- scan("Alobs_6KmDispDist_habAv.txt")
Alobs_6Km_DispDist_topo_attributes$habAv <- Alobs_6KmDispDist_habAv
Alobs_4KmDispDist_habAv <- scan("Alobs_4KmDispDist_habAv.txt")
Alobs_4Km_DispDist_topo_attributes$habAv <- Alobs_4KmDispDist_habAv
Alobs_8KmDispDist_habAv <- scan("Alobs_8KmDispDist_habAv.txt")
Alobs_8Km_DispDist_topo_attributes$habAv <- Alobs_8KmDispDist_habAv

# Epcal
Epcal_2KmDispDist_habAv <- scan("Epcal_2KmDispDist_habAv.txt")
Epcal_2Km_DispDist_topo_attributes$habAv <- Epcal_2KmDispDist_habAv
names(Epcal_2Km_DispDist_topo_attributes)
Epcal_300mDispDist_habAv <- scan("Epcal_300mDispDist_habAv.txt")
Epcal_300m_DispDist_topo_attributes$habAv <- Epcal_300mDispDist_habAv
names(Epcal_300m_DispDist_topo_attributes)
Epcal_10KmDispDist_habAv <- scan("Epcal_10KmDispDist_habAv.txt")
Epcal_10Km_DispDist_topo_attributes$habAv <- Epcal_10KmDispDist_habAv
names(Epcal_10Km_DispDist_topo_attributes)
Epcal_1KmDispDist_habAv <- scan("Epcal_1KmDispDist_habAv.txt")
Epcal_1Km_DispDist_topo_attributes$habAv <- Epcal_1KmDispDist_habAv
Epcal_6KmDispDist_habAv <- scan("Epcal_6KmDispDist_habAv.txt")
Epcal_6Km_DispDist_topo_attributes$habAv <- Epcal_6KmDispDist_habAv
Epcal_4KmDispDist_habAv <- scan("Epcal_4KmDispDist_habAv.txt")
Epcal_4Km_DispDist_topo_attributes$habAv <- Epcal_4KmDispDist_habAv
Epcal_8KmDispDist_habAv <- scan("Epcal_8KmDispDist_habAv.txt")
Epcal_8Km_DispDist_topo_attributes$habAv <- Epcal_8KmDispDist_habAv

# Peagg
Peagg_2KmDispDist_habAv <- scan("Peagg_2KmDispDist_habAv.txt")
Peagg_2Km_DispDist_topo_attributes$habAv <- Peagg_2KmDispDist_habAv
names(Peagg_2Km_DispDist_topo_attributes)
Peagg_300mDispDist_habAv <- scan("Peagg_300mDispDist_habAv.txt")
Peagg_300m_DispDist_topo_attributes$habAv <- Peagg_300mDispDist_habAv
names(Peagg_300m_DispDist_topo_attributes)
Peagg_10KmDispDist_habAv <- scan("Peagg_10KmDispDist_habAv.txt")
Peagg_10Km_DispDist_topo_attributes$habAv <- Peagg_10KmDispDist_habAv
names(Peagg_10Km_DispDist_topo_attributes)
Peagg_1KmDispDist_habAv <- scan("Peagg_1KmDispDist_habAv.txt")
Peagg_1Km_DispDist_topo_attributes$habAv <- Peagg_1KmDispDist_habAv
Peagg_6KmDispDist_habAv <- scan("Peagg_6KmDispDist_habAv.txt")
Peagg_6Km_DispDist_topo_attributes$habAv <- Peagg_6KmDispDist_habAv
Peagg_4KmDispDist_habAv <- scan("Peagg_4KmDispDist_habAv.txt")
Peagg_4Km_DispDist_topo_attributes$habAv <- Peagg_4KmDispDist_habAv
Peagg_8KmDispDist_habAv <- scan("Peagg_8KmDispDist_habAv.txt")
Peagg_8Km_DispDist_topo_attributes$habAv <- Peagg_8KmDispDist_habAv

# Perid
Perid_2KmDispDist_habAv <- scan("Perid_2KmDispDist_habAv.txt")
Perid_2Km_DispDist_topo_attributes$habAv <- Perid_2KmDispDist_habAv
names(Perid_2Km_DispDist_topo_attributes)
Perid_300mDispDist_habAv <- scan("Perid_300mDispDist_habAv.txt")
Perid_300m_DispDist_topo_attributes$habAv <- Perid_300mDispDist_habAv
names(Perid_300m_DispDist_topo_attributes)
Perid_10KmDispDist_habAv <- scan("Perid_10KmDispDist_habAv.txt")
Perid_10Km_DispDist_topo_attributes$habAv <- Perid_10KmDispDist_habAv
names(Perid_10Km_DispDist_topo_attributes)
Perid_1KmDispDist_habAv <- scan("Perid_1KmDispDist_habAv.txt")
Perid_1Km_DispDist_topo_attributes$habAv <- Perid_1KmDispDist_habAv
Perid_6KmDispDist_habAv <- scan("Perid_6KmDispDist_habAv.txt")
Perid_6Km_DispDist_topo_attributes$habAv <- Perid_6KmDispDist_habAv
Perid_4KmDispDist_habAv <- scan("Perid_4KmDispDist_habAv.txt")
Perid_4Km_DispDist_topo_attributes$habAv <- Perid_4KmDispDist_habAv
Perid_8KmDispDist_habAv <- scan("Perid_8KmDispDist_habAv.txt")
Perid_8Km_DispDist_topo_attributes$habAv <- Perid_8KmDispDist_habAv



#### Do a subset of the columns that only has the predictors and response for BRT fitting####
#save it still in the 'networks' folder
# Changed name of V_Ge_spe for each to V_s for all, to have it more generic

#V_s should be the proper name of vector V_Bo_var
#Bovar keeps denomination V_Bo_var instead of V_s in the 2km, 10km, and 300m models, it has no error because it was the template for the others
Bovar_2KmDispDist_stattest <- subset(Bovar_2Km_DispDist_topo_attributes, select=c("PatchID", "V_Bo_var", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Bovar_2KmDispDist_stattest)
write.csv(Bovar_2KmDispDist_stattest, file = "Bovar_2KmDispDist_stattest.csv")

Bovar_300mDispDist_stattest <- subset(Bovar_300m_DispDist_topo_attributes, select=c("PatchID", "V_Bo_var", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Bovar_300mDispDist_stattest)
write.csv(Bovar_300mDispDist_stattest, file = "Bovar_300mDispDist_stattest.csv")

Bovar_10KmDispDist_stattest <- subset(Bovar_10Km_DispDist_topo_attributes, select=c("PatchID", "V_Bo_var", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Bovar_10KmDispDist_stattest)
write.csv(Bovar_10KmDispDist_stattest, file = "Bovar_10KmDispDist_stattest.csv")

Bovar_1KmDispDist_stattest <- subset(Bovar_1Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Bovar_1KmDispDist_stattest)
write.csv(Bovar_1KmDispDist_stattest, file = "Bovar_1KmDispDist_stattest.csv")

Bovar_6KmDispDist_stattest <- subset(Bovar_6Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Bovar_6KmDispDist_stattest)
write.csv(Bovar_6KmDispDist_stattest, file = "Bovar_6KmDispDist_stattest.csv")

# Bovar_4KmDispDist_stattest <- subset(Bovar_4Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
# names(Bovar_4KmDispDist_stattest)
# write.csv(Bovar_4KmDispDist_stattest, file = "Bovar_4KmDispDist_stattest.csv")

Bovar_8KmDispDist_stattest <- subset(Bovar_8Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Bovar_8KmDispDist_stattest)
write.csv(Bovar_8KmDispDist_stattest, file = "Bovar_8KmDispDist_stattest.csv")


#Hyarb
Hyarb_2KmDispDist_stattest <- subset(Hyarb_2Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Hyarb_2KmDispDist_stattest)
write.csv(Hyarb_2KmDispDist_stattest, file = "Hyarb_2KmDispDist_stattest.csv")

Hyarb_300mDispDist_stattest <- subset(Hyarb_300m_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Hyarb_300mDispDist_stattest)
write.csv(Hyarb_300mDispDist_stattest, file = "Hyarb_300mDispDist_stattest.csv")

Hyarb_10KmDispDist_stattest <- subset(Hyarb_10Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Hyarb_10KmDispDist_stattest)
write.csv(Hyarb_10KmDispDist_stattest, file = "Hyarb_10KmDispDist_stattest.csv")

Hyarb_1KmDispDist_stattest <- subset(Hyarb_1Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Hyarb_1KmDispDist_stattest)
write.csv(Hyarb_1KmDispDist_stattest, file = "Hyarb_1KmDispDist_stattest.csv")

Hyarb_6KmDispDist_stattest <- subset(Hyarb_6Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Hyarb_6KmDispDist_stattest)
write.csv(Hyarb_6KmDispDist_stattest, file = "Hyarb_6KmDispDist_stattest.csv")

Hyarb_4KmDispDist_stattest <- subset(Hyarb_4Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Hyarb_4KmDispDist_stattest)
write.csv(Hyarb_4KmDispDist_stattest, file = "Hyarb_4KmDispDist_stattest.csv")

Hyarb_8KmDispDist_stattest <- subset(Hyarb_8Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Hyarb_8KmDispDist_stattest)
write.csv(Hyarb_8KmDispDist_stattest, file = "Hyarb_8KmDispDist_stattest.csv")


#Alobs
Alobs_2KmDispDist_stattest <- subset(Alobs_2Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Alobs_2KmDispDist_stattest)
write.csv(Alobs_2KmDispDist_stattest, file = "Alobs_2KmDispDist_stattest.csv")

Alobs_300mDispDist_stattest <- subset(Alobs_300m_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Alobs_300mDispDist_stattest)
write.csv(Alobs_300mDispDist_stattest, file = "Alobs_300mDispDist_stattest.csv")

Alobs_10KmDispDist_stattest <- subset(Alobs_10Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Alobs_10KmDispDist_stattest)
write.csv(Alobs_10KmDispDist_stattest, file = "Alobs_10KmDispDist_stattest.csv")

Alobs_1KmDispDist_stattest <- subset(Alobs_1Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Alobs_1KmDispDist_stattest)
write.csv(Alobs_1KmDispDist_stattest, file = "Alobs_1KmDispDist_stattest.csv")

Alobs_6KmDispDist_stattest <- subset(Alobs_6Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Alobs_6KmDispDist_stattest)
write.csv(Alobs_6KmDispDist_stattest, file = "Alobs_6KmDispDist_stattest.csv")

Alobs_4KmDispDist_stattest <- subset(Alobs_4Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Alobs_4KmDispDist_stattest)
write.csv(Alobs_4KmDispDist_stattest, file = "Alobs_4KmDispDist_stattest.csv")

Alobs_8KmDispDist_stattest <- subset(Alobs_8Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Alobs_8KmDispDist_stattest)
write.csv(Alobs_8KmDispDist_stattest, file = "Alobs_8KmDispDist_stattest.csv")


#Epcal
Epcal_2KmDispDist_stattest <- subset(Epcal_2Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Epcal_2KmDispDist_stattest)
write.csv(Epcal_2KmDispDist_stattest, file = "Epcal_2KmDispDist_stattest.csv")

Epcal_300mDispDist_stattest <- subset(Epcal_300m_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Epcal_300mDispDist_stattest)
write.csv(Epcal_300mDispDist_stattest, file = "Epcal_300mDispDist_stattest.csv")

Epcal_10KmDispDist_stattest <- subset(Epcal_10Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Epcal_10KmDispDist_stattest)
write.csv(Epcal_10KmDispDist_stattest, file = "Epcal_10KmDispDist_stattest.csv")

Epcal_1KmDispDist_stattest <- subset(Epcal_1Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Epcal_1KmDispDist_stattest)
write.csv(Epcal_1KmDispDist_stattest, file = "Epcal_1KmDispDist_stattest.csv")

Epcal_6KmDispDist_stattest <- subset(Epcal_6Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Epcal_6KmDispDist_stattest)
write.csv(Epcal_6KmDispDist_stattest, file = "Epcal_6KmDispDist_stattest.csv")

Epcal_4KmDispDist_stattest <- subset(Epcal_4Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Epcal_4KmDispDist_stattest)
write.csv(Epcal_4KmDispDist_stattest, file = "Epcal_4KmDispDist_stattest.csv")

Epcal_8KmDispDist_stattest <- subset(Epcal_8Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Epcal_8KmDispDist_stattest)
write.csv(Epcal_8KmDispDist_stattest, file = "Epcal_8KmDispDist_stattest.csv")


#Peagg
Peagg_2KmDispDist_stattest <- subset(Peagg_2Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Peagg_2KmDispDist_stattest)
write.csv(Peagg_2KmDispDist_stattest, file = "Peagg_2KmDispDist_stattest.csv")

Peagg_300mDispDist_stattest <- subset(Peagg_300m_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Peagg_300mDispDist_stattest)
write.csv(Peagg_300mDispDist_stattest, file = "Peagg_300mDispDist_stattest.csv")

Peagg_10KmDispDist_stattest <- subset(Peagg_10Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Peagg_10KmDispDist_stattest)
write.csv(Peagg_10KmDispDist_stattest, file = "Peagg_10KmDispDist_stattest.csv")

Peagg_1KmDispDist_stattest <- subset(Peagg_1Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Peagg_1KmDispDist_stattest)
write.csv(Peagg_1KmDispDist_stattest, file = "Peagg_1KmDispDist_stattest.csv")

Peagg_6KmDispDist_stattest <- subset(Peagg_6Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Peagg_6KmDispDist_stattest)
write.csv(Peagg_6KmDispDist_stattest, file = "Peagg_6KmDispDist_stattest.csv")

Peagg_4KmDispDist_stattest <- subset(Peagg_4Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Peagg_4KmDispDist_stattest)
write.csv(Peagg_4KmDispDist_stattest, file = "Peagg_4KmDispDist_stattest.csv")

Peagg_8KmDispDist_stattest <- subset(Peagg_8Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Peagg_8KmDispDist_stattest)
write.csv(Peagg_8KmDispDist_stattest, file = "Peagg_8KmDispDist_stattest.csv")


#Perid
Perid_2KmDispDist_stattest <- subset(Perid_2Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Perid_2KmDispDist_stattest)
write.csv(Perid_2KmDispDist_stattest, file = "Perid_2KmDispDist_stattest.csv")

Perid_300mDispDist_stattest <- subset(Perid_300m_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Perid_300mDispDist_stattest)
write.csv(Perid_300mDispDist_stattest, file = "Perid_300mDispDist_stattest.csv")

Perid_10KmDispDist_stattest <- subset(Perid_10Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Perid_10KmDispDist_stattest)
write.csv(Perid_10KmDispDist_stattest, file = "Perid_10KmDispDist_stattest.csv")

Perid_1KmDispDist_stattest <- subset(Perid_1Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Perid_1KmDispDist_stattest)
write.csv(Perid_1KmDispDist_stattest, file = "Perid_1KmDispDist_stattest.csv")

Perid_6KmDispDist_stattest <- subset(Perid_6Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Perid_6KmDispDist_stattest)
write.csv(Perid_6KmDispDist_stattest, file = "Perid_6KmDispDist_stattest.csv")

Perid_4KmDispDist_stattest <- subset(Perid_4Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Perid_4KmDispDist_stattest)
write.csv(Perid_4KmDispDist_stattest, file = "Perid_4KmDispDist_stattest.csv")

Perid_8KmDispDist_stattest <- subset(Perid_8Km_DispDist_topo_attributes, select=c("PatchID", "V_s", "V_t", "Patch_Area", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Perid_8KmDispDist_stattest)
write.csv(Perid_8KmDispDist_stattest, file = "Perid_8KmDispDist_stattest.csv")



##### Do & save model-fitting subsets (noNA, only 0s and 1s) ####
#### only pres_abs (no-NoData) df's #########################
setwd('C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs')
#Take out NA's, set threshold of absences, change name and value ####################

#Bovar
Bovar_2KmDispDist_woNA = Bovar_2KmDispDist_stattest
Bovar_2KmDispDist_woNA[Bovar_2KmDispDist_woNA$V_Bo_var==0 & Bovar_2KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Bovar_2KmDispDist_woNA = Bovar_2KmDispDist_woNA[!is.na(Bovar_2KmDispDist_woNA$pres_abs),]
write.csv(Bovar_2KmDispDist_woNA, file = "Bovar_2KmDispDist_stattest_woNA.csv")

Bovar_300mDispDist_woNA = Bovar_300mDispDist_stattest
Bovar_300mDispDist_woNA[Bovar_300mDispDist_woNA$V_Bo_var==0 & Bovar_300mDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Bovar_300mDispDist_woNA = Bovar_300mDispDist_woNA[!is.na(Bovar_300mDispDist_woNA$pres_abs),]
write.csv(Bovar_300mDispDist_woNA, file = "Bovar_300mDispDist_stattest_woNA.csv")

Bovar_10KmDispDist_woNA = Bovar_10KmDispDist_stattest
Bovar_10KmDispDist_woNA[Bovar_10KmDispDist_woNA$V_Bo_var==0 & Bovar_10KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Bovar_10KmDispDist_woNA = Bovar_10KmDispDist_woNA[!is.na(Bovar_10KmDispDist_woNA$pres_abs),]
write.csv(Bovar_10KmDispDist_woNA, file = "Bovar_10KmDispDist_stattest_woNA.csv")

Bovar_1KmDispDist_woNA = Bovar_1KmDispDist_stattest
Bovar_1KmDispDist_woNA[Bovar_1KmDispDist_woNA$V_Bo_var==0 & Bovar_1KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Bovar_1KmDispDist_woNA = Bovar_1KmDispDist_woNA[!is.na(Bovar_1KmDispDist_woNA$pres_abs),]
write.csv(Bovar_1KmDispDist_woNA, file = "Bovar_1KmDispDist_stattest_woNA.csv")

Bovar_6KmDispDist_woNA = Bovar_6KmDispDist_stattest
Bovar_6KmDispDist_woNA[Bovar_6KmDispDist_woNA$V_Bo_var==0 & Bovar_6KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Bovar_6KmDispDist_woNA = Bovar_6KmDispDist_woNA[!is.na(Bovar_6KmDispDist_woNA$pres_abs),]
write.csv(Bovar_6KmDispDist_woNA, file = "Bovar_6KmDispDist_stattest_woNA.csv")

# Bovar_4KmDispDist_woNA = Bovar_4KmDispDist_stattest
# Bovar_4KmDispDist_woNA[Bovar_4KmDispDist_woNA$V_Bo_var==0 & Bovar_4KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
# Bovar_4KmDispDist_woNA = Bovar_4KmDispDist_woNA[!is.na(Bovar_4KmDispDist_woNA$pres_abs),]
# write.csv(Bovar_4KmDispDist_woNA, file = "Bovar_4KmDispDist_stattest_woNA.csv")

Bovar_8KmDispDist_woNA = Bovar_8KmDispDist_stattest
Bovar_8KmDispDist_woNA[Bovar_8KmDispDist_woNA$V_Bo_var==0 & Bovar_8KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Bovar_8KmDispDist_woNA = Bovar_8KmDispDist_woNA[!is.na(Bovar_8KmDispDist_woNA$pres_abs),]
write.csv(Bovar_8KmDispDist_woNA, file = "Bovar_8KmDispDist_stattest_woNA.csv")


#Hyarb
Hyarb_2KmDispDist_woNA = Hyarb_2KmDispDist_stattest
Hyarb_2KmDispDist_woNA[Hyarb_2KmDispDist_woNA$V_s==0 & Hyarb_2KmDispDist_woNA$V_t >= 6, 'pres_abs'] = 0
Hyarb_2KmDispDist_woNA = Hyarb_2KmDispDist_woNA[!is.na(Hyarb_2KmDispDist_woNA$pres_abs),]
write.csv(Hyarb_2KmDispDist_woNA, file = "Hyarb_2KmDispDist_stattest_woNA.csv")

Hyarb_300mDispDist_woNA = Hyarb_300mDispDist_stattest
Hyarb_300mDispDist_woNA[Hyarb_300mDispDist_woNA$V_s==0 & Hyarb_300mDispDist_woNA$V_t >= 6, 'pres_abs'] = 0
Hyarb_300mDispDist_woNA = Hyarb_300mDispDist_woNA[!is.na(Hyarb_300mDispDist_woNA$pres_abs),]
write.csv(Hyarb_300mDispDist_woNA, file = "Hyarb_300mDispDist_stattest_woNA.csv")

Hyarb_10KmDispDist_woNA = Hyarb_10KmDispDist_stattest
Hyarb_10KmDispDist_woNA[Hyarb_10KmDispDist_woNA$V_s==0 & Hyarb_10KmDispDist_woNA$V_t >= 6, 'pres_abs'] = 0
Hyarb_10KmDispDist_woNA = Hyarb_10KmDispDist_woNA[!is.na(Hyarb_10KmDispDist_woNA$pres_abs),]
write.csv(Hyarb_10KmDispDist_woNA, file = "Hyarb_10KmDispDist_stattest_woNA.csv")

Hyarb_1KmDispDist_woNA = Hyarb_1KmDispDist_stattest
Hyarb_1KmDispDist_woNA[Hyarb_1KmDispDist_woNA$V_s==0 & Hyarb_1KmDispDist_woNA$V_t >= 6, 'pres_abs'] = 0
Hyarb_1KmDispDist_woNA = Hyarb_1KmDispDist_woNA[!is.na(Hyarb_1KmDispDist_woNA$pres_abs),]
write.csv(Hyarb_1KmDispDist_woNA, file = "Hyarb_1KmDispDist_stattest_woNA.csv")

Hyarb_6KmDispDist_woNA = Hyarb_6KmDispDist_stattest
Hyarb_6KmDispDist_woNA[Hyarb_6KmDispDist_woNA$V_s==0 & Hyarb_6KmDispDist_woNA$V_t >= 6, 'pres_abs'] = 0
Hyarb_6KmDispDist_woNA = Hyarb_6KmDispDist_woNA[!is.na(Hyarb_6KmDispDist_woNA$pres_abs),]
write.csv(Hyarb_6KmDispDist_woNA, file = "Hyarb_6KmDispDist_stattest_woNA.csv")

Hyarb_4KmDispDist_woNA = Hyarb_4KmDispDist_stattest
Hyarb_4KmDispDist_woNA[Hyarb_4KmDispDist_woNA$V_s==0 & Hyarb_4KmDispDist_woNA$V_t >= 6, 'pres_abs'] = 0
Hyarb_4KmDispDist_woNA = Hyarb_4KmDispDist_woNA[!is.na(Hyarb_4KmDispDist_woNA$pres_abs),]
write.csv(Hyarb_4KmDispDist_woNA, file = "Hyarb_4KmDispDist_stattest_woNA.csv")

Hyarb_8KmDispDist_woNA = Hyarb_8KmDispDist_stattest
Hyarb_8KmDispDist_woNA[Hyarb_8KmDispDist_woNA$V_s==0 & Hyarb_8KmDispDist_woNA$V_t >= 6, 'pres_abs'] = 0
Hyarb_8KmDispDist_woNA = Hyarb_8KmDispDist_woNA[!is.na(Hyarb_8KmDispDist_woNA$pres_abs),]
write.csv(Hyarb_8KmDispDist_woNA, file = "Hyarb_8KmDispDist_stattest_woNA.csv")


#Alobs
Alobs_2KmDispDist_woNA = Alobs_2KmDispDist_stattest
Alobs_2KmDispDist_woNA[Alobs_2KmDispDist_woNA$V_s==0 & Alobs_2KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Alobs_2KmDispDist_woNA = Alobs_2KmDispDist_woNA[!is.na(Alobs_2KmDispDist_woNA$pres_abs),]
write.csv(Alobs_2KmDispDist_woNA, file = "Alobs_2KmDispDist_stattest_woNA.csv")

Alobs_300mDispDist_woNA = Alobs_300mDispDist_stattest
Alobs_300mDispDist_woNA[Alobs_300mDispDist_woNA$V_s==0 & Alobs_300mDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Alobs_300mDispDist_woNA = Alobs_300mDispDist_woNA[!is.na(Alobs_300mDispDist_woNA$pres_abs),]
write.csv(Alobs_300mDispDist_woNA, file = "Alobs_300mDispDist_stattest_woNA.csv")

Alobs_10KmDispDist_woNA = Alobs_10KmDispDist_stattest
Alobs_10KmDispDist_woNA[Alobs_10KmDispDist_woNA$V_s==0 & Alobs_10KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Alobs_10KmDispDist_woNA = Alobs_10KmDispDist_woNA[!is.na(Alobs_10KmDispDist_woNA$pres_abs),]
write.csv(Alobs_10KmDispDist_woNA, file = "Alobs_10KmDispDist_stattest_woNA.csv")

Alobs_1KmDispDist_woNA = Alobs_1KmDispDist_stattest
Alobs_1KmDispDist_woNA[Alobs_1KmDispDist_woNA$V_s==0 & Alobs_1KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Alobs_1KmDispDist_woNA = Alobs_1KmDispDist_woNA[!is.na(Alobs_1KmDispDist_woNA$pres_abs),]
write.csv(Alobs_1KmDispDist_woNA, file = "Alobs_1KmDispDist_stattest_woNA.csv")

Alobs_6KmDispDist_woNA = Alobs_6KmDispDist_stattest
Alobs_6KmDispDist_woNA[Alobs_6KmDispDist_woNA$V_s==0 & Alobs_6KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Alobs_6KmDispDist_woNA = Alobs_6KmDispDist_woNA[!is.na(Alobs_6KmDispDist_woNA$pres_abs),]
write.csv(Alobs_6KmDispDist_woNA, file = "Alobs_6KmDispDist_stattest_woNA.csv")

Alobs_4KmDispDist_woNA = Alobs_4KmDispDist_stattest
Alobs_4KmDispDist_woNA[Alobs_4KmDispDist_woNA$V_s==0 & Alobs_4KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Alobs_4KmDispDist_woNA = Alobs_4KmDispDist_woNA[!is.na(Alobs_4KmDispDist_woNA$pres_abs),]
write.csv(Alobs_4KmDispDist_woNA, file = "Alobs_4KmDispDist_stattest_woNA.csv")

Alobs_8KmDispDist_woNA = Alobs_8KmDispDist_stattest
Alobs_8KmDispDist_woNA[Alobs_8KmDispDist_woNA$V_s==0 & Alobs_8KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Alobs_8KmDispDist_woNA = Alobs_8KmDispDist_woNA[!is.na(Alobs_8KmDispDist_woNA$pres_abs),]
write.csv(Alobs_8KmDispDist_woNA, file = "Alobs_8KmDispDist_stattest_woNA.csv")


#Epcal
Epcal_2KmDispDist_woNA = Epcal_2KmDispDist_stattest
Epcal_2KmDispDist_woNA[Epcal_2KmDispDist_woNA$V_s==0 & Epcal_2KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Epcal_2KmDispDist_woNA = Epcal_2KmDispDist_woNA[!is.na(Epcal_2KmDispDist_woNA$pres_abs),]
write.csv(Epcal_2KmDispDist_woNA, file = "Epcal_2KmDispDist_stattest_woNA.csv")

Epcal_300mDispDist_woNA = Epcal_300mDispDist_stattest
Epcal_300mDispDist_woNA[Epcal_300mDispDist_woNA$V_s==0 & Epcal_300mDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Epcal_300mDispDist_woNA = Epcal_300mDispDist_woNA[!is.na(Epcal_300mDispDist_woNA$pres_abs),]
write.csv(Epcal_300mDispDist_woNA, file = "Epcal_300mDispDist_stattest_woNA.csv")

Epcal_10KmDispDist_woNA = Epcal_10KmDispDist_stattest
Epcal_10KmDispDist_woNA[Epcal_10KmDispDist_woNA$V_s==0 & Epcal_10KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Epcal_10KmDispDist_woNA = Epcal_10KmDispDist_woNA[!is.na(Epcal_10KmDispDist_woNA$pres_abs),]
write.csv(Epcal_10KmDispDist_woNA, file = "Epcal_10KmDispDist_stattest_woNA.csv")

Epcal_1KmDispDist_woNA = Epcal_1KmDispDist_stattest
Epcal_1KmDispDist_woNA[Epcal_1KmDispDist_woNA$V_s==0 & Epcal_1KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Epcal_1KmDispDist_woNA = Epcal_1KmDispDist_woNA[!is.na(Epcal_1KmDispDist_woNA$pres_abs),]
write.csv(Epcal_1KmDispDist_woNA, file = "Epcal_1KmDispDist_stattest_woNA.csv")

Epcal_6KmDispDist_woNA = Epcal_6KmDispDist_stattest
Epcal_6KmDispDist_woNA[Epcal_6KmDispDist_woNA$V_s==0 & Epcal_6KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Epcal_6KmDispDist_woNA = Epcal_6KmDispDist_woNA[!is.na(Epcal_6KmDispDist_woNA$pres_abs),]
write.csv(Epcal_6KmDispDist_woNA, file = "Epcal_6KmDispDist_stattest_woNA.csv")

Epcal_4KmDispDist_woNA = Epcal_4KmDispDist_stattest
Epcal_4KmDispDist_woNA[Epcal_4KmDispDist_woNA$V_s==0 & Epcal_4KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Epcal_4KmDispDist_woNA = Epcal_4KmDispDist_woNA[!is.na(Epcal_4KmDispDist_woNA$pres_abs),]
write.csv(Epcal_4KmDispDist_woNA, file = "Epcal_4KmDispDist_stattest_woNA.csv")

Epcal_8KmDispDist_woNA = Epcal_8KmDispDist_stattest
Epcal_8KmDispDist_woNA[Epcal_8KmDispDist_woNA$V_s==0 & Epcal_8KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Epcal_8KmDispDist_woNA = Epcal_8KmDispDist_woNA[!is.na(Epcal_8KmDispDist_woNA$pres_abs),]
write.csv(Epcal_8KmDispDist_woNA, file = "Epcal_8KmDispDist_stattest_woNA.csv")


#Peagg 
Peagg_2KmDispDist_woNA = Peagg_2KmDispDist_stattest
Peagg_2KmDispDist_woNA[Peagg_2KmDispDist_woNA$V_s==0 & Peagg_2KmDispDist_woNA$V_t >= 3, 'pres_abs'] = 0
Peagg_2KmDispDist_woNA = Peagg_2KmDispDist_woNA[!is.na(Peagg_2KmDispDist_woNA$pres_abs),]
write.csv(Peagg_2KmDispDist_woNA, file = "Peagg_2KmDispDist_stattest_woNA.csv")

Peagg_300mDispDist_woNA = Peagg_300mDispDist_stattest
Peagg_300mDispDist_woNA[Peagg_300mDispDist_woNA$V_s==0 & Peagg_300mDispDist_woNA$V_t >= 3, 'pres_abs'] = 0
Peagg_300mDispDist_woNA = Peagg_300mDispDist_woNA[!is.na(Peagg_300mDispDist_woNA$pres_abs),]
write.csv(Peagg_300mDispDist_woNA, file = "Peagg_300mDispDist_stattest_woNA.csv")

Peagg_10KmDispDist_woNA = Peagg_10KmDispDist_stattest
Peagg_10KmDispDist_woNA[Peagg_10KmDispDist_woNA$V_s==0 & Peagg_10KmDispDist_woNA$V_t >= 3, 'pres_abs'] = 0
Peagg_10KmDispDist_woNA = Peagg_10KmDispDist_woNA[!is.na(Peagg_10KmDispDist_woNA$pres_abs),]
write.csv(Peagg_10KmDispDist_woNA, file = "Peagg_10KmDispDist_stattest_woNA.csv")

Peagg_1KmDispDist_woNA = Peagg_1KmDispDist_stattest
Peagg_1KmDispDist_woNA[Peagg_1KmDispDist_woNA$V_s==0 & Peagg_1KmDispDist_woNA$V_t >= 3, 'pres_abs'] = 0
Peagg_1KmDispDist_woNA = Peagg_1KmDispDist_woNA[!is.na(Peagg_1KmDispDist_woNA$pres_abs),]
write.csv(Peagg_1KmDispDist_woNA, file = "Peagg_1KmDispDist_stattest_woNA.csv")

Peagg_6KmDispDist_woNA = Peagg_6KmDispDist_stattest
Peagg_6KmDispDist_woNA[Peagg_6KmDispDist_woNA$V_s==0 & Peagg_6KmDispDist_woNA$V_t >= 3, 'pres_abs'] = 0
Peagg_6KmDispDist_woNA = Peagg_6KmDispDist_woNA[!is.na(Peagg_6KmDispDist_woNA$pres_abs),]
write.csv(Peagg_6KmDispDist_woNA, file = "Peagg_6KmDispDist_stattest_woNA.csv")

Peagg_4KmDispDist_woNA = Peagg_4KmDispDist_stattest
Peagg_4KmDispDist_woNA[Peagg_4KmDispDist_woNA$V_s==0 & Peagg_4KmDispDist_woNA$V_t >= 3, 'pres_abs'] = 0
Peagg_4KmDispDist_woNA = Peagg_4KmDispDist_woNA[!is.na(Peagg_4KmDispDist_woNA$pres_abs),]
write.csv(Peagg_4KmDispDist_woNA, file = "Peagg_4KmDispDist_stattest_woNA.csv")

Peagg_8KmDispDist_woNA = Peagg_8KmDispDist_stattest
Peagg_8KmDispDist_woNA[Peagg_8KmDispDist_woNA$V_s==0 & Peagg_8KmDispDist_woNA$V_t >= 3, 'pres_abs'] = 0
Peagg_8KmDispDist_woNA = Peagg_8KmDispDist_woNA[!is.na(Peagg_8KmDispDist_woNA$pres_abs),]
write.csv(Peagg_8KmDispDist_woNA, file = "Peagg_8KmDispDist_stattest_woNA.csv")


#Perid
Perid_2KmDispDist_woNA = Perid_2KmDispDist_stattest
Perid_2KmDispDist_woNA[Perid_2KmDispDist_woNA$V_s==0 & Perid_2KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Perid_2KmDispDist_woNA = Perid_2KmDispDist_woNA[!is.na(Perid_2KmDispDist_woNA$pres_abs),]
write.csv(Perid_2KmDispDist_woNA, file = "Perid_2KmDispDist_stattest_woNA.csv")

Perid_300mDispDist_woNA = Perid_300mDispDist_stattest
Perid_300mDispDist_woNA[Perid_300mDispDist_woNA$V_s==0 & Perid_300mDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Perid_300mDispDist_woNA = Perid_300mDispDist_woNA[!is.na(Perid_300mDispDist_woNA$pres_abs),]
write.csv(Perid_300mDispDist_woNA, file = "Perid_300mDispDist_stattest_woNA.csv")

Perid_10KmDispDist_woNA = Perid_10KmDispDist_stattest
Perid_10KmDispDist_woNA[Perid_10KmDispDist_woNA$V_s==0 & Perid_10KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Perid_10KmDispDist_woNA = Perid_10KmDispDist_woNA[!is.na(Perid_10KmDispDist_woNA$pres_abs),]
write.csv(Perid_10KmDispDist_woNA, file = "Perid_10KmDispDist_stattest_woNA.csv")

Perid_1KmDispDist_woNA = Perid_1KmDispDist_stattest
Perid_1KmDispDist_woNA[Perid_1KmDispDist_woNA$V_s==0 & Perid_1KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Perid_1KmDispDist_woNA = Perid_1KmDispDist_woNA[!is.na(Perid_1KmDispDist_woNA$pres_abs),]
write.csv(Perid_1KmDispDist_woNA, file = "Perid_1KmDispDist_stattest_woNA.csv")

Perid_6KmDispDist_woNA = Perid_6KmDispDist_stattest
Perid_6KmDispDist_woNA[Perid_6KmDispDist_woNA$V_s==0 & Perid_6KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Perid_6KmDispDist_woNA = Perid_6KmDispDist_woNA[!is.na(Perid_6KmDispDist_woNA$pres_abs),]
write.csv(Perid_6KmDispDist_woNA, file = "Perid_6KmDispDist_stattest_woNA.csv")

Perid_4KmDispDist_woNA = Perid_4KmDispDist_stattest
Perid_4KmDispDist_woNA[Perid_4KmDispDist_woNA$V_s==0 & Perid_4KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Perid_4KmDispDist_woNA = Perid_4KmDispDist_woNA[!is.na(Perid_4KmDispDist_woNA$pres_abs),]
write.csv(Perid_4KmDispDist_woNA, file = "Perid_4KmDispDist_stattest_woNA.csv")

Perid_8KmDispDist_woNA = Perid_8KmDispDist_stattest
Perid_8KmDispDist_woNA[Perid_8KmDispDist_woNA$V_s==0 & Perid_8KmDispDist_woNA$V_t >= 4, 'pres_abs'] = 0
Perid_8KmDispDist_woNA = Perid_8KmDispDist_woNA[!is.na(Perid_8KmDispDist_woNA$pres_abs),]
write.csv(Perid_8KmDispDist_woNA, file = "Perid_8KmDispDist_stattest_woNA.csv")




##################################################################################################################
#### Do BRT's for newly developed networks ####
##################################################################################################################

#Hyperparameters for all models
lr = 0.001
tc = 5
bf = 0.75
n_repeats = 100


################################################################################
##### Bovar ####################################################################


#### 2Km #######################################################################

# Create the output table that contains all the values of each of the runs
Bovar_2KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar

# Create vectors to record performance measures
Bovar_2KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Bovar_2KmDispDist_AUC_vec = vector() #

Bovar_2KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Bovar_2KmDispDist_HSI_vec = vector() 
Bovar_2KmDispDist_EgoSize_vec = vector() 
Bovar_2KmDispDist_strength_vec  = vector()
Bovar_2KmDispDist_deg_vec  = vector()
Bovar_2KmDispDist_habAv_vec  = vector()
Bovar_2KmDispDist_unw_b_c_vec  = vector()
Bovar_2KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Bovar_2KmDispDist_predict_mat = matrix(nrow = length(Bovar_2KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Bovar_2KmDispDist = gbm.step(data=Bovar_2KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Bovar_2KmDispDist)$var, imp = summary(gbm_mod_Bovar_2KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Bovar_2KmDispDist_HSI_vec = append(Bovar_2KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Bovar_2KmDispDist_EgoSize_vec = append(Bovar_2KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Bovar_2KmDispDist_strength_vec  = append(Bovar_2KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Bovar_2KmDispDist_deg_vec  = append(Bovar_2KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Bovar_2KmDispDist_habAv_vec  = append(Bovar_2KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Bovar_2KmDispDist_unw_b_c_vec  = append(Bovar_2KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Bovar_2KmDispDist_Patch_Area_vec  = append(Bovar_2KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Bovar_2KmDispDist_AUC_cv_vec = append(Bovar_2KmDispDist_AUC_cv_vec, gbm_mod_Bovar_2KmDispDist$cv.statistics$discrimination.mean)
  Bovar_2KmDispDist_AUC_vec = append(Bovar_2KmDispDist_AUC_vec, gbm_mod_Bovar_2KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Bovar_2KmDispDist$n.trees
  Bovar_2KmDispDist_nt_vec = append(Bovar_2KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Bovar_2KmDispDist_predict_mat[,i] = predict(gbm_mod_Bovar_2KmDispDist, Bovar_2KmDispDist_stattest, gbm_mod_Bovar_2KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Bovar_2KmDispDist_output_tab[,"AUC_train"] = Bovar_2KmDispDist_AUC_vec
Bovar_2KmDispDist_output_tab[,"AUC_cv"] = Bovar_2KmDispDist_AUC_cv_vec
Bovar_2KmDispDist_output_tab[,"ntrees"] = Bovar_2KmDispDist_nt_vec

#Var. importance columns
Bovar_2KmDispDist_output_tab[,"HSI_imp"] = Bovar_2KmDispDist_HSI_vec
Bovar_2KmDispDist_output_tab[,"EgoSize_imp"] = Bovar_2KmDispDist_EgoSize_vec
Bovar_2KmDispDist_output_tab[,"strength_imp"] = Bovar_2KmDispDist_strength_vec
Bovar_2KmDispDist_output_tab[,"deg_imp"] = Bovar_2KmDispDist_deg_vec
Bovar_2KmDispDist_output_tab[,"habAv_imp"] = Bovar_2KmDispDist_habAv_vec
Bovar_2KmDispDist_output_tab[,"unw_b_c_imp"] = Bovar_2KmDispDist_unw_b_c_vec
Bovar_2KmDispDist_output_tab[,"Patch_Area_imp"] = Bovar_2KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Bovar_2KmDispDist_HSI_tab = data.frame(imp = Bovar_2KmDispDist_HSI_vec)
Bovar_2KmDispDist_EgoSize_tab = data.frame(imp = Bovar_2KmDispDist_EgoSize_vec)
Bovar_2KmDispDist_strength_tab = data.frame(imp = Bovar_2KmDispDist_strength_vec)
Bovar_2KmDispDist_deg_tab = data.frame(imp = Bovar_2KmDispDist_deg_vec)
Bovar_2KmDispDist_habAv_tab = data.frame(imp = Bovar_2KmDispDist_habAv_vec)
Bovar_2KmDispDist_unw_b_c_tab = data.frame(imp = Bovar_2KmDispDist_unw_b_c_vec)
Bovar_2KmDispDist_Patch_Area_tab = data.frame(imp = Bovar_2KmDispDist_Patch_Area_vec)

Bovar_2KmDispDist_HSI_tab$variable = "HSI"
Bovar_2KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Bovar_2KmDispDist_strength_tab$variable = "Strength"
Bovar_2KmDispDist_deg_tab$variable = "Degree"
Bovar_2KmDispDist_habAv_tab$variable = "Hab. Av."
Bovar_2KmDispDist_unw_b_c_tab$variable = "B.C."
Bovar_2KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Bovar_2KmDispDist_AUC_cv_tab = data.frame(value = Bovar_2KmDispDist_AUC_cv_vec)
Bovar_2KmDispDist_AUC_train_tab = data.frame(value = Bovar_2KmDispDist_AUC_vec)
#Make label of model for plot
Bovar_2KmDispDist_AUC_cv_tab$model = "Bovar_2KmDispDist"
Bovar_2KmDispDist_AUC_train_tab$model = "Bovar_2KmDispDist"

#combine pred. vars. into new data frame 
Bovar_2KmDispDist_var_imp_tab = rbind(Bovar_2KmDispDist_HSI_tab,Bovar_2KmDispDist_EgoSize_tab,Bovar_2KmDispDist_strength_tab,Bovar_2KmDispDist_habAv_tab,Bovar_2KmDispDist_deg_tab,Bovar_2KmDispDist_unw_b_c_tab,Bovar_2KmDispDist_Patch_Area_tab)

ggplot(Bovar_2KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Bovar_2KmDispDist_var_imp_tab$imp~Bovar_2KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "2 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Bovar_2KmDispDist_HSI_vec)
mean(Bovar_2KmDispDist_EgoSize_vec)
mean(Bovar_2KmDispDist_strength_vec)
mean(Bovar_2KmDispDist_deg_vec)
mean(Bovar_2KmDispDist_habAv_vec)
mean(Bovar_2KmDispDist_unw_b_c_vec)
mean(Bovar_2KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Bovar_2KmDispDist_output_tab$AUC_cv)
summary(Bovar_2KmDispDist_output_tab$AUC_train)


#### 300m ######################################################################

# Create the output table that contains all the values of each of the runs
Bovar_300mDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar

# Create vectors to record performance measures
Bovar_300mDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Bovar_300mDispDist_AUC_vec = vector() #

Bovar_300mDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Bovar_300mDispDist_HSI_vec = vector() 
Bovar_300mDispDist_EgoSize_vec = vector() 
Bovar_300mDispDist_strength_vec  = vector()
Bovar_300mDispDist_deg_vec  = vector()
Bovar_300mDispDist_habAv_vec  = vector()
Bovar_300mDispDist_unw_b_c_vec  = vector()
Bovar_300mDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Bovar_300mDispDist_predict_mat = matrix(nrow = length(Bovar_300mDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Bovar_300mDispDist = gbm.step(data=Bovar_300mDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Bovar_300mDispDist)$var, imp = summary(gbm_mod_Bovar_300mDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Bovar_300mDispDist_HSI_vec = append(Bovar_300mDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Bovar_300mDispDist_EgoSize_vec = append(Bovar_300mDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Bovar_300mDispDist_strength_vec  = append(Bovar_300mDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Bovar_300mDispDist_deg_vec  = append(Bovar_300mDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Bovar_300mDispDist_habAv_vec  = append(Bovar_300mDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Bovar_300mDispDist_unw_b_c_vec  = append(Bovar_300mDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Bovar_300mDispDist_Patch_Area_vec  = append(Bovar_300mDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Bovar_300mDispDist_AUC_cv_vec = append(Bovar_300mDispDist_AUC_cv_vec, gbm_mod_Bovar_300mDispDist$cv.statistics$discrimination.mean)
  Bovar_300mDispDist_AUC_vec = append(Bovar_300mDispDist_AUC_vec, gbm_mod_Bovar_300mDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Bovar_300mDispDist$n.trees
  Bovar_300mDispDist_nt_vec = append(Bovar_300mDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Bovar_300mDispDist_predict_mat[,i] = predict(gbm_mod_Bovar_300mDispDist, Bovar_300mDispDist_stattest, gbm_mod_Bovar_300mDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Bovar_300mDispDist_output_tab[,"AUC_train"] = Bovar_300mDispDist_AUC_vec
Bovar_300mDispDist_output_tab[,"AUC_cv"] = Bovar_300mDispDist_AUC_cv_vec
Bovar_300mDispDist_output_tab[,"ntrees"] = Bovar_300mDispDist_nt_vec

#Var. importance columns
Bovar_300mDispDist_output_tab[,"HSI_imp"] = Bovar_300mDispDist_HSI_vec
Bovar_300mDispDist_output_tab[,"EgoSize_imp"] = Bovar_300mDispDist_EgoSize_vec
Bovar_300mDispDist_output_tab[,"strength_imp"] = Bovar_300mDispDist_strength_vec
Bovar_300mDispDist_output_tab[,"deg_imp"] = Bovar_300mDispDist_deg_vec
Bovar_300mDispDist_output_tab[,"habAv_imp"] = Bovar_300mDispDist_habAv_vec
Bovar_300mDispDist_output_tab[,"unw_b_c_imp"] = Bovar_300mDispDist_unw_b_c_vec
Bovar_300mDispDist_output_tab[,"Patch_Area_imp"] = Bovar_300mDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Bovar_300mDispDist_HSI_tab = data.frame(imp = Bovar_300mDispDist_HSI_vec)
Bovar_300mDispDist_EgoSize_tab = data.frame(imp = Bovar_300mDispDist_EgoSize_vec)
Bovar_300mDispDist_strength_tab = data.frame(imp = Bovar_300mDispDist_strength_vec)
Bovar_300mDispDist_deg_tab = data.frame(imp = Bovar_300mDispDist_deg_vec)
Bovar_300mDispDist_habAv_tab = data.frame(imp = Bovar_300mDispDist_habAv_vec)
Bovar_300mDispDist_unw_b_c_tab = data.frame(imp = Bovar_300mDispDist_unw_b_c_vec)
Bovar_300mDispDist_Patch_Area_tab = data.frame(imp = Bovar_300mDispDist_Patch_Area_vec)

Bovar_300mDispDist_HSI_tab$variable = "HSI"
Bovar_300mDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Bovar_300mDispDist_strength_tab$variable = "Strength"
Bovar_300mDispDist_deg_tab$variable = "Degree"
Bovar_300mDispDist_habAv_tab$variable = "Hab. Av."
Bovar_300mDispDist_unw_b_c_tab$variable = "B.C."
Bovar_300mDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Bovar_300mDispDist_AUC_cv_tab = data.frame(value = Bovar_300mDispDist_AUC_cv_vec)
Bovar_300mDispDist_AUC_train_tab = data.frame(value = Bovar_300mDispDist_AUC_vec)
#Make label of model for plot
Bovar_300mDispDist_AUC_cv_tab$model = "Bovar_300mDispDist"
Bovar_300mDispDist_AUC_train_tab$model = "Bovar_300mDispDist"

#combine pred. vars. into new data frame 
Bovar_300mDispDist_var_imp_tab = rbind(Bovar_300mDispDist_HSI_tab,Bovar_300mDispDist_EgoSize_tab,Bovar_300mDispDist_strength_tab,Bovar_300mDispDist_habAv_tab,Bovar_300mDispDist_deg_tab,Bovar_300mDispDist_unw_b_c_tab,Bovar_300mDispDist_Patch_Area_tab)

ggplot(Bovar_300mDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Bovar_300mDispDist_var_imp_tab$imp~Bovar_300mDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "300 m maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Bovar_300mDispDist_HSI_vec)
mean(Bovar_300mDispDist_EgoSize_vec)
mean(Bovar_300mDispDist_strength_vec)
mean(Bovar_300mDispDist_deg_vec)
mean(Bovar_300mDispDist_habAv_vec)
mean(Bovar_300mDispDist_unw_b_c_vec)
mean(Bovar_300mDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Bovar_300mDispDist_output_tab$AUC_cv)
summary(Bovar_300mDispDist_output_tab$AUC_train)


#### 10Km ######################################################################

# Create the output table that contains all the values of each of the runs
Bovar_10KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar

# Create vectors to record performance measures
Bovar_10KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Bovar_10KmDispDist_AUC_vec = vector() #

Bovar_10KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Bovar_10KmDispDist_HSI_vec = vector() 
Bovar_10KmDispDist_EgoSize_vec = vector() 
Bovar_10KmDispDist_strength_vec  = vector()
Bovar_10KmDispDist_deg_vec  = vector()
Bovar_10KmDispDist_habAv_vec  = vector()
Bovar_10KmDispDist_unw_b_c_vec  = vector()
Bovar_10KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Bovar_10KmDispDist_predict_mat = matrix(nrow = length(Bovar_10KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Bovar_10KmDispDist = gbm.step(data=Bovar_10KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Bovar_10KmDispDist)$var, imp = summary(gbm_mod_Bovar_10KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Bovar_10KmDispDist_HSI_vec = append(Bovar_10KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Bovar_10KmDispDist_EgoSize_vec = append(Bovar_10KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Bovar_10KmDispDist_strength_vec  = append(Bovar_10KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Bovar_10KmDispDist_deg_vec  = append(Bovar_10KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Bovar_10KmDispDist_habAv_vec  = append(Bovar_10KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Bovar_10KmDispDist_unw_b_c_vec  = append(Bovar_10KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Bovar_10KmDispDist_Patch_Area_vec  = append(Bovar_10KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Bovar_10KmDispDist_AUC_cv_vec = append(Bovar_10KmDispDist_AUC_cv_vec, gbm_mod_Bovar_10KmDispDist$cv.statistics$discrimination.mean)
  Bovar_10KmDispDist_AUC_vec = append(Bovar_10KmDispDist_AUC_vec, gbm_mod_Bovar_10KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Bovar_10KmDispDist$n.trees
  Bovar_10KmDispDist_nt_vec = append(Bovar_10KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Bovar_10KmDispDist_predict_mat[,i] = predict(gbm_mod_Bovar_10KmDispDist, Bovar_10KmDispDist_stattest, gbm_mod_Bovar_10KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Bovar_10KmDispDist_output_tab[,"AUC_train"] = Bovar_10KmDispDist_AUC_vec
Bovar_10KmDispDist_output_tab[,"AUC_cv"] = Bovar_10KmDispDist_AUC_cv_vec
Bovar_10KmDispDist_output_tab[,"ntrees"] = Bovar_10KmDispDist_nt_vec

#Var. importance columns
Bovar_10KmDispDist_output_tab[,"HSI_imp"] = Bovar_10KmDispDist_HSI_vec
Bovar_10KmDispDist_output_tab[,"EgoSize_imp"] = Bovar_10KmDispDist_EgoSize_vec
Bovar_10KmDispDist_output_tab[,"strength_imp"] = Bovar_10KmDispDist_strength_vec
Bovar_10KmDispDist_output_tab[,"deg_imp"] = Bovar_10KmDispDist_deg_vec
Bovar_10KmDispDist_output_tab[,"habAv_imp"] = Bovar_10KmDispDist_habAv_vec
Bovar_10KmDispDist_output_tab[,"unw_b_c_imp"] = Bovar_10KmDispDist_unw_b_c_vec
Bovar_10KmDispDist_output_tab[,"Patch_Area_imp"] = Bovar_10KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Bovar_10KmDispDist_HSI_tab = data.frame(imp = Bovar_10KmDispDist_HSI_vec)
Bovar_10KmDispDist_EgoSize_tab = data.frame(imp = Bovar_10KmDispDist_EgoSize_vec)
Bovar_10KmDispDist_strength_tab = data.frame(imp = Bovar_10KmDispDist_strength_vec)
Bovar_10KmDispDist_deg_tab = data.frame(imp = Bovar_10KmDispDist_deg_vec)
Bovar_10KmDispDist_habAv_tab = data.frame(imp = Bovar_10KmDispDist_habAv_vec)
Bovar_10KmDispDist_unw_b_c_tab = data.frame(imp = Bovar_10KmDispDist_unw_b_c_vec)
Bovar_10KmDispDist_Patch_Area_tab = data.frame(imp = Bovar_10KmDispDist_Patch_Area_vec)

Bovar_10KmDispDist_HSI_tab$variable = "HSI"
Bovar_10KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Bovar_10KmDispDist_strength_tab$variable = "Strength"
Bovar_10KmDispDist_deg_tab$variable = "Degree"
Bovar_10KmDispDist_habAv_tab$variable = "Hab. Av."
Bovar_10KmDispDist_unw_b_c_tab$variable = "B.C."
Bovar_10KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Bovar_10KmDispDist_AUC_cv_tab = data.frame(value = Bovar_10KmDispDist_AUC_cv_vec)
Bovar_10KmDispDist_AUC_train_tab = data.frame(value = Bovar_10KmDispDist_AUC_vec)
#Make label of model for plot
Bovar_10KmDispDist_AUC_cv_tab$model = "Bovar_10KmDispDist"
Bovar_10KmDispDist_AUC_train_tab$model = "Bovar_10KmDispDist"

#combine pred. vars. into new data frame 
Bovar_10KmDispDist_var_imp_tab = rbind(Bovar_10KmDispDist_HSI_tab,Bovar_10KmDispDist_EgoSize_tab,Bovar_10KmDispDist_strength_tab,Bovar_10KmDispDist_habAv_tab,Bovar_10KmDispDist_deg_tab,Bovar_10KmDispDist_unw_b_c_tab,Bovar_10KmDispDist_Patch_Area_tab)

ggplot(Bovar_10KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Bovar_10KmDispDist_var_imp_tab$imp~Bovar_10KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "10 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Bovar_10KmDispDist_HSI_vec)
mean(Bovar_10KmDispDist_EgoSize_vec)
mean(Bovar_10KmDispDist_strength_vec)
mean(Bovar_10KmDispDist_deg_vec)
mean(Bovar_10KmDispDist_habAv_vec)
mean(Bovar_10KmDispDist_unw_b_c_vec)
mean(Bovar_10KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Bovar_10KmDispDist_output_tab$AUC_cv)
summary(Bovar_10KmDispDist_output_tab$AUC_train)


####################################### 1Km ####################################

# Create the output table that contains all the values of each of the runs
Bovar_1KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar

# Create vectors to record performance measures
Bovar_1KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Bovar_1KmDispDist_AUC_vec = vector() #

Bovar_1KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Bovar_1KmDispDist_HSI_vec = vector() 
Bovar_1KmDispDist_EgoSize_vec = vector() 
Bovar_1KmDispDist_strength_vec  = vector()
Bovar_1KmDispDist_deg_vec  = vector()
Bovar_1KmDispDist_habAv_vec  = vector()
Bovar_1KmDispDist_unw_b_c_vec  = vector()
Bovar_1KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Bovar_1KmDispDist_predict_mat = matrix(nrow = length(Bovar_1KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Bovar_1KmDispDist = gbm.step(data=Bovar_1KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Bovar_1KmDispDist)$var, imp = summary(gbm_mod_Bovar_1KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Bovar_1KmDispDist_HSI_vec = append(Bovar_1KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Bovar_1KmDispDist_EgoSize_vec = append(Bovar_1KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Bovar_1KmDispDist_strength_vec  = append(Bovar_1KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Bovar_1KmDispDist_deg_vec  = append(Bovar_1KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Bovar_1KmDispDist_habAv_vec  = append(Bovar_1KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Bovar_1KmDispDist_unw_b_c_vec  = append(Bovar_1KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Bovar_1KmDispDist_Patch_Area_vec  = append(Bovar_1KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Bovar_1KmDispDist_AUC_cv_vec = append(Bovar_1KmDispDist_AUC_cv_vec, gbm_mod_Bovar_1KmDispDist$cv.statistics$discrimination.mean)
  Bovar_1KmDispDist_AUC_vec = append(Bovar_1KmDispDist_AUC_vec, gbm_mod_Bovar_1KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Bovar_1KmDispDist$n.trees
  Bovar_1KmDispDist_nt_vec = append(Bovar_1KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Bovar_1KmDispDist_predict_mat[,i] = predict(gbm_mod_Bovar_1KmDispDist, Bovar_1KmDispDist_stattest, gbm_mod_Bovar_1KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Bovar_1KmDispDist_output_tab[,"AUC_train"] = Bovar_1KmDispDist_AUC_vec
Bovar_1KmDispDist_output_tab[,"AUC_cv"] = Bovar_1KmDispDist_AUC_cv_vec
Bovar_1KmDispDist_output_tab[,"ntrees"] = Bovar_1KmDispDist_nt_vec

#Var. importance columns
Bovar_1KmDispDist_output_tab[,"HSI_imp"] = Bovar_1KmDispDist_HSI_vec
Bovar_1KmDispDist_output_tab[,"EgoSize_imp"] = Bovar_1KmDispDist_EgoSize_vec
Bovar_1KmDispDist_output_tab[,"strength_imp"] = Bovar_1KmDispDist_strength_vec
Bovar_1KmDispDist_output_tab[,"deg_imp"] = Bovar_1KmDispDist_deg_vec
Bovar_1KmDispDist_output_tab[,"habAv_imp"] = Bovar_1KmDispDist_habAv_vec
Bovar_1KmDispDist_output_tab[,"unw_b_c_imp"] = Bovar_1KmDispDist_unw_b_c_vec
Bovar_1KmDispDist_output_tab[,"Patch_Area_imp"] = Bovar_1KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Bovar_1KmDispDist_HSI_tab = data.frame(imp = Bovar_1KmDispDist_HSI_vec)
Bovar_1KmDispDist_EgoSize_tab = data.frame(imp = Bovar_1KmDispDist_EgoSize_vec)
Bovar_1KmDispDist_strength_tab = data.frame(imp = Bovar_1KmDispDist_strength_vec)
Bovar_1KmDispDist_deg_tab = data.frame(imp = Bovar_1KmDispDist_deg_vec)
Bovar_1KmDispDist_habAv_tab = data.frame(imp = Bovar_1KmDispDist_habAv_vec)
Bovar_1KmDispDist_unw_b_c_tab = data.frame(imp = Bovar_1KmDispDist_unw_b_c_vec)
Bovar_1KmDispDist_Patch_Area_tab = data.frame(imp = Bovar_1KmDispDist_Patch_Area_vec)

Bovar_1KmDispDist_HSI_tab$variable = "HSI"
Bovar_1KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Bovar_1KmDispDist_strength_tab$variable = "Strength"
Bovar_1KmDispDist_deg_tab$variable = "Degree"
Bovar_1KmDispDist_habAv_tab$variable = "Hab. Av."
Bovar_1KmDispDist_unw_b_c_tab$variable = "B.C."
Bovar_1KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Bovar_1KmDispDist_AUC_cv_tab = data.frame(value = Bovar_1KmDispDist_AUC_cv_vec)
Bovar_1KmDispDist_AUC_train_tab = data.frame(value = Bovar_1KmDispDist_AUC_vec)
#Make label of model for plot
Bovar_1KmDispDist_AUC_cv_tab$model = "Bovar_1KmDispDist"
Bovar_1KmDispDist_AUC_train_tab$model = "Bovar_1KmDispDist"

#combine pred. vars. into new data frame 
Bovar_1KmDispDist_var_imp_tab = rbind(Bovar_1KmDispDist_HSI_tab,Bovar_1KmDispDist_EgoSize_tab,Bovar_1KmDispDist_strength_tab,Bovar_1KmDispDist_habAv_tab,Bovar_1KmDispDist_deg_tab,Bovar_1KmDispDist_unw_b_c_tab,Bovar_1KmDispDist_Patch_Area_tab)

ggplot(Bovar_1KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Bovar_1KmDispDist_var_imp_tab$imp~Bovar_1KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "1 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Bovar_1KmDispDist_HSI_vec)
mean(Bovar_1KmDispDist_EgoSize_vec)
mean(Bovar_1KmDispDist_strength_vec)
mean(Bovar_1KmDispDist_deg_vec)
mean(Bovar_1KmDispDist_habAv_vec)
mean(Bovar_1KmDispDist_unw_b_c_vec)
mean(Bovar_1KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Bovar_1KmDispDist_output_tab$AUC_cv)
summary(Bovar_1KmDispDist_output_tab$AUC_train)


####################################### 6Km ####################################

# Create the output table that contains all the values of each of the runs
Bovar_6KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar

# Create vectors to record performance measures
Bovar_6KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Bovar_6KmDispDist_AUC_vec = vector() #

Bovar_6KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Bovar_6KmDispDist_HSI_vec = vector() 
Bovar_6KmDispDist_EgoSize_vec = vector() 
Bovar_6KmDispDist_strength_vec  = vector()
Bovar_6KmDispDist_deg_vec  = vector()
Bovar_6KmDispDist_habAv_vec  = vector()
Bovar_6KmDispDist_unw_b_c_vec  = vector()
Bovar_6KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Bovar_6KmDispDist_predict_mat = matrix(nrow = length(Bovar_6KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Bovar_6KmDispDist = gbm.step(data=Bovar_6KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Bovar_6KmDispDist)$var, imp = summary(gbm_mod_Bovar_6KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Bovar_6KmDispDist_HSI_vec = append(Bovar_6KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Bovar_6KmDispDist_EgoSize_vec = append(Bovar_6KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Bovar_6KmDispDist_strength_vec  = append(Bovar_6KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Bovar_6KmDispDist_deg_vec  = append(Bovar_6KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Bovar_6KmDispDist_habAv_vec  = append(Bovar_6KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Bovar_6KmDispDist_unw_b_c_vec  = append(Bovar_6KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Bovar_6KmDispDist_Patch_Area_vec  = append(Bovar_6KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Bovar_6KmDispDist_AUC_cv_vec = append(Bovar_6KmDispDist_AUC_cv_vec, gbm_mod_Bovar_6KmDispDist$cv.statistics$discrimination.mean)
  Bovar_6KmDispDist_AUC_vec = append(Bovar_6KmDispDist_AUC_vec, gbm_mod_Bovar_6KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Bovar_6KmDispDist$n.trees
  Bovar_6KmDispDist_nt_vec = append(Bovar_6KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Bovar_6KmDispDist_predict_mat[,i] = predict(gbm_mod_Bovar_6KmDispDist, Bovar_6KmDispDist_stattest, gbm_mod_Bovar_6KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Bovar_6KmDispDist_output_tab[,"AUC_train"] = Bovar_6KmDispDist_AUC_vec
Bovar_6KmDispDist_output_tab[,"AUC_cv"] = Bovar_6KmDispDist_AUC_cv_vec
Bovar_6KmDispDist_output_tab[,"ntrees"] = Bovar_6KmDispDist_nt_vec

#Var. importance columns
Bovar_6KmDispDist_output_tab[,"HSI_imp"] = Bovar_6KmDispDist_HSI_vec
Bovar_6KmDispDist_output_tab[,"EgoSize_imp"] = Bovar_6KmDispDist_EgoSize_vec
Bovar_6KmDispDist_output_tab[,"strength_imp"] = Bovar_6KmDispDist_strength_vec
Bovar_6KmDispDist_output_tab[,"deg_imp"] = Bovar_6KmDispDist_deg_vec
Bovar_6KmDispDist_output_tab[,"habAv_imp"] = Bovar_6KmDispDist_habAv_vec
Bovar_6KmDispDist_output_tab[,"unw_b_c_imp"] = Bovar_6KmDispDist_unw_b_c_vec
Bovar_6KmDispDist_output_tab[,"Patch_Area_imp"] = Bovar_6KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Bovar_6KmDispDist_HSI_tab = data.frame(imp = Bovar_6KmDispDist_HSI_vec)
Bovar_6KmDispDist_EgoSize_tab = data.frame(imp = Bovar_6KmDispDist_EgoSize_vec)
Bovar_6KmDispDist_strength_tab = data.frame(imp = Bovar_6KmDispDist_strength_vec)
Bovar_6KmDispDist_deg_tab = data.frame(imp = Bovar_6KmDispDist_deg_vec)
Bovar_6KmDispDist_habAv_tab = data.frame(imp = Bovar_6KmDispDist_habAv_vec)
Bovar_6KmDispDist_unw_b_c_tab = data.frame(imp = Bovar_6KmDispDist_unw_b_c_vec)
Bovar_6KmDispDist_Patch_Area_tab = data.frame(imp = Bovar_6KmDispDist_Patch_Area_vec)

Bovar_6KmDispDist_HSI_tab$variable = "HSI"
Bovar_6KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Bovar_6KmDispDist_strength_tab$variable = "Strength"
Bovar_6KmDispDist_deg_tab$variable = "Degree"
Bovar_6KmDispDist_habAv_tab$variable = "Hab. Av."
Bovar_6KmDispDist_unw_b_c_tab$variable = "B.C."
Bovar_6KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Bovar_6KmDispDist_AUC_cv_tab = data.frame(value = Bovar_6KmDispDist_AUC_cv_vec)
Bovar_6KmDispDist_AUC_train_tab = data.frame(value = Bovar_6KmDispDist_AUC_vec)
#Make label of model for plot
Bovar_6KmDispDist_AUC_cv_tab$model = "Bovar_6KmDispDist"
Bovar_6KmDispDist_AUC_train_tab$model = "Bovar_6KmDispDist"

#combine pred. vars. into new data frame 
Bovar_6KmDispDist_var_imp_tab = rbind(Bovar_6KmDispDist_HSI_tab,Bovar_6KmDispDist_EgoSize_tab,Bovar_6KmDispDist_strength_tab,Bovar_6KmDispDist_habAv_tab,Bovar_6KmDispDist_deg_tab,Bovar_6KmDispDist_unw_b_c_tab,Bovar_6KmDispDist_Patch_Area_tab)

ggplot(Bovar_6KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Bovar_6KmDispDist_var_imp_tab$imp~Bovar_6KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "6 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Bovar_6KmDispDist_HSI_vec)
mean(Bovar_6KmDispDist_EgoSize_vec)
mean(Bovar_6KmDispDist_strength_vec)
mean(Bovar_6KmDispDist_deg_vec)
mean(Bovar_6KmDispDist_habAv_vec)
mean(Bovar_6KmDispDist_unw_b_c_vec)
mean(Bovar_6KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Bovar_6KmDispDist_output_tab$AUC_cv)
summary(Bovar_6KmDispDist_output_tab$AUC_train)


####################################### 4Km
# Create the output table that contains all the values of each of the runs
# Bovar_4KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar
# 
# # Create vectors to record performance measures
# Bovar_4KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
# Bovar_4KmDispDist_AUC_vec = vector() #
# 
# Bovar_4KmDispDist_nt_vec = vector() #Number of trees
# 
# #Create vectors with all the values of a certain predictor along all the runs
# Bovar_4KmDispDist_HSI_vec = vector() 
# Bovar_4KmDispDist_EgoSize_vec = vector() 
# Bovar_4KmDispDist_strength_vec  = vector()
# Bovar_4KmDispDist_deg_vec  = vector()
# Bovar_4KmDispDist_habAv_vec  = vector()
# Bovar_4KmDispDist_unw_b_c_vec  = vector()
# Bovar_4KmDispDist_Patch_Area_vec  = vector()
# 
# #Discrete Prediction of occurrence state for all patches
# Bovar_4KmDispDist_predict_mat = matrix(nrow = length(Bovar_4KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)
# 
# ### Loop that goes exactly for 100 iterations, to get distributions 
# for(i in c(1:n_repeats)){
#   
#   #Perform gbm step to set number of trees, no cross-validation.
#   gbm_mod_Bovar_4KmDispDist = gbm.step(data=Bovar_4KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
#   
#   #data frame of variable importance, to fill the vectors of the model var. importance scores
#   var_imp = data.frame(var = summary(gbm_mod_Bovar_4KmDispDist)$var, imp = summary(gbm_mod_Bovar_4KmDispDist)$rel.inf)
#   
#   #Write var_imp results to the vectors
#   Bovar_4KmDispDist_HSI_vec = append(Bovar_4KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
#   Bovar_4KmDispDist_EgoSize_vec = append(Bovar_4KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
#   Bovar_4KmDispDist_strength_vec  = append(Bovar_4KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
#   Bovar_4KmDispDist_deg_vec  = append(Bovar_4KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
#   Bovar_4KmDispDist_habAv_vec  = append(Bovar_4KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
#   Bovar_4KmDispDist_unw_b_c_vec  = append(Bovar_4KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
#   Bovar_4KmDispDist_Patch_Area_vec  = append(Bovar_4KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
#   
#   #Write the AUC & CV-AUC of this model to a vector.
#   Bovar_4KmDispDist_AUC_cv_vec = append(Bovar_4KmDispDist_AUC_cv_vec, gbm_mod_Bovar_4KmDispDist$cv.statistics$discrimination.mean)
#   Bovar_4KmDispDist_AUC_vec = append(Bovar_4KmDispDist_AUC_vec, gbm_mod_Bovar_4KmDispDist$self.statistics$discrimination)
#   
#   #Write the number of trees to a vector
#   nt = gbm_mod_Bovar_4KmDispDist$n.trees
#   Bovar_4KmDispDist_nt_vec = append(Bovar_4KmDispDist_nt_vec, nt)
#   print(nt)
#   
#   #write the continuous prediction over all the patches
#   Bovar_4KmDispDist_predict_mat[,i] = predict(gbm_mod_Bovar_4KmDispDist, Bovar_4KmDispDist_stattest, gbm_mod_Bovar_4KmDispDist$n.trees, type = "response", single.tree = FALSE)
#   
#   print(paste("Finished:",i,"/",n_repeats,sep = ""))
# }
# 
# Bovar_4KmDispDist_output_tab[,"AUC_train"] = Bovar_4KmDispDist_AUC_vec
# Bovar_4KmDispDist_output_tab[,"AUC_cv"] = Bovar_4KmDispDist_AUC_cv_vec
# Bovar_4KmDispDist_output_tab[,"ntrees"] = Bovar_4KmDispDist_nt_vec
# 
# #Var. importance columns
# Bovar_4KmDispDist_output_tab[,"HSI_imp"] = Bovar_4KmDispDist_HSI_vec
# Bovar_4KmDispDist_output_tab[,"EgoSize_imp"] = Bovar_4KmDispDist_EgoSize_vec
# Bovar_4KmDispDist_output_tab[,"strength_imp"] = Bovar_4KmDispDist_strength_vec
# Bovar_4KmDispDist_output_tab[,"deg_imp"] = Bovar_4KmDispDist_deg_vec
# Bovar_4KmDispDist_output_tab[,"habAv_imp"] = Bovar_4KmDispDist_habAv_vec
# Bovar_4KmDispDist_output_tab[,"unw_b_c_imp"] = Bovar_4KmDispDist_unw_b_c_vec
# Bovar_4KmDispDist_output_tab[,"Patch_Area_imp"] = Bovar_4KmDispDist_Patch_Area_vec
# #is.data.frame(output_tab)
# 
# # Make a dataframe for plotting overlaying histrograms in R
# Bovar_4KmDispDist_HSI_tab = data.frame(imp = Bovar_4KmDispDist_HSI_vec)
# Bovar_4KmDispDist_EgoSize_tab = data.frame(imp = Bovar_4KmDispDist_EgoSize_vec)
# Bovar_4KmDispDist_strength_tab = data.frame(imp = Bovar_4KmDispDist_strength_vec)
# Bovar_4KmDispDist_deg_tab = data.frame(imp = Bovar_4KmDispDist_deg_vec)
# Bovar_4KmDispDist_habAv_tab = data.frame(imp = Bovar_4KmDispDist_habAv_vec)
# Bovar_4KmDispDist_unw_b_c_tab = data.frame(imp = Bovar_4KmDispDist_unw_b_c_vec)
# Bovar_4KmDispDist_Patch_Area_tab = data.frame(imp = Bovar_4KmDispDist_Patch_Area_vec)
# 
# Bovar_4KmDispDist_HSI_tab$variable = "HSI"
# Bovar_4KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
# Bovar_4KmDispDist_strength_tab$variable = "Strength"
# Bovar_4KmDispDist_deg_tab$variable = "Degree"
# Bovar_4KmDispDist_habAv_tab$variable = "Hab. Av."
# Bovar_4KmDispDist_unw_b_c_tab$variable = "B.C."
# Bovar_4KmDispDist_Patch_Area_tab$variable = "Patch Area"
# 
# #Reserve also measures in df to do overlaying histograms comparing performance between models
# Bovar_4KmDispDist_AUC_cv_tab = data.frame(value = Bovar_4KmDispDist_AUC_cv_vec)
# Bovar_4KmDispDist_AUC_train_tab = data.frame(value = Bovar_4KmDispDist_AUC_vec)
# #Make label of model for plot
# Bovar_4KmDispDist_AUC_cv_tab$model = "Bovar_4KmDispDist"
# Bovar_4KmDispDist_AUC_train_tab$model = "Bovar_4KmDispDist"
# 
# #combine pred. vars. into new data frame 
# Bovar_4KmDispDist_var_imp_tab = rbind(Bovar_4KmDispDist_HSI_tab,Bovar_4KmDispDist_EgoSize_tab,Bovar_4KmDispDist_strength_tab,Bovar_4KmDispDist_habAv_tab,Bovar_4KmDispDist_deg_tab,Bovar_4KmDispDist_unw_b_c_tab,Bovar_4KmDispDist_Patch_Area_tab)
# 
# ggplot(Bovar_4KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
# boxplot(Bovar_4KmDispDist_var_imp_tab$imp~Bovar_4KmDispDist_var_imp_tab$variable,
# xlab = NULL, ylab= "Variable importance", 
# main = "4 km maximum dispersal distance",
# cex.axis = 1.25, cex.lab = 1.2)
# 
# #Get mean var. importance of all of the vars. 
# mean(Bovar_4KmDispDist_HSI_vec)
# mean(Bovar_4KmDispDist_EgoSize_vec)
# mean(Bovar_4KmDispDist_strength_vec)
# mean(Bovar_4KmDispDist_deg_vec)
# mean(Bovar_4KmDispDist_habAv_vec)
# mean(Bovar_4KmDispDist_unw_b_c_vec)
# mean(Bovar_4KmDispDist_Patch_Area_vec)
# 
# ##Check distr. of measures of prediction accuracy 
# summary(Bovar_4KmDispDist_output_tab$AUC_cv)
# summary(Bovar_4KmDispDist_output_tab$AUC_train)


####################################### 8Km 
# Create the output table that contains all the values of each of the runs
Bovar_8KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar

# Create vectors to record performance measures
Bovar_8KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Bovar_8KmDispDist_AUC_vec = vector() #

Bovar_8KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Bovar_8KmDispDist_HSI_vec = vector() 
Bovar_8KmDispDist_EgoSize_vec = vector() 
Bovar_8KmDispDist_strength_vec  = vector()
Bovar_8KmDispDist_deg_vec  = vector()
Bovar_8KmDispDist_habAv_vec  = vector()
Bovar_8KmDispDist_unw_b_c_vec  = vector()
Bovar_8KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Bovar_8KmDispDist_predict_mat = matrix(nrow = length(Bovar_8KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Bovar_8KmDispDist = gbm.step(data=Bovar_8KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Bovar_8KmDispDist)$var, imp = summary(gbm_mod_Bovar_8KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Bovar_8KmDispDist_HSI_vec = append(Bovar_8KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Bovar_8KmDispDist_EgoSize_vec = append(Bovar_8KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Bovar_8KmDispDist_strength_vec  = append(Bovar_8KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Bovar_8KmDispDist_deg_vec  = append(Bovar_8KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Bovar_8KmDispDist_habAv_vec  = append(Bovar_8KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Bovar_8KmDispDist_unw_b_c_vec  = append(Bovar_8KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Bovar_8KmDispDist_Patch_Area_vec  = append(Bovar_8KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Bovar_8KmDispDist_AUC_cv_vec = append(Bovar_8KmDispDist_AUC_cv_vec, gbm_mod_Bovar_8KmDispDist$cv.statistics$discrimination.mean)
  Bovar_8KmDispDist_AUC_vec = append(Bovar_8KmDispDist_AUC_vec, gbm_mod_Bovar_8KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Bovar_8KmDispDist$n.trees
  Bovar_8KmDispDist_nt_vec = append(Bovar_8KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Bovar_8KmDispDist_predict_mat[,i] = predict(gbm_mod_Bovar_8KmDispDist, Bovar_8KmDispDist_stattest, gbm_mod_Bovar_8KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Bovar_8KmDispDist_output_tab[,"AUC_train"] = Bovar_8KmDispDist_AUC_vec
Bovar_8KmDispDist_output_tab[,"AUC_cv"] = Bovar_8KmDispDist_AUC_cv_vec
Bovar_8KmDispDist_output_tab[,"ntrees"] = Bovar_8KmDispDist_nt_vec

#Var. importance columns
Bovar_8KmDispDist_output_tab[,"HSI_imp"] = Bovar_8KmDispDist_HSI_vec
Bovar_8KmDispDist_output_tab[,"EgoSize_imp"] = Bovar_8KmDispDist_EgoSize_vec
Bovar_8KmDispDist_output_tab[,"strength_imp"] = Bovar_8KmDispDist_strength_vec
Bovar_8KmDispDist_output_tab[,"deg_imp"] = Bovar_8KmDispDist_deg_vec
Bovar_8KmDispDist_output_tab[,"habAv_imp"] = Bovar_8KmDispDist_habAv_vec
Bovar_8KmDispDist_output_tab[,"unw_b_c_imp"] = Bovar_8KmDispDist_unw_b_c_vec
Bovar_8KmDispDist_output_tab[,"Patch_Area_imp"] = Bovar_8KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Bovar_8KmDispDist_HSI_tab = data.frame(imp = Bovar_8KmDispDist_HSI_vec)
Bovar_8KmDispDist_EgoSize_tab = data.frame(imp = Bovar_8KmDispDist_EgoSize_vec)
Bovar_8KmDispDist_strength_tab = data.frame(imp = Bovar_8KmDispDist_strength_vec)
Bovar_8KmDispDist_deg_tab = data.frame(imp = Bovar_8KmDispDist_deg_vec)
Bovar_8KmDispDist_habAv_tab = data.frame(imp = Bovar_8KmDispDist_habAv_vec)
Bovar_8KmDispDist_unw_b_c_tab = data.frame(imp = Bovar_8KmDispDist_unw_b_c_vec)
Bovar_8KmDispDist_Patch_Area_tab = data.frame(imp = Bovar_8KmDispDist_Patch_Area_vec)

Bovar_8KmDispDist_HSI_tab$variable = "HSI"
Bovar_8KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Bovar_8KmDispDist_strength_tab$variable = "Strength"
Bovar_8KmDispDist_deg_tab$variable = "Degree"
Bovar_8KmDispDist_habAv_tab$variable = "Hab. Av."
Bovar_8KmDispDist_unw_b_c_tab$variable = "B.C."
Bovar_8KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Bovar_8KmDispDist_AUC_cv_tab = data.frame(value = Bovar_8KmDispDist_AUC_cv_vec)
Bovar_8KmDispDist_AUC_train_tab = data.frame(value = Bovar_8KmDispDist_AUC_vec)
#Make label of model for plot
Bovar_8KmDispDist_AUC_cv_tab$model = "Bovar_8KmDispDist"
Bovar_8KmDispDist_AUC_train_tab$model = "Bovar_8KmDispDist"

#combine pred. vars. into new data frame 
Bovar_8KmDispDist_var_imp_tab = rbind(Bovar_8KmDispDist_HSI_tab,Bovar_8KmDispDist_EgoSize_tab,Bovar_8KmDispDist_strength_tab,Bovar_8KmDispDist_habAv_tab,Bovar_8KmDispDist_deg_tab,Bovar_8KmDispDist_unw_b_c_tab,Bovar_8KmDispDist_Patch_Area_tab)

ggplot(Bovar_8KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Bovar_8KmDispDist_var_imp_tab$imp~Bovar_8KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "8 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Bovar_8KmDispDist_HSI_vec)
mean(Bovar_8KmDispDist_EgoSize_vec)
mean(Bovar_8KmDispDist_strength_vec)
mean(Bovar_8KmDispDist_deg_vec)
mean(Bovar_8KmDispDist_habAv_vec)
mean(Bovar_8KmDispDist_unw_b_c_vec)
mean(Bovar_8KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Bovar_8KmDispDist_output_tab$AUC_cv)
summary(Bovar_8KmDispDist_output_tab$AUC_train)


####################################################################################################
#### Compare scores between networks w/d0 variations of the same species ###########################
### Bovar ##################
setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
#Import evaluation df's of original run with species-specific dispersal distance
Bovar_DefaultDispDist_AUC_cv_tab <- read.csv("Bovar_DefaultDispDist_AUC_cv_tab.csv")
Bovar_DefaultDispDist_AUC_train_tab <- read.csv("Bovar_DefaultDispDist_AUC_train_tab.csv")
Bovar_DefaultDispDist_var_imp_tab <- read.csv("Bovar_DefaultDispDist_var_imp_tab.csv")
Bovar_DefaultDispDist_output_tab <- read.csv("Bovar_DefaultDispDist_BRToutput_tab.csv")
Bovar_DefaultDispDist_predict_df <- read.csv("Bovar_DefaultDispDist_predict_df.csv")

head(Bovar_DefaultDispDist_AUC_cv_tab)
Bovar_DefaultDispDist_AUC_cv_tab$X <- NULL
Bovar_DefaultDispDist_AUC_train_tab$X <- NULL

#### Make a dataframe for plotting overlaying histograms in R
###cv AUC
Bovar_AUC_cv_tab = rbind(Bovar_DefaultDispDist_AUC_cv_tab, Bovar_300mDispDist_AUC_cv_tab, Bovar_1KmDispDist_AUC_cv_tab, Bovar_2KmDispDist_AUC_cv_tab, 
                         # Bovar_4KmDispDist_AUC_cv_tab, 
                         Bovar_6KmDispDist_AUC_cv_tab, Bovar_8KmDispDist_AUC_cv_tab, Bovar_10KmDispDist_AUC_cv_tab)
#Change order of factors 
Bovar_AUC_cv_tab$model<- as.factor(Bovar_AUC_cv_tab$model)
levels(Bovar_AUC_cv_tab$model)
Bovar_AUC_cv_tab$model<-factor(Bovar_AUC_cv_tab$model, levels=c("Bovar_300mDispDist", "Bovar_1KmDispDist", "Bovar_2KmDispDist", 
                                                                "Bovar", "Bovar_6KmDispDist", "Bovar_8KmDispDist", "Bovar_10KmDispDist"))

#Plot
ggplot(Bovar_AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Bovar_AUC_cv_tab$value~Bovar_AUC_cv_tab$model, ylab= "Cross-validated AUC")

### train AUC
Bovar_AUC_train_tab = rbind(Bovar_DefaultDispDist_AUC_train_tab, Bovar_300mDispDist_AUC_train_tab, Bovar_1KmDispDist_AUC_train_tab, Bovar_2KmDispDist_AUC_train_tab, 
                            #Bovar_4KmDispDist_AUC_train_tab, 
                            Bovar_6KmDispDist_AUC_train_tab, Bovar_8KmDispDist_AUC_train_tab, Bovar_10KmDispDist_AUC_train_tab)
#Change order of factors 
levels(Bovar_AUC_train_tab$model)
Bovar_AUC_train_tab$model<-factor(Bovar_AUC_train_tab$model, levels=c("Bovar_300mDispDist", "Bovar_1KmDispDist", "Bovar_2KmDispDist", 
                                                                      "Bovar", "Bovar_6KmDispDist", "Bovar_8KmDispDist", "Bovar_10KmDispDist"))
#Plot
ggplot(Bovar_AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Bovar_AUC_train_tab$value~Bovar_AUC_train_tab$model, ylab= "Training AUC")

summary(Bovar_DefaultDispDist_output_tab$AUC_cv)
summary(Bovar_300mDispDist_output_tab$AUC_cv)
summary(Bovar_1KmDispDist_output_tab$AUC_cv)
summary(Bovar_2KmDispDist_output_tab$AUC_cv)
# summary(Bovar_4KmDispDist_output_tab$AUC_cv)
summary(Bovar_6KmDispDist_output_tab$AUC_cv)
summary(Bovar_8KmDispDist_output_tab$AUC_cv)
summary(Bovar_10KmDispDist_output_tab$AUC_cv)

summary(Bovar_DefaultDispDist_output_tab$AUC_train)
summary(Bovar_300mDispDist_output_tab$AUC_train)
summary(Bovar_1KmDispDist_output_tab$AUC_train)
summary(Bovar_2KmDispDist_output_tab$AUC_train)
# summary(Bovar_4KmDispDist_output_tab$AUC_train)
summary(Bovar_6KmDispDist_output_tab$AUC_train)
summary(Bovar_8KmDispDist_output_tab$AUC_train)
summary(Bovar_10KmDispDist_output_tab$AUC_train)

summary(Bovar_DefaultDispDist_output_tab$ntrees)
summary(Bovar_300mDispDist_output_tab$ntrees)
summary(Bovar_1KmDispDist_output_tab$ntrees)
summary(Bovar_2KmDispDist_output_tab$ntrees)
# summary(Bovar_4KmDispDist_output_tab$ntrees)
summary(Bovar_6KmDispDist_output_tab$ntrees)
summary(Bovar_8KmDispDist_output_tab$ntrees)
summary(Bovar_10KmDispDist_output_tab$ntrees)

# ###Get avg_#trees over the AUC>=.75 threshold
# Bovar_mean_acceptable_ntrees <- mean(Bovar_output_tab$ntrees[Bovar_output_tab$AUC_cv>=0.75])
# Bovar_mean_acceptable_ntrees
# noTopo_mean_acceptable_ntrees <- mean(noTopo_output_tab$ntrees[noTopo_output_tab$AUC_cv>=0.75])
# noTopo_mean_acceptable_ntrees
# Total_mean_acceptable_ntrees <- mean(Bovar_mean_acceptable_ntrees,noTopo_mean_acceptable_ntrees)
# Total_mean_acceptable_ntrees



################################################################################
##### Hyarb ####################################################################


#### 2Km #######################################################################

# Create the output table that contains all the values of each of the runs
Hyarb_2KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb


# Create vectors to record performance measures
Hyarb_2KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Hyarb_2KmDispDist_AUC_vec = vector() # training AUC
Hyarb_2KmDispDist_nt_vec = vector() #Number of trees


#Create vectors with all the values of a certain predictor along all the runs
Hyarb_2KmDispDist_HSI_vec = vector() 
Hyarb_2KmDispDist_EgoSize_vec = vector() 
Hyarb_2KmDispDist_strength_vec  = vector()
Hyarb_2KmDispDist_deg_vec  = vector()
Hyarb_2KmDispDist_habAv_vec  = vector()
Hyarb_2KmDispDist_unw_b_c_vec  = vector()
Hyarb_2KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Hyarb_2KmDispDist_predict_mat = matrix(nrow = length(Hyarb_2KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Hyarb_2KmDispDist = gbm.step(data=Hyarb_2KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Hyarb_2KmDispDist)$var, imp = summary(gbm_mod_Hyarb_2KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Hyarb_2KmDispDist_HSI_vec = append(Hyarb_2KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Hyarb_2KmDispDist_EgoSize_vec = append(Hyarb_2KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Hyarb_2KmDispDist_strength_vec  = append(Hyarb_2KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Hyarb_2KmDispDist_deg_vec  = append(Hyarb_2KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Hyarb_2KmDispDist_habAv_vec  = append(Hyarb_2KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Hyarb_2KmDispDist_unw_b_c_vec  = append(Hyarb_2KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Hyarb_2KmDispDist_Patch_Area_vec  = append(Hyarb_2KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Hyarb_2KmDispDist_AUC_cv_vec = append(Hyarb_2KmDispDist_AUC_cv_vec, gbm_mod_Hyarb_2KmDispDist$cv.statistics$discrimination.mean)
  Hyarb_2KmDispDist_AUC_vec = append(Hyarb_2KmDispDist_AUC_vec, gbm_mod_Hyarb_2KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Hyarb_2KmDispDist$n.trees
  Hyarb_2KmDispDist_nt_vec = append(Hyarb_2KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Hyarb_2KmDispDist_predict_mat[,i] = predict(gbm_mod_Hyarb_2KmDispDist, Hyarb_2KmDispDist_stattest, gbm_mod_Hyarb_2KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Hyarb_2KmDispDist_output_tab[,"AUC_train"] = Hyarb_2KmDispDist_AUC_vec
Hyarb_2KmDispDist_output_tab[,"AUC_cv"] = Hyarb_2KmDispDist_AUC_cv_vec
Hyarb_2KmDispDist_output_tab[,"ntrees"] = Hyarb_2KmDispDist_nt_vec

#Var. importance columns
Hyarb_2KmDispDist_output_tab[,"HSI_imp"] = Hyarb_2KmDispDist_HSI_vec
Hyarb_2KmDispDist_output_tab[,"EgoSize_imp"] = Hyarb_2KmDispDist_EgoSize_vec
Hyarb_2KmDispDist_output_tab[,"strength_imp"] = Hyarb_2KmDispDist_strength_vec
Hyarb_2KmDispDist_output_tab[,"deg_imp"] = Hyarb_2KmDispDist_deg_vec
Hyarb_2KmDispDist_output_tab[,"habAv_imp"] = Hyarb_2KmDispDist_habAv_vec
Hyarb_2KmDispDist_output_tab[,"unw_b_c_imp"] = Hyarb_2KmDispDist_unw_b_c_vec
Hyarb_2KmDispDist_output_tab[,"Patch_Area_imp"] = Hyarb_2KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Hyarb_2KmDispDist_HSI_tab = data.frame(imp = Hyarb_2KmDispDist_HSI_vec)
Hyarb_2KmDispDist_EgoSize_tab = data.frame(imp = Hyarb_2KmDispDist_EgoSize_vec)
Hyarb_2KmDispDist_strength_tab = data.frame(imp = Hyarb_2KmDispDist_strength_vec)
Hyarb_2KmDispDist_deg_tab = data.frame(imp = Hyarb_2KmDispDist_deg_vec)
Hyarb_2KmDispDist_habAv_tab = data.frame(imp = Hyarb_2KmDispDist_habAv_vec)
Hyarb_2KmDispDist_unw_b_c_tab = data.frame(imp = Hyarb_2KmDispDist_unw_b_c_vec)
Hyarb_2KmDispDist_Patch_Area_tab = data.frame(imp = Hyarb_2KmDispDist_Patch_Area_vec)

Hyarb_2KmDispDist_HSI_tab$variable = "HSI"
Hyarb_2KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Hyarb_2KmDispDist_strength_tab$variable = "Strength"
Hyarb_2KmDispDist_deg_tab$variable = "Degree"
Hyarb_2KmDispDist_habAv_tab$variable = "Hab. Av."
Hyarb_2KmDispDist_unw_b_c_tab$variable = "B.C."
Hyarb_2KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Hyarb_2KmDispDist_AUC_cv_tab = data.frame(value = Hyarb_2KmDispDist_AUC_cv_vec)
Hyarb_2KmDispDist_AUC_train_tab = data.frame(value = Hyarb_2KmDispDist_AUC_vec)
#Make label of model for plot
Hyarb_2KmDispDist_AUC_cv_tab$model = "Hyarb_2KmDispDist"
Hyarb_2KmDispDist_AUC_train_tab$model = "Hyarb_2KmDispDist"

#combine pred. vars. into new data frame 
Hyarb_2KmDispDist_var_imp_tab = rbind(Hyarb_2KmDispDist_HSI_tab,Hyarb_2KmDispDist_EgoSize_tab,Hyarb_2KmDispDist_strength_tab,Hyarb_2KmDispDist_habAv_tab,Hyarb_2KmDispDist_deg_tab,Hyarb_2KmDispDist_unw_b_c_tab,Hyarb_2KmDispDist_Patch_Area_tab)

ggplot(Hyarb_2KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Hyarb_2KmDispDist_var_imp_tab$imp~Hyarb_2KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "2 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Hyarb_2KmDispDist_HSI_vec)
mean(Hyarb_2KmDispDist_EgoSize_vec)
mean(Hyarb_2KmDispDist_strength_vec)
mean(Hyarb_2KmDispDist_deg_vec)
mean(Hyarb_2KmDispDist_habAv_vec)
mean(Hyarb_2KmDispDist_unw_b_c_vec)
mean(Hyarb_2KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Hyarb_2KmDispDist_output_tab$AUC_cv)
summary(Hyarb_2KmDispDist_output_tab$AUC_train)


#### 300m ######################################################################

# Create the output table that contains all the values of each of the runs
Hyarb_300mDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb

# Create vectors to record performance measures
Hyarb_300mDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Hyarb_300mDispDist_AUC_vec = vector() #

Hyarb_300mDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Hyarb_300mDispDist_HSI_vec = vector() 
Hyarb_300mDispDist_EgoSize_vec = vector() 
Hyarb_300mDispDist_strength_vec  = vector()
Hyarb_300mDispDist_deg_vec  = vector()
Hyarb_300mDispDist_habAv_vec  = vector()
Hyarb_300mDispDist_unw_b_c_vec  = vector()
Hyarb_300mDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Hyarb_300mDispDist_predict_mat = matrix(nrow = length(Hyarb_300mDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Hyarb_300mDispDist = gbm.step(data=Hyarb_300mDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Hyarb_300mDispDist)$var, imp = summary(gbm_mod_Hyarb_300mDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Hyarb_300mDispDist_HSI_vec = append(Hyarb_300mDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Hyarb_300mDispDist_EgoSize_vec = append(Hyarb_300mDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Hyarb_300mDispDist_strength_vec  = append(Hyarb_300mDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Hyarb_300mDispDist_deg_vec  = append(Hyarb_300mDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Hyarb_300mDispDist_habAv_vec  = append(Hyarb_300mDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Hyarb_300mDispDist_unw_b_c_vec  = append(Hyarb_300mDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Hyarb_300mDispDist_Patch_Area_vec  = append(Hyarb_300mDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Hyarb_300mDispDist_AUC_cv_vec = append(Hyarb_300mDispDist_AUC_cv_vec, gbm_mod_Hyarb_300mDispDist$cv.statistics$discrimination.mean)
  Hyarb_300mDispDist_AUC_vec = append(Hyarb_300mDispDist_AUC_vec, gbm_mod_Hyarb_300mDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Hyarb_300mDispDist$n.trees
  Hyarb_300mDispDist_nt_vec = append(Hyarb_300mDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Hyarb_300mDispDist_predict_mat[,i] = predict(gbm_mod_Hyarb_300mDispDist, Hyarb_300mDispDist_stattest, gbm_mod_Hyarb_300mDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Hyarb_300mDispDist_output_tab[,"AUC_train"] = Hyarb_300mDispDist_AUC_vec
Hyarb_300mDispDist_output_tab[,"AUC_cv"] = Hyarb_300mDispDist_AUC_cv_vec
Hyarb_300mDispDist_output_tab[,"ntrees"] = Hyarb_300mDispDist_nt_vec

#Var. importance columns
Hyarb_300mDispDist_output_tab[,"HSI_imp"] = Hyarb_300mDispDist_HSI_vec
Hyarb_300mDispDist_output_tab[,"EgoSize_imp"] = Hyarb_300mDispDist_EgoSize_vec
Hyarb_300mDispDist_output_tab[,"strength_imp"] = Hyarb_300mDispDist_strength_vec
Hyarb_300mDispDist_output_tab[,"deg_imp"] = Hyarb_300mDispDist_deg_vec
Hyarb_300mDispDist_output_tab[,"habAv_imp"] = Hyarb_300mDispDist_habAv_vec
Hyarb_300mDispDist_output_tab[,"unw_b_c_imp"] = Hyarb_300mDispDist_unw_b_c_vec
Hyarb_300mDispDist_output_tab[,"Patch_Area_imp"] = Hyarb_300mDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Hyarb_300mDispDist_HSI_tab = data.frame(imp = Hyarb_300mDispDist_HSI_vec)
Hyarb_300mDispDist_EgoSize_tab = data.frame(imp = Hyarb_300mDispDist_EgoSize_vec)
Hyarb_300mDispDist_strength_tab = data.frame(imp = Hyarb_300mDispDist_strength_vec)
Hyarb_300mDispDist_deg_tab = data.frame(imp = Hyarb_300mDispDist_deg_vec)
Hyarb_300mDispDist_habAv_tab = data.frame(imp = Hyarb_300mDispDist_habAv_vec)
Hyarb_300mDispDist_unw_b_c_tab = data.frame(imp = Hyarb_300mDispDist_unw_b_c_vec)
Hyarb_300mDispDist_Patch_Area_tab = data.frame(imp = Hyarb_300mDispDist_Patch_Area_vec)

Hyarb_300mDispDist_HSI_tab$variable = "HSI"
Hyarb_300mDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Hyarb_300mDispDist_strength_tab$variable = "Strength"
Hyarb_300mDispDist_deg_tab$variable = "Degree"
Hyarb_300mDispDist_habAv_tab$variable = "Hab. Av."
Hyarb_300mDispDist_unw_b_c_tab$variable = "B.C."
Hyarb_300mDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Hyarb_300mDispDist_AUC_cv_tab = data.frame(value = Hyarb_300mDispDist_AUC_cv_vec)
Hyarb_300mDispDist_AUC_train_tab = data.frame(value = Hyarb_300mDispDist_AUC_vec)
#Make label of model for plot
Hyarb_300mDispDist_AUC_cv_tab$model = "Hyarb_300mDispDist"
Hyarb_300mDispDist_AUC_train_tab$model = "Hyarb_300mDispDist"

#combine pred. vars. into new data frame 
Hyarb_300mDispDist_var_imp_tab = rbind(Hyarb_300mDispDist_HSI_tab,Hyarb_300mDispDist_EgoSize_tab,Hyarb_300mDispDist_strength_tab,Hyarb_300mDispDist_habAv_tab,Hyarb_300mDispDist_deg_tab,Hyarb_300mDispDist_unw_b_c_tab,Hyarb_300mDispDist_Patch_Area_tab)

ggplot(Hyarb_300mDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Hyarb_300mDispDist_var_imp_tab$imp~Hyarb_300mDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "300 m maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Hyarb_300mDispDist_HSI_vec)
mean(Hyarb_300mDispDist_EgoSize_vec)
mean(Hyarb_300mDispDist_strength_vec)
mean(Hyarb_300mDispDist_deg_vec)
mean(Hyarb_300mDispDist_habAv_vec)
mean(Hyarb_300mDispDist_unw_b_c_vec)
mean(Hyarb_300mDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Hyarb_300mDispDist_output_tab$AUC_cv)
summary(Hyarb_300mDispDist_output_tab$AUC_train)


#### 10Km ######################################################################

# Create the output table that contains all the values of each of the runs
Hyarb_10KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb

# Create vectors to record performance measures
Hyarb_10KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Hyarb_10KmDispDist_AUC_vec = vector() #

Hyarb_10KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Hyarb_10KmDispDist_HSI_vec = vector() 
Hyarb_10KmDispDist_EgoSize_vec = vector() 
Hyarb_10KmDispDist_strength_vec  = vector()
Hyarb_10KmDispDist_deg_vec  = vector()
Hyarb_10KmDispDist_habAv_vec  = vector()
Hyarb_10KmDispDist_unw_b_c_vec  = vector()
Hyarb_10KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Hyarb_10KmDispDist_predict_mat = matrix(nrow = length(Hyarb_10KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Hyarb_10KmDispDist = gbm.step(data=Hyarb_10KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Hyarb_10KmDispDist)$var, imp = summary(gbm_mod_Hyarb_10KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Hyarb_10KmDispDist_HSI_vec = append(Hyarb_10KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Hyarb_10KmDispDist_EgoSize_vec = append(Hyarb_10KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Hyarb_10KmDispDist_strength_vec  = append(Hyarb_10KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Hyarb_10KmDispDist_deg_vec  = append(Hyarb_10KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Hyarb_10KmDispDist_habAv_vec  = append(Hyarb_10KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Hyarb_10KmDispDist_unw_b_c_vec  = append(Hyarb_10KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Hyarb_10KmDispDist_Patch_Area_vec  = append(Hyarb_10KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Hyarb_10KmDispDist_AUC_cv_vec = append(Hyarb_10KmDispDist_AUC_cv_vec, gbm_mod_Hyarb_10KmDispDist$cv.statistics$discrimination.mean)
  Hyarb_10KmDispDist_AUC_vec = append(Hyarb_10KmDispDist_AUC_vec, gbm_mod_Hyarb_10KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Hyarb_10KmDispDist$n.trees
  Hyarb_10KmDispDist_nt_vec = append(Hyarb_10KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Hyarb_10KmDispDist_predict_mat[,i] = predict(gbm_mod_Hyarb_10KmDispDist, Hyarb_10KmDispDist_stattest, gbm_mod_Hyarb_10KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Hyarb_10KmDispDist_output_tab[,"AUC_train"] = Hyarb_10KmDispDist_AUC_vec
Hyarb_10KmDispDist_output_tab[,"AUC_cv"] = Hyarb_10KmDispDist_AUC_cv_vec
Hyarb_10KmDispDist_output_tab[,"ntrees"] = Hyarb_10KmDispDist_nt_vec

#Var. importance columns
Hyarb_10KmDispDist_output_tab[,"HSI_imp"] = Hyarb_10KmDispDist_HSI_vec
Hyarb_10KmDispDist_output_tab[,"EgoSize_imp"] = Hyarb_10KmDispDist_EgoSize_vec
Hyarb_10KmDispDist_output_tab[,"strength_imp"] = Hyarb_10KmDispDist_strength_vec
Hyarb_10KmDispDist_output_tab[,"deg_imp"] = Hyarb_10KmDispDist_deg_vec
Hyarb_10KmDispDist_output_tab[,"habAv_imp"] = Hyarb_10KmDispDist_habAv_vec
Hyarb_10KmDispDist_output_tab[,"unw_b_c_imp"] = Hyarb_10KmDispDist_unw_b_c_vec
Hyarb_10KmDispDist_output_tab[,"Patch_Area_imp"] = Hyarb_10KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Hyarb_10KmDispDist_HSI_tab = data.frame(imp = Hyarb_10KmDispDist_HSI_vec)
Hyarb_10KmDispDist_EgoSize_tab = data.frame(imp = Hyarb_10KmDispDist_EgoSize_vec)
Hyarb_10KmDispDist_strength_tab = data.frame(imp = Hyarb_10KmDispDist_strength_vec)
Hyarb_10KmDispDist_deg_tab = data.frame(imp = Hyarb_10KmDispDist_deg_vec)
Hyarb_10KmDispDist_habAv_tab = data.frame(imp = Hyarb_10KmDispDist_habAv_vec)
Hyarb_10KmDispDist_unw_b_c_tab = data.frame(imp = Hyarb_10KmDispDist_unw_b_c_vec)
Hyarb_10KmDispDist_Patch_Area_tab = data.frame(imp = Hyarb_10KmDispDist_Patch_Area_vec)

Hyarb_10KmDispDist_HSI_tab$variable = "HSI"
Hyarb_10KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Hyarb_10KmDispDist_strength_tab$variable = "Strength"
Hyarb_10KmDispDist_deg_tab$variable = "Degree"
Hyarb_10KmDispDist_habAv_tab$variable = "Hab. Av."
Hyarb_10KmDispDist_unw_b_c_tab$variable = "B.C."
Hyarb_10KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Hyarb_10KmDispDist_AUC_cv_tab = data.frame(value = Hyarb_10KmDispDist_AUC_cv_vec)
Hyarb_10KmDispDist_AUC_train_tab = data.frame(value = Hyarb_10KmDispDist_AUC_vec)
#Make label of model for plot
Hyarb_10KmDispDist_AUC_cv_tab$model = "Hyarb_10KmDispDist"
Hyarb_10KmDispDist_AUC_train_tab$model = "Hyarb_10KmDispDist"

#combine pred. vars. into new data frame 
Hyarb_10KmDispDist_var_imp_tab = rbind(Hyarb_10KmDispDist_HSI_tab,Hyarb_10KmDispDist_EgoSize_tab,Hyarb_10KmDispDist_strength_tab,Hyarb_10KmDispDist_habAv_tab,Hyarb_10KmDispDist_deg_tab,Hyarb_10KmDispDist_unw_b_c_tab,Hyarb_10KmDispDist_Patch_Area_tab)

ggplot(Hyarb_10KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Hyarb_10KmDispDist_var_imp_tab$imp~Hyarb_10KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "10 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Hyarb_10KmDispDist_HSI_vec)
mean(Hyarb_10KmDispDist_EgoSize_vec)
mean(Hyarb_10KmDispDist_strength_vec)
mean(Hyarb_10KmDispDist_deg_vec)
mean(Hyarb_10KmDispDist_habAv_vec)
mean(Hyarb_10KmDispDist_unw_b_c_vec)
mean(Hyarb_10KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Hyarb_10KmDispDist_output_tab$AUC_cv)
summary(Hyarb_10KmDispDist_output_tab$AUC_train)


####################################### 1Km ####################################

# Create the output table that contains all the values of each of the runs
Hyarb_1KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb

# Create vectors to record performance measures
Hyarb_1KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Hyarb_1KmDispDist_AUC_vec = vector() #

Hyarb_1KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Hyarb_1KmDispDist_HSI_vec = vector() 
Hyarb_1KmDispDist_EgoSize_vec = vector() 
Hyarb_1KmDispDist_strength_vec  = vector()
Hyarb_1KmDispDist_deg_vec  = vector()
Hyarb_1KmDispDist_habAv_vec  = vector()
Hyarb_1KmDispDist_unw_b_c_vec  = vector()
Hyarb_1KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Hyarb_1KmDispDist_predict_mat = matrix(nrow = length(Hyarb_1KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Hyarb_1KmDispDist = gbm.step(data=Hyarb_1KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Hyarb_1KmDispDist)$var, imp = summary(gbm_mod_Hyarb_1KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Hyarb_1KmDispDist_HSI_vec = append(Hyarb_1KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Hyarb_1KmDispDist_EgoSize_vec = append(Hyarb_1KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Hyarb_1KmDispDist_strength_vec  = append(Hyarb_1KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Hyarb_1KmDispDist_deg_vec  = append(Hyarb_1KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Hyarb_1KmDispDist_habAv_vec  = append(Hyarb_1KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Hyarb_1KmDispDist_unw_b_c_vec  = append(Hyarb_1KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Hyarb_1KmDispDist_Patch_Area_vec  = append(Hyarb_1KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Hyarb_1KmDispDist_AUC_cv_vec = append(Hyarb_1KmDispDist_AUC_cv_vec, gbm_mod_Hyarb_1KmDispDist$cv.statistics$discrimination.mean)
  Hyarb_1KmDispDist_AUC_vec = append(Hyarb_1KmDispDist_AUC_vec, gbm_mod_Hyarb_1KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Hyarb_1KmDispDist$n.trees
  Hyarb_1KmDispDist_nt_vec = append(Hyarb_1KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Hyarb_1KmDispDist_predict_mat[,i] = predict(gbm_mod_Hyarb_1KmDispDist, Hyarb_1KmDispDist_stattest, gbm_mod_Hyarb_1KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Hyarb_1KmDispDist_output_tab[,"AUC_train"] = Hyarb_1KmDispDist_AUC_vec
Hyarb_1KmDispDist_output_tab[,"AUC_cv"] = Hyarb_1KmDispDist_AUC_cv_vec
Hyarb_1KmDispDist_output_tab[,"ntrees"] = Hyarb_1KmDispDist_nt_vec

#Var. importance columns
Hyarb_1KmDispDist_output_tab[,"HSI_imp"] = Hyarb_1KmDispDist_HSI_vec
Hyarb_1KmDispDist_output_tab[,"EgoSize_imp"] = Hyarb_1KmDispDist_EgoSize_vec
Hyarb_1KmDispDist_output_tab[,"strength_imp"] = Hyarb_1KmDispDist_strength_vec
Hyarb_1KmDispDist_output_tab[,"deg_imp"] = Hyarb_1KmDispDist_deg_vec
Hyarb_1KmDispDist_output_tab[,"habAv_imp"] = Hyarb_1KmDispDist_habAv_vec
Hyarb_1KmDispDist_output_tab[,"unw_b_c_imp"] = Hyarb_1KmDispDist_unw_b_c_vec
Hyarb_1KmDispDist_output_tab[,"Patch_Area_imp"] = Hyarb_1KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Hyarb_1KmDispDist_HSI_tab = data.frame(imp = Hyarb_1KmDispDist_HSI_vec)
Hyarb_1KmDispDist_EgoSize_tab = data.frame(imp = Hyarb_1KmDispDist_EgoSize_vec)
Hyarb_1KmDispDist_strength_tab = data.frame(imp = Hyarb_1KmDispDist_strength_vec)
Hyarb_1KmDispDist_deg_tab = data.frame(imp = Hyarb_1KmDispDist_deg_vec)
Hyarb_1KmDispDist_habAv_tab = data.frame(imp = Hyarb_1KmDispDist_habAv_vec)
Hyarb_1KmDispDist_unw_b_c_tab = data.frame(imp = Hyarb_1KmDispDist_unw_b_c_vec)
Hyarb_1KmDispDist_Patch_Area_tab = data.frame(imp = Hyarb_1KmDispDist_Patch_Area_vec)

Hyarb_1KmDispDist_HSI_tab$variable = "HSI"
Hyarb_1KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Hyarb_1KmDispDist_strength_tab$variable = "Strength"
Hyarb_1KmDispDist_deg_tab$variable = "Degree"
Hyarb_1KmDispDist_habAv_tab$variable = "Hab. Av."
Hyarb_1KmDispDist_unw_b_c_tab$variable = "B.C."
Hyarb_1KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Hyarb_1KmDispDist_AUC_cv_tab = data.frame(value = Hyarb_1KmDispDist_AUC_cv_vec)
Hyarb_1KmDispDist_AUC_train_tab = data.frame(value = Hyarb_1KmDispDist_AUC_vec)
#Make label of model for plot
Hyarb_1KmDispDist_AUC_cv_tab$model = "Hyarb_1KmDispDist"
Hyarb_1KmDispDist_AUC_train_tab$model = "Hyarb_1KmDispDist"

#combine pred. vars. into new data frame 
Hyarb_1KmDispDist_var_imp_tab = rbind(Hyarb_1KmDispDist_HSI_tab,Hyarb_1KmDispDist_EgoSize_tab,Hyarb_1KmDispDist_strength_tab,Hyarb_1KmDispDist_habAv_tab,Hyarb_1KmDispDist_deg_tab,Hyarb_1KmDispDist_unw_b_c_tab,Hyarb_1KmDispDist_Patch_Area_tab)

ggplot(Hyarb_1KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Hyarb_1KmDispDist_var_imp_tab$imp~Hyarb_1KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "1 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Hyarb_1KmDispDist_HSI_vec)
mean(Hyarb_1KmDispDist_EgoSize_vec)
mean(Hyarb_1KmDispDist_strength_vec)
mean(Hyarb_1KmDispDist_deg_vec)
mean(Hyarb_1KmDispDist_habAv_vec)
mean(Hyarb_1KmDispDist_unw_b_c_vec)
mean(Hyarb_1KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Hyarb_1KmDispDist_output_tab$AUC_cv)
summary(Hyarb_1KmDispDist_output_tab$AUC_train)


####################################### 6Km ####################################

# Create the output table that contains all the values of each of the runs
Hyarb_6KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb

# Create vectors to record performance measures
Hyarb_6KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Hyarb_6KmDispDist_AUC_vec = vector() #

Hyarb_6KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Hyarb_6KmDispDist_HSI_vec = vector() 
Hyarb_6KmDispDist_EgoSize_vec = vector() 
Hyarb_6KmDispDist_strength_vec  = vector()
Hyarb_6KmDispDist_deg_vec  = vector()
Hyarb_6KmDispDist_habAv_vec  = vector()
Hyarb_6KmDispDist_unw_b_c_vec  = vector()
Hyarb_6KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Hyarb_6KmDispDist_predict_mat = matrix(nrow = length(Hyarb_6KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Hyarb_6KmDispDist = gbm.step(data=Hyarb_6KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Hyarb_6KmDispDist)$var, imp = summary(gbm_mod_Hyarb_6KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Hyarb_6KmDispDist_HSI_vec = append(Hyarb_6KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Hyarb_6KmDispDist_EgoSize_vec = append(Hyarb_6KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Hyarb_6KmDispDist_strength_vec  = append(Hyarb_6KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Hyarb_6KmDispDist_deg_vec  = append(Hyarb_6KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Hyarb_6KmDispDist_habAv_vec  = append(Hyarb_6KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Hyarb_6KmDispDist_unw_b_c_vec  = append(Hyarb_6KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Hyarb_6KmDispDist_Patch_Area_vec  = append(Hyarb_6KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Hyarb_6KmDispDist_AUC_cv_vec = append(Hyarb_6KmDispDist_AUC_cv_vec, gbm_mod_Hyarb_6KmDispDist$cv.statistics$discrimination.mean)
  Hyarb_6KmDispDist_AUC_vec = append(Hyarb_6KmDispDist_AUC_vec, gbm_mod_Hyarb_6KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Hyarb_6KmDispDist$n.trees
  Hyarb_6KmDispDist_nt_vec = append(Hyarb_6KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Hyarb_6KmDispDist_predict_mat[,i] = predict(gbm_mod_Hyarb_6KmDispDist, Hyarb_6KmDispDist_stattest, gbm_mod_Hyarb_6KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Hyarb_6KmDispDist_output_tab[,"AUC_train"] = Hyarb_6KmDispDist_AUC_vec
Hyarb_6KmDispDist_output_tab[,"AUC_cv"] = Hyarb_6KmDispDist_AUC_cv_vec
Hyarb_6KmDispDist_output_tab[,"ntrees"] = Hyarb_6KmDispDist_nt_vec

#Var. importance columns
Hyarb_6KmDispDist_output_tab[,"HSI_imp"] = Hyarb_6KmDispDist_HSI_vec
Hyarb_6KmDispDist_output_tab[,"EgoSize_imp"] = Hyarb_6KmDispDist_EgoSize_vec
Hyarb_6KmDispDist_output_tab[,"strength_imp"] = Hyarb_6KmDispDist_strength_vec
Hyarb_6KmDispDist_output_tab[,"deg_imp"] = Hyarb_6KmDispDist_deg_vec
Hyarb_6KmDispDist_output_tab[,"habAv_imp"] = Hyarb_6KmDispDist_habAv_vec
Hyarb_6KmDispDist_output_tab[,"unw_b_c_imp"] = Hyarb_6KmDispDist_unw_b_c_vec
Hyarb_6KmDispDist_output_tab[,"Patch_Area_imp"] = Hyarb_6KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Hyarb_6KmDispDist_HSI_tab = data.frame(imp = Hyarb_6KmDispDist_HSI_vec)
Hyarb_6KmDispDist_EgoSize_tab = data.frame(imp = Hyarb_6KmDispDist_EgoSize_vec)
Hyarb_6KmDispDist_strength_tab = data.frame(imp = Hyarb_6KmDispDist_strength_vec)
Hyarb_6KmDispDist_deg_tab = data.frame(imp = Hyarb_6KmDispDist_deg_vec)
Hyarb_6KmDispDist_habAv_tab = data.frame(imp = Hyarb_6KmDispDist_habAv_vec)
Hyarb_6KmDispDist_unw_b_c_tab = data.frame(imp = Hyarb_6KmDispDist_unw_b_c_vec)
Hyarb_6KmDispDist_Patch_Area_tab = data.frame(imp = Hyarb_6KmDispDist_Patch_Area_vec)

Hyarb_6KmDispDist_HSI_tab$variable = "HSI"
Hyarb_6KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Hyarb_6KmDispDist_strength_tab$variable = "Strength"
Hyarb_6KmDispDist_deg_tab$variable = "Degree"
Hyarb_6KmDispDist_habAv_tab$variable = "Hab. Av."
Hyarb_6KmDispDist_unw_b_c_tab$variable = "B.C."
Hyarb_6KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Hyarb_6KmDispDist_AUC_cv_tab = data.frame(value = Hyarb_6KmDispDist_AUC_cv_vec)
Hyarb_6KmDispDist_AUC_train_tab = data.frame(value = Hyarb_6KmDispDist_AUC_vec)
#Make label of model for plot
Hyarb_6KmDispDist_AUC_cv_tab$model = "Hyarb_6KmDispDist"
Hyarb_6KmDispDist_AUC_train_tab$model = "Hyarb_6KmDispDist"

#combine pred. vars. into new data frame 
Hyarb_6KmDispDist_var_imp_tab = rbind(Hyarb_6KmDispDist_HSI_tab,Hyarb_6KmDispDist_EgoSize_tab,Hyarb_6KmDispDist_strength_tab,Hyarb_6KmDispDist_habAv_tab,Hyarb_6KmDispDist_deg_tab,Hyarb_6KmDispDist_unw_b_c_tab,Hyarb_6KmDispDist_Patch_Area_tab)

ggplot(Hyarb_6KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Hyarb_6KmDispDist_var_imp_tab$imp~Hyarb_6KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "6 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Hyarb_6KmDispDist_HSI_vec)
mean(Hyarb_6KmDispDist_EgoSize_vec)
mean(Hyarb_6KmDispDist_strength_vec)
mean(Hyarb_6KmDispDist_deg_vec)
mean(Hyarb_6KmDispDist_habAv_vec)
mean(Hyarb_6KmDispDist_unw_b_c_vec)
mean(Hyarb_6KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Hyarb_6KmDispDist_output_tab$AUC_cv)
summary(Hyarb_6KmDispDist_output_tab$AUC_train)


####################################### 4Km ####################################

# Create the output table that contains all the values of each of the runs
Hyarb_4KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb

# Create vectors to record performance measures
Hyarb_4KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Hyarb_4KmDispDist_AUC_vec = vector() #

Hyarb_4KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Hyarb_4KmDispDist_HSI_vec = vector()
Hyarb_4KmDispDist_EgoSize_vec = vector()
Hyarb_4KmDispDist_strength_vec  = vector()
Hyarb_4KmDispDist_deg_vec  = vector()
Hyarb_4KmDispDist_habAv_vec  = vector()
Hyarb_4KmDispDist_unw_b_c_vec  = vector()
Hyarb_4KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Hyarb_4KmDispDist_predict_mat = matrix(nrow = length(Hyarb_4KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Hyarb_4KmDispDist = gbm.step(data=Hyarb_4KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Hyarb_4KmDispDist)$var, imp = summary(gbm_mod_Hyarb_4KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Hyarb_4KmDispDist_HSI_vec = append(Hyarb_4KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Hyarb_4KmDispDist_EgoSize_vec = append(Hyarb_4KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Hyarb_4KmDispDist_strength_vec  = append(Hyarb_4KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Hyarb_4KmDispDist_deg_vec  = append(Hyarb_4KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Hyarb_4KmDispDist_habAv_vec  = append(Hyarb_4KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Hyarb_4KmDispDist_unw_b_c_vec  = append(Hyarb_4KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Hyarb_4KmDispDist_Patch_Area_vec  = append(Hyarb_4KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Hyarb_4KmDispDist_AUC_cv_vec = append(Hyarb_4KmDispDist_AUC_cv_vec, gbm_mod_Hyarb_4KmDispDist$cv.statistics$discrimination.mean)
  Hyarb_4KmDispDist_AUC_vec = append(Hyarb_4KmDispDist_AUC_vec, gbm_mod_Hyarb_4KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Hyarb_4KmDispDist$n.trees
  Hyarb_4KmDispDist_nt_vec = append(Hyarb_4KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Hyarb_4KmDispDist_predict_mat[,i] = predict(gbm_mod_Hyarb_4KmDispDist, Hyarb_4KmDispDist_stattest, gbm_mod_Hyarb_4KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Hyarb_4KmDispDist_output_tab[,"AUC_train"] = Hyarb_4KmDispDist_AUC_vec
Hyarb_4KmDispDist_output_tab[,"AUC_cv"] = Hyarb_4KmDispDist_AUC_cv_vec
Hyarb_4KmDispDist_output_tab[,"ntrees"] = Hyarb_4KmDispDist_nt_vec

#Var. importance columns
Hyarb_4KmDispDist_output_tab[,"HSI_imp"] = Hyarb_4KmDispDist_HSI_vec
Hyarb_4KmDispDist_output_tab[,"EgoSize_imp"] = Hyarb_4KmDispDist_EgoSize_vec
Hyarb_4KmDispDist_output_tab[,"strength_imp"] = Hyarb_4KmDispDist_strength_vec
Hyarb_4KmDispDist_output_tab[,"deg_imp"] = Hyarb_4KmDispDist_deg_vec
Hyarb_4KmDispDist_output_tab[,"habAv_imp"] = Hyarb_4KmDispDist_habAv_vec
Hyarb_4KmDispDist_output_tab[,"unw_b_c_imp"] = Hyarb_4KmDispDist_unw_b_c_vec
Hyarb_4KmDispDist_output_tab[,"Patch_Area_imp"] = Hyarb_4KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Hyarb_4KmDispDist_HSI_tab = data.frame(imp = Hyarb_4KmDispDist_HSI_vec)
Hyarb_4KmDispDist_EgoSize_tab = data.frame(imp = Hyarb_4KmDispDist_EgoSize_vec)
Hyarb_4KmDispDist_strength_tab = data.frame(imp = Hyarb_4KmDispDist_strength_vec)
Hyarb_4KmDispDist_deg_tab = data.frame(imp = Hyarb_4KmDispDist_deg_vec)
Hyarb_4KmDispDist_habAv_tab = data.frame(imp = Hyarb_4KmDispDist_habAv_vec)
Hyarb_4KmDispDist_unw_b_c_tab = data.frame(imp = Hyarb_4KmDispDist_unw_b_c_vec)
Hyarb_4KmDispDist_Patch_Area_tab = data.frame(imp = Hyarb_4KmDispDist_Patch_Area_vec)

Hyarb_4KmDispDist_HSI_tab$variable = "HSI"
Hyarb_4KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Hyarb_4KmDispDist_strength_tab$variable = "Strength"
Hyarb_4KmDispDist_deg_tab$variable = "Degree"
Hyarb_4KmDispDist_habAv_tab$variable = "Hab. Av."
Hyarb_4KmDispDist_unw_b_c_tab$variable = "B.C."
Hyarb_4KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Hyarb_4KmDispDist_AUC_cv_tab = data.frame(value = Hyarb_4KmDispDist_AUC_cv_vec)
Hyarb_4KmDispDist_AUC_train_tab = data.frame(value = Hyarb_4KmDispDist_AUC_vec)
#Make label of model for plot
Hyarb_4KmDispDist_AUC_cv_tab$model = "Hyarb_4KmDispDist"
Hyarb_4KmDispDist_AUC_train_tab$model = "Hyarb_4KmDispDist"

#combine pred. vars. into new data frame
Hyarb_4KmDispDist_var_imp_tab = rbind(Hyarb_4KmDispDist_HSI_tab,Hyarb_4KmDispDist_EgoSize_tab,Hyarb_4KmDispDist_strength_tab,Hyarb_4KmDispDist_habAv_tab,Hyarb_4KmDispDist_deg_tab,Hyarb_4KmDispDist_unw_b_c_tab,Hyarb_4KmDispDist_Patch_Area_tab)

ggplot(Hyarb_4KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Hyarb_4KmDispDist_var_imp_tab$imp~Hyarb_4KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "4 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Hyarb_4KmDispDist_HSI_vec)
mean(Hyarb_4KmDispDist_EgoSize_vec)
mean(Hyarb_4KmDispDist_strength_vec)
mean(Hyarb_4KmDispDist_deg_vec)
mean(Hyarb_4KmDispDist_habAv_vec)
mean(Hyarb_4KmDispDist_unw_b_c_vec)
mean(Hyarb_4KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Hyarb_4KmDispDist_output_tab$AUC_cv)
summary(Hyarb_4KmDispDist_output_tab$AUC_train)


####################################### 8Km ####################################

# Create the output table that contains all the values of each of the runs
Hyarb_8KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb

# Create vectors to record performance measures
Hyarb_8KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Hyarb_8KmDispDist_AUC_vec = vector() #

Hyarb_8KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Hyarb_8KmDispDist_HSI_vec = vector() 
Hyarb_8KmDispDist_EgoSize_vec = vector() 
Hyarb_8KmDispDist_strength_vec  = vector()
Hyarb_8KmDispDist_deg_vec  = vector()
Hyarb_8KmDispDist_habAv_vec  = vector()
Hyarb_8KmDispDist_unw_b_c_vec  = vector()
Hyarb_8KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Hyarb_8KmDispDist_predict_mat = matrix(nrow = length(Hyarb_8KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Hyarb_8KmDispDist = gbm.step(data=Hyarb_8KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Hyarb_8KmDispDist)$var, imp = summary(gbm_mod_Hyarb_8KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Hyarb_8KmDispDist_HSI_vec = append(Hyarb_8KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Hyarb_8KmDispDist_EgoSize_vec = append(Hyarb_8KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Hyarb_8KmDispDist_strength_vec  = append(Hyarb_8KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Hyarb_8KmDispDist_deg_vec  = append(Hyarb_8KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Hyarb_8KmDispDist_habAv_vec  = append(Hyarb_8KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Hyarb_8KmDispDist_unw_b_c_vec  = append(Hyarb_8KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Hyarb_8KmDispDist_Patch_Area_vec  = append(Hyarb_8KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Hyarb_8KmDispDist_AUC_cv_vec = append(Hyarb_8KmDispDist_AUC_cv_vec, gbm_mod_Hyarb_8KmDispDist$cv.statistics$discrimination.mean)
  Hyarb_8KmDispDist_AUC_vec = append(Hyarb_8KmDispDist_AUC_vec, gbm_mod_Hyarb_8KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Hyarb_8KmDispDist$n.trees
  Hyarb_8KmDispDist_nt_vec = append(Hyarb_8KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Hyarb_8KmDispDist_predict_mat[,i] = predict(gbm_mod_Hyarb_8KmDispDist, Hyarb_8KmDispDist_stattest, gbm_mod_Hyarb_8KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Hyarb_8KmDispDist_output_tab[,"AUC_train"] = Hyarb_8KmDispDist_AUC_vec
Hyarb_8KmDispDist_output_tab[,"AUC_cv"] = Hyarb_8KmDispDist_AUC_cv_vec
Hyarb_8KmDispDist_output_tab[,"ntrees"] = Hyarb_8KmDispDist_nt_vec

#Var. importance columns
Hyarb_8KmDispDist_output_tab[,"HSI_imp"] = Hyarb_8KmDispDist_HSI_vec
Hyarb_8KmDispDist_output_tab[,"EgoSize_imp"] = Hyarb_8KmDispDist_EgoSize_vec
Hyarb_8KmDispDist_output_tab[,"strength_imp"] = Hyarb_8KmDispDist_strength_vec
Hyarb_8KmDispDist_output_tab[,"deg_imp"] = Hyarb_8KmDispDist_deg_vec
Hyarb_8KmDispDist_output_tab[,"habAv_imp"] = Hyarb_8KmDispDist_habAv_vec
Hyarb_8KmDispDist_output_tab[,"unw_b_c_imp"] = Hyarb_8KmDispDist_unw_b_c_vec
Hyarb_8KmDispDist_output_tab[,"Patch_Area_imp"] = Hyarb_8KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Hyarb_8KmDispDist_HSI_tab = data.frame(imp = Hyarb_8KmDispDist_HSI_vec)
Hyarb_8KmDispDist_EgoSize_tab = data.frame(imp = Hyarb_8KmDispDist_EgoSize_vec)
Hyarb_8KmDispDist_strength_tab = data.frame(imp = Hyarb_8KmDispDist_strength_vec)
Hyarb_8KmDispDist_deg_tab = data.frame(imp = Hyarb_8KmDispDist_deg_vec)
Hyarb_8KmDispDist_habAv_tab = data.frame(imp = Hyarb_8KmDispDist_habAv_vec)
Hyarb_8KmDispDist_unw_b_c_tab = data.frame(imp = Hyarb_8KmDispDist_unw_b_c_vec)
Hyarb_8KmDispDist_Patch_Area_tab = data.frame(imp = Hyarb_8KmDispDist_Patch_Area_vec)

Hyarb_8KmDispDist_HSI_tab$variable = "HSI"
Hyarb_8KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Hyarb_8KmDispDist_strength_tab$variable = "Strength"
Hyarb_8KmDispDist_deg_tab$variable = "Degree"
Hyarb_8KmDispDist_habAv_tab$variable = "Hab. Av."
Hyarb_8KmDispDist_unw_b_c_tab$variable = "B.C."
Hyarb_8KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Hyarb_8KmDispDist_AUC_cv_tab = data.frame(value = Hyarb_8KmDispDist_AUC_cv_vec)
Hyarb_8KmDispDist_AUC_train_tab = data.frame(value = Hyarb_8KmDispDist_AUC_vec)
#Make label of model for plot
Hyarb_8KmDispDist_AUC_cv_tab$model = "Hyarb_8KmDispDist"
Hyarb_8KmDispDist_AUC_train_tab$model = "Hyarb_8KmDispDist"

#combine pred. vars. into new data frame 
Hyarb_8KmDispDist_var_imp_tab = rbind(Hyarb_8KmDispDist_HSI_tab,Hyarb_8KmDispDist_EgoSize_tab,Hyarb_8KmDispDist_strength_tab,Hyarb_8KmDispDist_habAv_tab,Hyarb_8KmDispDist_deg_tab,Hyarb_8KmDispDist_unw_b_c_tab,Hyarb_8KmDispDist_Patch_Area_tab)

ggplot(Hyarb_8KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Hyarb_8KmDispDist_var_imp_tab$imp~Hyarb_8KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "8 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Hyarb_8KmDispDist_HSI_vec)
mean(Hyarb_8KmDispDist_EgoSize_vec)
mean(Hyarb_8KmDispDist_strength_vec)
mean(Hyarb_8KmDispDist_deg_vec)
mean(Hyarb_8KmDispDist_habAv_vec)
mean(Hyarb_8KmDispDist_unw_b_c_vec)
mean(Hyarb_8KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Hyarb_8KmDispDist_output_tab$AUC_cv)
summary(Hyarb_8KmDispDist_output_tab$AUC_train)


####################################################################################################
#### Compare scores between networks w/d0 variations of the same species ###########################
### Hyarb ##################

#Import evaluation df's of original run with species-specific dispersal distance
Hyarb_DefaultDispDist_AUC_cv_tab <- read.csv("Hyarb_DefaultDispDist_AUC_cv_tab.csv")
Hyarb_DefaultDispDist_AUC_train_tab <- read.csv("Hyarb_DefaultDispDist_AUC_train_tab.csv")
Hyarb_DefaultDispDist_var_imp_tab <- read.csv("Hyarb_DefaultDispDist_var_imp_tab.csv")
Hyarb_DefaultDispDist_output_tab <- read.csv("Hyarb_DefaultDispDist_BRToutput_tab.csv")
Hyarb_DefaultDispDist_predict_df <- read.csv("Hyarb_DefaultDispDist_predict_df.csv")

head(Hyarb_DefaultDispDist_AUC_cv_tab)
Hyarb_DefaultDispDist_AUC_cv_tab$X <- NULL
Hyarb_DefaultDispDist_AUC_train_tab$X <- NULL

#### Make a dataframe for plotting overlaying histograms in R
### cv AUC
Hyarb_AUC_cv_tab = rbind(Hyarb_DefaultDispDist_AUC_cv_tab, Hyarb_300mDispDist_AUC_cv_tab, Hyarb_1KmDispDist_AUC_cv_tab, Hyarb_2KmDispDist_AUC_cv_tab, 
                         Hyarb_4KmDispDist_AUC_cv_tab,
                         Hyarb_6KmDispDist_AUC_cv_tab, Hyarb_8KmDispDist_AUC_cv_tab, Hyarb_10KmDispDist_AUC_cv_tab)
#Change order of factors 
Hyarb_AUC_cv_tab$model<- as.factor(Hyarb_AUC_cv_tab$model)
levels(Hyarb_AUC_cv_tab$model)
Hyarb_AUC_cv_tab$model<-factor(Hyarb_AUC_cv_tab$model, levels=c("Hyarb_300mDispDist", "Hyarb_1KmDispDist", "Hyarb_2KmDispDist", 
                                                                "Hyarb", "Hyarb_4KmDispDist", "Hyarb_6KmDispDist", 
                                                                "Hyarb_8KmDispDist", "Hyarb_10KmDispDist"))

#Plot
ggplot(Hyarb_AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Hyarb_AUC_cv_tab$value~Hyarb_AUC_cv_tab$model, ylab= "Cross-validated AUC")


### train AUC
Hyarb_AUC_train_tab = rbind(Hyarb_DefaultDispDist_AUC_train_tab, Hyarb_300mDispDist_AUC_train_tab, Hyarb_1KmDispDist_AUC_train_tab, Hyarb_2KmDispDist_AUC_train_tab, 
                            Hyarb_4KmDispDist_AUC_train_tab,
                            Hyarb_6KmDispDist_AUC_train_tab, Hyarb_8KmDispDist_AUC_train_tab, Hyarb_10KmDispDist_AUC_train_tab)
#Change order of factors 
levels(Hyarb_AUC_train_tab$model)
Hyarb_AUC_train_tab$model<-factor(Hyarb_AUC_train_tab$model, levels=c("Hyarb_300mDispDist", "Hyarb_1KmDispDist", "Hyarb_2KmDispDist", 
                                                                      "Hyarb", "Hyarb_4KmDispDist", "Hyarb_6KmDispDist", 
                                                                      "Hyarb_8KmDispDist", "Hyarb_10KmDispDist"))
#Plot
ggplot(Hyarb_AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Hyarb_AUC_train_tab$value~Hyarb_AUC_train_tab$model, ylab= "Training AUC")

summary(Hyarb_DefaultDispDist_output_tab$AUC_cv)
summary(Hyarb_300mDispDist_output_tab$AUC_cv)
summary(Hyarb_1KmDispDist_output_tab$AUC_cv)
summary(Hyarb_2KmDispDist_output_tab$AUC_cv)
summary(Hyarb_4KmDispDist_output_tab$AUC_cv)
summary(Hyarb_6KmDispDist_output_tab$AUC_cv)
summary(Hyarb_8KmDispDist_output_tab$AUC_cv)
summary(Hyarb_10KmDispDist_output_tab$AUC_cv)

summary(Hyarb_DefaultDispDist_output_tab$AUC_train)
summary(Hyarb_300mDispDist_output_tab$AUC_train)
summary(Hyarb_1KmDispDist_output_tab$AUC_train)
summary(Hyarb_2KmDispDist_output_tab$AUC_train)
summary(Hyarb_4KmDispDist_output_tab$AUC_train)
summary(Hyarb_6KmDispDist_output_tab$AUC_train)
summary(Hyarb_8KmDispDist_output_tab$AUC_train)
summary(Hyarb_10KmDispDist_output_tab$AUC_train)

summary(Hyarb_DefaultDispDist_output_tab$ntrees)
summary(Hyarb_300mDispDist_output_tab$ntrees)
summary(Hyarb_1KmDispDist_output_tab$ntrees)
summary(Hyarb_2KmDispDist_output_tab$ntrees)
summary(Hyarb_4KmDispDist_output_tab$ntrees)
summary(Hyarb_6KmDispDist_output_tab$ntrees)
summary(Hyarb_8KmDispDist_output_tab$ntrees)
summary(Hyarb_10KmDispDist_output_tab$ntrees)



################################################################################
##### Alobs ####################################################################


#### 2Km #######################################################################

# Create the output table that contains all the values of each of the runs
Alobs_2KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Alobs

# Create vectors to record performance measures
Alobs_2KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Alobs_2KmDispDist_AUC_vec = vector() #

Alobs_2KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Alobs_2KmDispDist_HSI_vec = vector()
Alobs_2KmDispDist_EgoSize_vec = vector()
Alobs_2KmDispDist_strength_vec  = vector()
Alobs_2KmDispDist_deg_vec  = vector()
Alobs_2KmDispDist_habAv_vec  = vector()
Alobs_2KmDispDist_unw_b_c_vec  = vector()
Alobs_2KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Alobs_2KmDispDist_predict_mat = matrix(nrow = length(Alobs_2KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Alobs_2KmDispDist = gbm.step(data=Alobs_2KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Alobs_2KmDispDist)$var, imp = summary(gbm_mod_Alobs_2KmDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Alobs_2KmDispDist_HSI_vec = append(Alobs_2KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Alobs_2KmDispDist_EgoSize_vec = append(Alobs_2KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Alobs_2KmDispDist_strength_vec  = append(Alobs_2KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Alobs_2KmDispDist_deg_vec  = append(Alobs_2KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Alobs_2KmDispDist_habAv_vec  = append(Alobs_2KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Alobs_2KmDispDist_unw_b_c_vec  = append(Alobs_2KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Alobs_2KmDispDist_Patch_Area_vec  = append(Alobs_2KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Alobs_2KmDispDist_AUC_cv_vec = append(Alobs_2KmDispDist_AUC_cv_vec, gbm_mod_Alobs_2KmDispDist$cv.statistics$discrimination.mean)
  Alobs_2KmDispDist_AUC_vec = append(Alobs_2KmDispDist_AUC_vec, gbm_mod_Alobs_2KmDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Alobs_2KmDispDist$n.trees
  Alobs_2KmDispDist_nt_vec = append(Alobs_2KmDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Alobs_2KmDispDist_predict_mat[,i] = predict(gbm_mod_Alobs_2KmDispDist, Alobs_2KmDispDist_stattest, gbm_mod_Alobs_2KmDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Alobs_2KmDispDist_output_tab[,"AUC_train"] = Alobs_2KmDispDist_AUC_vec
Alobs_2KmDispDist_output_tab[,"AUC_cv"] = Alobs_2KmDispDist_AUC_cv_vec
Alobs_2KmDispDist_output_tab[,"ntrees"] = Alobs_2KmDispDist_nt_vec

#Var. importance columns
Alobs_2KmDispDist_output_tab[,"HSI_imp"] = Alobs_2KmDispDist_HSI_vec
Alobs_2KmDispDist_output_tab[,"EgoSize_imp"] = Alobs_2KmDispDist_EgoSize_vec
Alobs_2KmDispDist_output_tab[,"strength_imp"] = Alobs_2KmDispDist_strength_vec
Alobs_2KmDispDist_output_tab[,"deg_imp"] = Alobs_2KmDispDist_deg_vec
Alobs_2KmDispDist_output_tab[,"habAv_imp"] = Alobs_2KmDispDist_habAv_vec
Alobs_2KmDispDist_output_tab[,"unw_b_c_imp"] = Alobs_2KmDispDist_unw_b_c_vec
Alobs_2KmDispDist_output_tab[,"Patch_Area_imp"] = Alobs_2KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Alobs_2KmDispDist_HSI_tab = data.frame(imp = Alobs_2KmDispDist_HSI_vec)
Alobs_2KmDispDist_EgoSize_tab = data.frame(imp = Alobs_2KmDispDist_EgoSize_vec)
Alobs_2KmDispDist_strength_tab = data.frame(imp = Alobs_2KmDispDist_strength_vec)
Alobs_2KmDispDist_deg_tab = data.frame(imp = Alobs_2KmDispDist_deg_vec)
Alobs_2KmDispDist_habAv_tab = data.frame(imp = Alobs_2KmDispDist_habAv_vec)
Alobs_2KmDispDist_unw_b_c_tab = data.frame(imp = Alobs_2KmDispDist_unw_b_c_vec)
Alobs_2KmDispDist_Patch_Area_tab = data.frame(imp = Alobs_2KmDispDist_Patch_Area_vec)

Alobs_2KmDispDist_HSI_tab$variable = "HSI"
Alobs_2KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Alobs_2KmDispDist_strength_tab$variable = "Strength"
Alobs_2KmDispDist_deg_tab$variable = "Degree"
Alobs_2KmDispDist_habAv_tab$variable = "Hab. Av."
Alobs_2KmDispDist_unw_b_c_tab$variable = "B.C."
Alobs_2KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Alobs_2KmDispDist_AUC_cv_tab = data.frame(value = Alobs_2KmDispDist_AUC_cv_vec)
Alobs_2KmDispDist_AUC_train_tab = data.frame(value = Alobs_2KmDispDist_AUC_vec)
#Make label of model for plot
Alobs_2KmDispDist_AUC_cv_tab$model = "Alobs_2KmDispDist"
Alobs_2KmDispDist_AUC_train_tab$model = "Alobs_2KmDispDist"

#combine pred. vars. into new data frame
Alobs_2KmDispDist_var_imp_tab = rbind(Alobs_2KmDispDist_HSI_tab,Alobs_2KmDispDist_EgoSize_tab,Alobs_2KmDispDist_strength_tab,Alobs_2KmDispDist_habAv_tab,Alobs_2KmDispDist_deg_tab,Alobs_2KmDispDist_unw_b_c_tab,Alobs_2KmDispDist_Patch_Area_tab)

ggplot(Alobs_2KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)

boxplot(Alobs_2KmDispDist_var_imp_tab$imp~Alobs_2KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "2 km maximum dispersal distance", 
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Alobs_2KmDispDist_HSI_vec)
mean(Alobs_2KmDispDist_EgoSize_vec)
mean(Alobs_2KmDispDist_strength_vec)
mean(Alobs_2KmDispDist_deg_vec)
mean(Alobs_2KmDispDist_habAv_vec)
mean(Alobs_2KmDispDist_unw_b_c_vec)
mean(Alobs_2KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Alobs_2KmDispDist_output_tab$AUC_cv)
summary(Alobs_2KmDispDist_output_tab$AUC_train)


# #### 300m ####################################################################

# Create the output table that contains all the values of each of the runs
Alobs_300mDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Alobs

# Create vectors to record performance measures
Alobs_300mDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Alobs_300mDispDist_AUC_vec = vector() #

Alobs_300mDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Alobs_300mDispDist_HSI_vec = vector()
Alobs_300mDispDist_EgoSize_vec = vector()
Alobs_300mDispDist_strength_vec  = vector()
Alobs_300mDispDist_deg_vec  = vector()
Alobs_300mDispDist_habAv_vec  = vector()
Alobs_300mDispDist_unw_b_c_vec  = vector()
Alobs_300mDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Alobs_300mDispDist_predict_mat = matrix(nrow = length(Alobs_300mDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Alobs_300mDispDist = gbm.step(data=Alobs_300mDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Alobs_300mDispDist)$var, imp = summary(gbm_mod_Alobs_300mDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Alobs_300mDispDist_HSI_vec = append(Alobs_300mDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Alobs_300mDispDist_EgoSize_vec = append(Alobs_300mDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Alobs_300mDispDist_strength_vec  = append(Alobs_300mDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Alobs_300mDispDist_deg_vec  = append(Alobs_300mDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Alobs_300mDispDist_habAv_vec  = append(Alobs_300mDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Alobs_300mDispDist_unw_b_c_vec  = append(Alobs_300mDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Alobs_300mDispDist_Patch_Area_vec  = append(Alobs_300mDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Alobs_300mDispDist_AUC_cv_vec = append(Alobs_300mDispDist_AUC_cv_vec, gbm_mod_Alobs_300mDispDist$cv.statistics$discrimination.mean)
  Alobs_300mDispDist_AUC_vec = append(Alobs_300mDispDist_AUC_vec, gbm_mod_Alobs_300mDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Alobs_300mDispDist$n.trees
  Alobs_300mDispDist_nt_vec = append(Alobs_300mDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Alobs_300mDispDist_predict_mat[,i] = predict(gbm_mod_Alobs_300mDispDist, Alobs_300mDispDist_stattest, gbm_mod_Alobs_300mDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Alobs_300mDispDist_output_tab[,"AUC_train"] = Alobs_300mDispDist_AUC_vec
Alobs_300mDispDist_output_tab[,"AUC_cv"] = Alobs_300mDispDist_AUC_cv_vec
Alobs_300mDispDist_output_tab[,"ntrees"] = Alobs_300mDispDist_nt_vec

#Var. importance columns
Alobs_300mDispDist_output_tab[,"HSI_imp"] = Alobs_300mDispDist_HSI_vec
Alobs_300mDispDist_output_tab[,"EgoSize_imp"] = Alobs_300mDispDist_EgoSize_vec
Alobs_300mDispDist_output_tab[,"strength_imp"] = Alobs_300mDispDist_strength_vec
Alobs_300mDispDist_output_tab[,"deg_imp"] = Alobs_300mDispDist_deg_vec
Alobs_300mDispDist_output_tab[,"habAv_imp"] = Alobs_300mDispDist_habAv_vec
Alobs_300mDispDist_output_tab[,"unw_b_c_imp"] = Alobs_300mDispDist_unw_b_c_vec
Alobs_300mDispDist_output_tab[,"Patch_Area_imp"] = Alobs_300mDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Alobs_300mDispDist_HSI_tab = data.frame(imp = Alobs_300mDispDist_HSI_vec)
Alobs_300mDispDist_EgoSize_tab = data.frame(imp = Alobs_300mDispDist_EgoSize_vec)
Alobs_300mDispDist_strength_tab = data.frame(imp = Alobs_300mDispDist_strength_vec)
Alobs_300mDispDist_deg_tab = data.frame(imp = Alobs_300mDispDist_deg_vec)
Alobs_300mDispDist_habAv_tab = data.frame(imp = Alobs_300mDispDist_habAv_vec)
Alobs_300mDispDist_unw_b_c_tab = data.frame(imp = Alobs_300mDispDist_unw_b_c_vec)
Alobs_300mDispDist_Patch_Area_tab = data.frame(imp = Alobs_300mDispDist_Patch_Area_vec)

Alobs_300mDispDist_HSI_tab$variable = "HSI"
Alobs_300mDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Alobs_300mDispDist_strength_tab$variable = "Strength"
Alobs_300mDispDist_deg_tab$variable = "Degree"
Alobs_300mDispDist_habAv_tab$variable = "Hab. Av."
Alobs_300mDispDist_unw_b_c_tab$variable = "B.C."
Alobs_300mDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Alobs_300mDispDist_AUC_cv_tab = data.frame(value = Alobs_300mDispDist_AUC_cv_vec)
Alobs_300mDispDist_AUC_train_tab = data.frame(value = Alobs_300mDispDist_AUC_vec)
#Make label of model for plot
Alobs_300mDispDist_AUC_cv_tab$model = "Alobs_300mDispDist"
Alobs_300mDispDist_AUC_train_tab$model = "Alobs_300mDispDist"

#combine pred. vars. into new data frame
Alobs_300mDispDist_var_imp_tab = rbind(Alobs_300mDispDist_HSI_tab,Alobs_300mDispDist_EgoSize_tab,Alobs_300mDispDist_strength_tab,Alobs_300mDispDist_habAv_tab,Alobs_300mDispDist_deg_tab,Alobs_300mDispDist_unw_b_c_tab,Alobs_300mDispDist_Patch_Area_tab)

ggplot(Alobs_300mDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Alobs_300mDispDist_var_imp_tab$imp~Alobs_300mDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "300 m maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Alobs_300mDispDist_HSI_vec)
mean(Alobs_300mDispDist_EgoSize_vec)
mean(Alobs_300mDispDist_strength_vec)
mean(Alobs_300mDispDist_deg_vec)
mean(Alobs_300mDispDist_habAv_vec)
mean(Alobs_300mDispDist_unw_b_c_vec)
mean(Alobs_300mDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Alobs_300mDispDist_output_tab$AUC_cv)
summary(Alobs_300mDispDist_output_tab$AUC_train)


# #### 10Km ####################################################################

# Create the output table that contains all the values of each of the runs
Alobs_10KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Alobs

# Create vectors to record performance measures
Alobs_10KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Alobs_10KmDispDist_AUC_vec = vector() #

Alobs_10KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Alobs_10KmDispDist_HSI_vec = vector()
Alobs_10KmDispDist_EgoSize_vec = vector()
Alobs_10KmDispDist_strength_vec  = vector()
Alobs_10KmDispDist_deg_vec  = vector()
Alobs_10KmDispDist_habAv_vec  = vector()
Alobs_10KmDispDist_unw_b_c_vec  = vector()
Alobs_10KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Alobs_10KmDispDist_predict_mat = matrix(nrow = length(Alobs_10KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Alobs_10KmDispDist = gbm.step(data=Alobs_10KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Alobs_10KmDispDist)$var, imp = summary(gbm_mod_Alobs_10KmDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Alobs_10KmDispDist_HSI_vec = append(Alobs_10KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Alobs_10KmDispDist_EgoSize_vec = append(Alobs_10KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Alobs_10KmDispDist_strength_vec  = append(Alobs_10KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Alobs_10KmDispDist_deg_vec  = append(Alobs_10KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Alobs_10KmDispDist_habAv_vec  = append(Alobs_10KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Alobs_10KmDispDist_unw_b_c_vec  = append(Alobs_10KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Alobs_10KmDispDist_Patch_Area_vec  = append(Alobs_10KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Alobs_10KmDispDist_AUC_cv_vec = append(Alobs_10KmDispDist_AUC_cv_vec, gbm_mod_Alobs_10KmDispDist$cv.statistics$discrimination.mean)
  Alobs_10KmDispDist_AUC_vec = append(Alobs_10KmDispDist_AUC_vec, gbm_mod_Alobs_10KmDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Alobs_10KmDispDist$n.trees
  Alobs_10KmDispDist_nt_vec = append(Alobs_10KmDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Alobs_10KmDispDist_predict_mat[,i] = predict(gbm_mod_Alobs_10KmDispDist, Alobs_10KmDispDist_stattest, gbm_mod_Alobs_10KmDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Alobs_10KmDispDist_output_tab[,"AUC_train"] = Alobs_10KmDispDist_AUC_vec
Alobs_10KmDispDist_output_tab[,"AUC_cv"] = Alobs_10KmDispDist_AUC_cv_vec
Alobs_10KmDispDist_output_tab[,"ntrees"] = Alobs_10KmDispDist_nt_vec

#Var. importance columns
Alobs_10KmDispDist_output_tab[,"HSI_imp"] = Alobs_10KmDispDist_HSI_vec
Alobs_10KmDispDist_output_tab[,"EgoSize_imp"] = Alobs_10KmDispDist_EgoSize_vec
Alobs_10KmDispDist_output_tab[,"strength_imp"] = Alobs_10KmDispDist_strength_vec
Alobs_10KmDispDist_output_tab[,"deg_imp"] = Alobs_10KmDispDist_deg_vec
Alobs_10KmDispDist_output_tab[,"habAv_imp"] = Alobs_10KmDispDist_habAv_vec
Alobs_10KmDispDist_output_tab[,"unw_b_c_imp"] = Alobs_10KmDispDist_unw_b_c_vec
Alobs_10KmDispDist_output_tab[,"Patch_Area_imp"] = Alobs_10KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Alobs_10KmDispDist_HSI_tab = data.frame(imp = Alobs_10KmDispDist_HSI_vec)
Alobs_10KmDispDist_EgoSize_tab = data.frame(imp = Alobs_10KmDispDist_EgoSize_vec)
Alobs_10KmDispDist_strength_tab = data.frame(imp = Alobs_10KmDispDist_strength_vec)
Alobs_10KmDispDist_deg_tab = data.frame(imp = Alobs_10KmDispDist_deg_vec)
Alobs_10KmDispDist_habAv_tab = data.frame(imp = Alobs_10KmDispDist_habAv_vec)
Alobs_10KmDispDist_unw_b_c_tab = data.frame(imp = Alobs_10KmDispDist_unw_b_c_vec)
Alobs_10KmDispDist_Patch_Area_tab = data.frame(imp = Alobs_10KmDispDist_Patch_Area_vec)

Alobs_10KmDispDist_HSI_tab$variable = "HSI"
Alobs_10KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Alobs_10KmDispDist_strength_tab$variable = "Strength"
Alobs_10KmDispDist_deg_tab$variable = "Degree"
Alobs_10KmDispDist_habAv_tab$variable = "Hab. Av."
Alobs_10KmDispDist_unw_b_c_tab$variable = "B.C."
Alobs_10KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Alobs_10KmDispDist_AUC_cv_tab = data.frame(value = Alobs_10KmDispDist_AUC_cv_vec)
Alobs_10KmDispDist_AUC_train_tab = data.frame(value = Alobs_10KmDispDist_AUC_vec)
#Make label of model for plot
Alobs_10KmDispDist_AUC_cv_tab$model = "Alobs_10KmDispDist"
Alobs_10KmDispDist_AUC_train_tab$model = "Alobs_10KmDispDist"

#combine pred. vars. into new data frame
Alobs_10KmDispDist_var_imp_tab = rbind(Alobs_10KmDispDist_HSI_tab,Alobs_10KmDispDist_EgoSize_tab,Alobs_10KmDispDist_strength_tab,Alobs_10KmDispDist_habAv_tab,Alobs_10KmDispDist_deg_tab,Alobs_10KmDispDist_unw_b_c_tab,Alobs_10KmDispDist_Patch_Area_tab)

ggplot(Alobs_10KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Alobs_10KmDispDist_var_imp_tab$imp~Alobs_10KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "10 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Alobs_10KmDispDist_HSI_vec)
mean(Alobs_10KmDispDist_EgoSize_vec)
mean(Alobs_10KmDispDist_strength_vec)
mean(Alobs_10KmDispDist_deg_vec)
mean(Alobs_10KmDispDist_habAv_vec)
mean(Alobs_10KmDispDist_unw_b_c_vec)
mean(Alobs_10KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Alobs_10KmDispDist_output_tab$AUC_cv)
summary(Alobs_10KmDispDist_output_tab$AUC_train)


####################################### 1Km ####################################

# Create the output table that contains all the values of each of the runs
Alobs_1KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Alobs

# Create vectors to record performance measures
Alobs_1KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Alobs_1KmDispDist_AUC_vec = vector() #

Alobs_1KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Alobs_1KmDispDist_HSI_vec = vector() 
Alobs_1KmDispDist_EgoSize_vec = vector() 
Alobs_1KmDispDist_strength_vec  = vector()
Alobs_1KmDispDist_deg_vec  = vector()
Alobs_1KmDispDist_habAv_vec  = vector()
Alobs_1KmDispDist_unw_b_c_vec  = vector()
Alobs_1KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Alobs_1KmDispDist_predict_mat = matrix(nrow = length(Alobs_1KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Alobs_1KmDispDist = gbm.step(data=Alobs_1KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Alobs_1KmDispDist)$var, imp = summary(gbm_mod_Alobs_1KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Alobs_1KmDispDist_HSI_vec = append(Alobs_1KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Alobs_1KmDispDist_EgoSize_vec = append(Alobs_1KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Alobs_1KmDispDist_strength_vec  = append(Alobs_1KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Alobs_1KmDispDist_deg_vec  = append(Alobs_1KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Alobs_1KmDispDist_habAv_vec  = append(Alobs_1KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Alobs_1KmDispDist_unw_b_c_vec  = append(Alobs_1KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Alobs_1KmDispDist_Patch_Area_vec  = append(Alobs_1KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Alobs_1KmDispDist_AUC_cv_vec = append(Alobs_1KmDispDist_AUC_cv_vec, gbm_mod_Alobs_1KmDispDist$cv.statistics$discrimination.mean)
  Alobs_1KmDispDist_AUC_vec = append(Alobs_1KmDispDist_AUC_vec, gbm_mod_Alobs_1KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Alobs_1KmDispDist$n.trees
  Alobs_1KmDispDist_nt_vec = append(Alobs_1KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Alobs_1KmDispDist_predict_mat[,i] = predict(gbm_mod_Alobs_1KmDispDist, Alobs_1KmDispDist_stattest, gbm_mod_Alobs_1KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Alobs_1KmDispDist_output_tab[,"AUC_train"] = Alobs_1KmDispDist_AUC_vec
Alobs_1KmDispDist_output_tab[,"AUC_cv"] = Alobs_1KmDispDist_AUC_cv_vec
Alobs_1KmDispDist_output_tab[,"ntrees"] = Alobs_1KmDispDist_nt_vec

#Var. importance columns
Alobs_1KmDispDist_output_tab[,"HSI_imp"] = Alobs_1KmDispDist_HSI_vec
Alobs_1KmDispDist_output_tab[,"EgoSize_imp"] = Alobs_1KmDispDist_EgoSize_vec
Alobs_1KmDispDist_output_tab[,"strength_imp"] = Alobs_1KmDispDist_strength_vec
Alobs_1KmDispDist_output_tab[,"deg_imp"] = Alobs_1KmDispDist_deg_vec
Alobs_1KmDispDist_output_tab[,"habAv_imp"] = Alobs_1KmDispDist_habAv_vec
Alobs_1KmDispDist_output_tab[,"unw_b_c_imp"] = Alobs_1KmDispDist_unw_b_c_vec
Alobs_1KmDispDist_output_tab[,"Patch_Area_imp"] = Alobs_1KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Alobs_1KmDispDist_HSI_tab = data.frame(imp = Alobs_1KmDispDist_HSI_vec)
Alobs_1KmDispDist_EgoSize_tab = data.frame(imp = Alobs_1KmDispDist_EgoSize_vec)
Alobs_1KmDispDist_strength_tab = data.frame(imp = Alobs_1KmDispDist_strength_vec)
Alobs_1KmDispDist_deg_tab = data.frame(imp = Alobs_1KmDispDist_deg_vec)
Alobs_1KmDispDist_habAv_tab = data.frame(imp = Alobs_1KmDispDist_habAv_vec)
Alobs_1KmDispDist_unw_b_c_tab = data.frame(imp = Alobs_1KmDispDist_unw_b_c_vec)
Alobs_1KmDispDist_Patch_Area_tab = data.frame(imp = Alobs_1KmDispDist_Patch_Area_vec)

Alobs_1KmDispDist_HSI_tab$variable = "HSI"
Alobs_1KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Alobs_1KmDispDist_strength_tab$variable = "Strength"
Alobs_1KmDispDist_deg_tab$variable = "Degree"
Alobs_1KmDispDist_habAv_tab$variable = "Hab. Av."
Alobs_1KmDispDist_unw_b_c_tab$variable = "B.C."
Alobs_1KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Alobs_1KmDispDist_AUC_cv_tab = data.frame(value = Alobs_1KmDispDist_AUC_cv_vec)
Alobs_1KmDispDist_AUC_train_tab = data.frame(value = Alobs_1KmDispDist_AUC_vec)
#Make label of model for plot
Alobs_1KmDispDist_AUC_cv_tab$model = "Alobs_1KmDispDist"
Alobs_1KmDispDist_AUC_train_tab$model = "Alobs_1KmDispDist"

#combine pred. vars. into new data frame 
Alobs_1KmDispDist_var_imp_tab = rbind(Alobs_1KmDispDist_HSI_tab,Alobs_1KmDispDist_EgoSize_tab,Alobs_1KmDispDist_strength_tab,Alobs_1KmDispDist_habAv_tab,Alobs_1KmDispDist_deg_tab,Alobs_1KmDispDist_unw_b_c_tab,Alobs_1KmDispDist_Patch_Area_tab)

ggplot(Alobs_1KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Alobs_1KmDispDist_var_imp_tab$imp~Alobs_1KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "1 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Alobs_1KmDispDist_HSI_vec)
mean(Alobs_1KmDispDist_EgoSize_vec)
mean(Alobs_1KmDispDist_strength_vec)
mean(Alobs_1KmDispDist_deg_vec)
mean(Alobs_1KmDispDist_habAv_vec)
mean(Alobs_1KmDispDist_unw_b_c_vec)
mean(Alobs_1KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Alobs_1KmDispDist_output_tab$AUC_cv)
summary(Alobs_1KmDispDist_output_tab$AUC_train)


####################################### 6Km ####################################

# Create the output table that contains all the values of each of the runs
Alobs_6KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Alobs

# Create vectors to record performance measures
Alobs_6KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Alobs_6KmDispDist_AUC_vec = vector() #

Alobs_6KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Alobs_6KmDispDist_HSI_vec = vector() 
Alobs_6KmDispDist_EgoSize_vec = vector() 
Alobs_6KmDispDist_strength_vec  = vector()
Alobs_6KmDispDist_deg_vec  = vector()
Alobs_6KmDispDist_habAv_vec  = vector()
Alobs_6KmDispDist_unw_b_c_vec  = vector()
Alobs_6KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Alobs_6KmDispDist_predict_mat = matrix(nrow = length(Alobs_6KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Alobs_6KmDispDist = gbm.step(data=Alobs_6KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Alobs_6KmDispDist)$var, imp = summary(gbm_mod_Alobs_6KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Alobs_6KmDispDist_HSI_vec = append(Alobs_6KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Alobs_6KmDispDist_EgoSize_vec = append(Alobs_6KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Alobs_6KmDispDist_strength_vec  = append(Alobs_6KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Alobs_6KmDispDist_deg_vec  = append(Alobs_6KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Alobs_6KmDispDist_habAv_vec  = append(Alobs_6KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Alobs_6KmDispDist_unw_b_c_vec  = append(Alobs_6KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Alobs_6KmDispDist_Patch_Area_vec  = append(Alobs_6KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Alobs_6KmDispDist_AUC_cv_vec = append(Alobs_6KmDispDist_AUC_cv_vec, gbm_mod_Alobs_6KmDispDist$cv.statistics$discrimination.mean)
  Alobs_6KmDispDist_AUC_vec = append(Alobs_6KmDispDist_AUC_vec, gbm_mod_Alobs_6KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Alobs_6KmDispDist$n.trees
  Alobs_6KmDispDist_nt_vec = append(Alobs_6KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Alobs_6KmDispDist_predict_mat[,i] = predict(gbm_mod_Alobs_6KmDispDist, Alobs_6KmDispDist_stattest, gbm_mod_Alobs_6KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Alobs_6KmDispDist_output_tab[,"AUC_train"] = Alobs_6KmDispDist_AUC_vec
Alobs_6KmDispDist_output_tab[,"AUC_cv"] = Alobs_6KmDispDist_AUC_cv_vec
Alobs_6KmDispDist_output_tab[,"ntrees"] = Alobs_6KmDispDist_nt_vec

#Var. importance columns
Alobs_6KmDispDist_output_tab[,"HSI_imp"] = Alobs_6KmDispDist_HSI_vec
Alobs_6KmDispDist_output_tab[,"EgoSize_imp"] = Alobs_6KmDispDist_EgoSize_vec
Alobs_6KmDispDist_output_tab[,"strength_imp"] = Alobs_6KmDispDist_strength_vec
Alobs_6KmDispDist_output_tab[,"deg_imp"] = Alobs_6KmDispDist_deg_vec
Alobs_6KmDispDist_output_tab[,"habAv_imp"] = Alobs_6KmDispDist_habAv_vec
Alobs_6KmDispDist_output_tab[,"unw_b_c_imp"] = Alobs_6KmDispDist_unw_b_c_vec
Alobs_6KmDispDist_output_tab[,"Patch_Area_imp"] = Alobs_6KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Alobs_6KmDispDist_HSI_tab = data.frame(imp = Alobs_6KmDispDist_HSI_vec)
Alobs_6KmDispDist_EgoSize_tab = data.frame(imp = Alobs_6KmDispDist_EgoSize_vec)
Alobs_6KmDispDist_strength_tab = data.frame(imp = Alobs_6KmDispDist_strength_vec)
Alobs_6KmDispDist_deg_tab = data.frame(imp = Alobs_6KmDispDist_deg_vec)
Alobs_6KmDispDist_habAv_tab = data.frame(imp = Alobs_6KmDispDist_habAv_vec)
Alobs_6KmDispDist_unw_b_c_tab = data.frame(imp = Alobs_6KmDispDist_unw_b_c_vec)
Alobs_6KmDispDist_Patch_Area_tab = data.frame(imp = Alobs_6KmDispDist_Patch_Area_vec)

Alobs_6KmDispDist_HSI_tab$variable = "HSI"
Alobs_6KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Alobs_6KmDispDist_strength_tab$variable = "Strength"
Alobs_6KmDispDist_deg_tab$variable = "Degree"
Alobs_6KmDispDist_habAv_tab$variable = "Hab. Av."
Alobs_6KmDispDist_unw_b_c_tab$variable = "B.C."
Alobs_6KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Alobs_6KmDispDist_AUC_cv_tab = data.frame(value = Alobs_6KmDispDist_AUC_cv_vec)
Alobs_6KmDispDist_AUC_train_tab = data.frame(value = Alobs_6KmDispDist_AUC_vec)
#Make label of model for plot
Alobs_6KmDispDist_AUC_cv_tab$model = "Alobs_6KmDispDist"
Alobs_6KmDispDist_AUC_train_tab$model = "Alobs_6KmDispDist"

#combine pred. vars. into new data frame 
Alobs_6KmDispDist_var_imp_tab = rbind(Alobs_6KmDispDist_HSI_tab,Alobs_6KmDispDist_EgoSize_tab,Alobs_6KmDispDist_strength_tab,Alobs_6KmDispDist_habAv_tab,Alobs_6KmDispDist_deg_tab,Alobs_6KmDispDist_unw_b_c_tab,Alobs_6KmDispDist_Patch_Area_tab)

ggplot(Alobs_6KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Alobs_6KmDispDist_var_imp_tab$imp~Alobs_6KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "6 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Alobs_6KmDispDist_HSI_vec)
mean(Alobs_6KmDispDist_EgoSize_vec)
mean(Alobs_6KmDispDist_strength_vec)
mean(Alobs_6KmDispDist_deg_vec)
mean(Alobs_6KmDispDist_habAv_vec)
mean(Alobs_6KmDispDist_unw_b_c_vec)
mean(Alobs_6KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Alobs_6KmDispDist_output_tab$AUC_cv)
summary(Alobs_6KmDispDist_output_tab$AUC_train)


####################################### 4Km ####################################

# Create the output table that contains all the values of each of the runs
Alobs_4KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Alobs

# Create vectors to record performance measures
Alobs_4KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Alobs_4KmDispDist_AUC_vec = vector() #

Alobs_4KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Alobs_4KmDispDist_HSI_vec = vector()
Alobs_4KmDispDist_EgoSize_vec = vector()
Alobs_4KmDispDist_strength_vec  = vector()
Alobs_4KmDispDist_deg_vec  = vector()
Alobs_4KmDispDist_habAv_vec  = vector()
Alobs_4KmDispDist_unw_b_c_vec  = vector()
Alobs_4KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Alobs_4KmDispDist_predict_mat = matrix(nrow = length(Alobs_4KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Alobs_4KmDispDist = gbm.step(data=Alobs_4KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Alobs_4KmDispDist)$var, imp = summary(gbm_mod_Alobs_4KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Alobs_4KmDispDist_HSI_vec = append(Alobs_4KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Alobs_4KmDispDist_EgoSize_vec = append(Alobs_4KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Alobs_4KmDispDist_strength_vec  = append(Alobs_4KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Alobs_4KmDispDist_deg_vec  = append(Alobs_4KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Alobs_4KmDispDist_habAv_vec  = append(Alobs_4KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Alobs_4KmDispDist_unw_b_c_vec  = append(Alobs_4KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Alobs_4KmDispDist_Patch_Area_vec  = append(Alobs_4KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Alobs_4KmDispDist_AUC_cv_vec = append(Alobs_4KmDispDist_AUC_cv_vec, gbm_mod_Alobs_4KmDispDist$cv.statistics$discrimination.mean)
  Alobs_4KmDispDist_AUC_vec = append(Alobs_4KmDispDist_AUC_vec, gbm_mod_Alobs_4KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Alobs_4KmDispDist$n.trees
  Alobs_4KmDispDist_nt_vec = append(Alobs_4KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Alobs_4KmDispDist_predict_mat[,i] = predict(gbm_mod_Alobs_4KmDispDist, Alobs_4KmDispDist_stattest, gbm_mod_Alobs_4KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Alobs_4KmDispDist_output_tab[,"AUC_train"] = Alobs_4KmDispDist_AUC_vec
Alobs_4KmDispDist_output_tab[,"AUC_cv"] = Alobs_4KmDispDist_AUC_cv_vec
Alobs_4KmDispDist_output_tab[,"ntrees"] = Alobs_4KmDispDist_nt_vec

#Var. importance columns
Alobs_4KmDispDist_output_tab[,"HSI_imp"] = Alobs_4KmDispDist_HSI_vec
Alobs_4KmDispDist_output_tab[,"EgoSize_imp"] = Alobs_4KmDispDist_EgoSize_vec
Alobs_4KmDispDist_output_tab[,"strength_imp"] = Alobs_4KmDispDist_strength_vec
Alobs_4KmDispDist_output_tab[,"deg_imp"] = Alobs_4KmDispDist_deg_vec
Alobs_4KmDispDist_output_tab[,"habAv_imp"] = Alobs_4KmDispDist_habAv_vec
Alobs_4KmDispDist_output_tab[,"unw_b_c_imp"] = Alobs_4KmDispDist_unw_b_c_vec
Alobs_4KmDispDist_output_tab[,"Patch_Area_imp"] = Alobs_4KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Alobs_4KmDispDist_HSI_tab = data.frame(imp = Alobs_4KmDispDist_HSI_vec)
Alobs_4KmDispDist_EgoSize_tab = data.frame(imp = Alobs_4KmDispDist_EgoSize_vec)
Alobs_4KmDispDist_strength_tab = data.frame(imp = Alobs_4KmDispDist_strength_vec)
Alobs_4KmDispDist_deg_tab = data.frame(imp = Alobs_4KmDispDist_deg_vec)
Alobs_4KmDispDist_habAv_tab = data.frame(imp = Alobs_4KmDispDist_habAv_vec)
Alobs_4KmDispDist_unw_b_c_tab = data.frame(imp = Alobs_4KmDispDist_unw_b_c_vec)
Alobs_4KmDispDist_Patch_Area_tab = data.frame(imp = Alobs_4KmDispDist_Patch_Area_vec)

Alobs_4KmDispDist_HSI_tab$variable = "HSI"
Alobs_4KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Alobs_4KmDispDist_strength_tab$variable = "Strength"
Alobs_4KmDispDist_deg_tab$variable = "Degree"
Alobs_4KmDispDist_habAv_tab$variable = "Hab. Av."
Alobs_4KmDispDist_unw_b_c_tab$variable = "B.C."
Alobs_4KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Alobs_4KmDispDist_AUC_cv_tab = data.frame(value = Alobs_4KmDispDist_AUC_cv_vec)
Alobs_4KmDispDist_AUC_train_tab = data.frame(value = Alobs_4KmDispDist_AUC_vec)
#Make label of model for plot
Alobs_4KmDispDist_AUC_cv_tab$model = "Alobs_4KmDispDist"
Alobs_4KmDispDist_AUC_train_tab$model = "Alobs_4KmDispDist"

#combine pred. vars. into new data frame
Alobs_4KmDispDist_var_imp_tab = rbind(Alobs_4KmDispDist_HSI_tab,Alobs_4KmDispDist_EgoSize_tab,Alobs_4KmDispDist_strength_tab,Alobs_4KmDispDist_habAv_tab,Alobs_4KmDispDist_deg_tab,Alobs_4KmDispDist_unw_b_c_tab,Alobs_4KmDispDist_Patch_Area_tab)

ggplot(Alobs_4KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Alobs_4KmDispDist_var_imp_tab$imp~Alobs_4KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "4 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Alobs_4KmDispDist_HSI_vec)
mean(Alobs_4KmDispDist_EgoSize_vec)
mean(Alobs_4KmDispDist_strength_vec)
mean(Alobs_4KmDispDist_deg_vec)
mean(Alobs_4KmDispDist_habAv_vec)
mean(Alobs_4KmDispDist_unw_b_c_vec)
mean(Alobs_4KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Alobs_4KmDispDist_output_tab$AUC_cv)
summary(Alobs_4KmDispDist_output_tab$AUC_train)


####################################### 8Km ####################################

# Create the output table that contains all the values of each of the runs
Alobs_8KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Alobs

# Create vectors to record performance measures
Alobs_8KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Alobs_8KmDispDist_AUC_vec = vector() #

Alobs_8KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Alobs_8KmDispDist_HSI_vec = vector() 
Alobs_8KmDispDist_EgoSize_vec = vector() 
Alobs_8KmDispDist_strength_vec  = vector()
Alobs_8KmDispDist_deg_vec  = vector()
Alobs_8KmDispDist_habAv_vec  = vector()
Alobs_8KmDispDist_unw_b_c_vec  = vector()
Alobs_8KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Alobs_8KmDispDist_predict_mat = matrix(nrow = length(Alobs_8KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Alobs_8KmDispDist = gbm.step(data=Alobs_8KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Alobs_8KmDispDist)$var, imp = summary(gbm_mod_Alobs_8KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Alobs_8KmDispDist_HSI_vec = append(Alobs_8KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Alobs_8KmDispDist_EgoSize_vec = append(Alobs_8KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Alobs_8KmDispDist_strength_vec  = append(Alobs_8KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Alobs_8KmDispDist_deg_vec  = append(Alobs_8KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Alobs_8KmDispDist_habAv_vec  = append(Alobs_8KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Alobs_8KmDispDist_unw_b_c_vec  = append(Alobs_8KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Alobs_8KmDispDist_Patch_Area_vec  = append(Alobs_8KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Alobs_8KmDispDist_AUC_cv_vec = append(Alobs_8KmDispDist_AUC_cv_vec, gbm_mod_Alobs_8KmDispDist$cv.statistics$discrimination.mean)
  Alobs_8KmDispDist_AUC_vec = append(Alobs_8KmDispDist_AUC_vec, gbm_mod_Alobs_8KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Alobs_8KmDispDist$n.trees
  Alobs_8KmDispDist_nt_vec = append(Alobs_8KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Alobs_8KmDispDist_predict_mat[,i] = predict(gbm_mod_Alobs_8KmDispDist, Alobs_8KmDispDist_stattest, gbm_mod_Alobs_8KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Alobs_8KmDispDist_output_tab[,"AUC_train"] = Alobs_8KmDispDist_AUC_vec
Alobs_8KmDispDist_output_tab[,"AUC_cv"] = Alobs_8KmDispDist_AUC_cv_vec
Alobs_8KmDispDist_output_tab[,"ntrees"] = Alobs_8KmDispDist_nt_vec

#Var. importance columns
Alobs_8KmDispDist_output_tab[,"HSI_imp"] = Alobs_8KmDispDist_HSI_vec
Alobs_8KmDispDist_output_tab[,"EgoSize_imp"] = Alobs_8KmDispDist_EgoSize_vec
Alobs_8KmDispDist_output_tab[,"strength_imp"] = Alobs_8KmDispDist_strength_vec
Alobs_8KmDispDist_output_tab[,"deg_imp"] = Alobs_8KmDispDist_deg_vec
Alobs_8KmDispDist_output_tab[,"habAv_imp"] = Alobs_8KmDispDist_habAv_vec
Alobs_8KmDispDist_output_tab[,"unw_b_c_imp"] = Alobs_8KmDispDist_unw_b_c_vec
Alobs_8KmDispDist_output_tab[,"Patch_Area_imp"] = Alobs_8KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Alobs_8KmDispDist_HSI_tab = data.frame(imp = Alobs_8KmDispDist_HSI_vec)
Alobs_8KmDispDist_EgoSize_tab = data.frame(imp = Alobs_8KmDispDist_EgoSize_vec)
Alobs_8KmDispDist_strength_tab = data.frame(imp = Alobs_8KmDispDist_strength_vec)
Alobs_8KmDispDist_deg_tab = data.frame(imp = Alobs_8KmDispDist_deg_vec)
Alobs_8KmDispDist_habAv_tab = data.frame(imp = Alobs_8KmDispDist_habAv_vec)
Alobs_8KmDispDist_unw_b_c_tab = data.frame(imp = Alobs_8KmDispDist_unw_b_c_vec)
Alobs_8KmDispDist_Patch_Area_tab = data.frame(imp = Alobs_8KmDispDist_Patch_Area_vec)

Alobs_8KmDispDist_HSI_tab$variable = "HSI"
Alobs_8KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Alobs_8KmDispDist_strength_tab$variable = "Strength"
Alobs_8KmDispDist_deg_tab$variable = "Degree"
Alobs_8KmDispDist_habAv_tab$variable = "Hab. Av."
Alobs_8KmDispDist_unw_b_c_tab$variable = "B.C."
Alobs_8KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Alobs_8KmDispDist_AUC_cv_tab = data.frame(value = Alobs_8KmDispDist_AUC_cv_vec)
Alobs_8KmDispDist_AUC_train_tab = data.frame(value = Alobs_8KmDispDist_AUC_vec)
#Make label of model for plot
Alobs_8KmDispDist_AUC_cv_tab$model = "Alobs_8KmDispDist"
Alobs_8KmDispDist_AUC_train_tab$model = "Alobs_8KmDispDist"

#combine pred. vars. into new data frame 
Alobs_8KmDispDist_var_imp_tab = rbind(Alobs_8KmDispDist_HSI_tab,Alobs_8KmDispDist_EgoSize_tab,Alobs_8KmDispDist_strength_tab,Alobs_8KmDispDist_habAv_tab,Alobs_8KmDispDist_deg_tab,Alobs_8KmDispDist_unw_b_c_tab,Alobs_8KmDispDist_Patch_Area_tab)

ggplot(Alobs_8KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Alobs_8KmDispDist_var_imp_tab$imp~Alobs_8KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "8 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Alobs_8KmDispDist_HSI_vec)
mean(Alobs_8KmDispDist_EgoSize_vec)
mean(Alobs_8KmDispDist_strength_vec)
mean(Alobs_8KmDispDist_deg_vec)
mean(Alobs_8KmDispDist_habAv_vec)
mean(Alobs_8KmDispDist_unw_b_c_vec)
mean(Alobs_8KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Alobs_8KmDispDist_output_tab$AUC_cv)
summary(Alobs_8KmDispDist_output_tab$AUC_train)


####################################################################################################
#### Compare scores between networks w/d0 variations of the same species ###########################
### Alobs ##################

#Import evaluation df's of original run with species-specific dispersal distance
Alobs_DefaultDispDist_AUC_cv_tab <- read.csv("Alobs_DefaultDispDist_AUC_cv_tab.csv")
Alobs_DefaultDispDist_AUC_train_tab <- read.csv("Alobs_DefaultDispDist_AUC_train_tab.csv")
Alobs_DefaultDispDist_var_imp_tab <- read.csv("Alobs_DefaultDispDist_var_imp_tab.csv")
Alobs_DefaultDispDist_output_tab <- read.csv("Alobs_DefaultDispDist_BRToutput_tab.csv")
Alobs_DefaultDispDist_predict_df <- read.csv("Alobs_DefaultDispDist_predict_df.csv")

head(Alobs_DefaultDispDist_AUC_cv_tab)
Alobs_DefaultDispDist_AUC_cv_tab$X <- NULL
Alobs_DefaultDispDist_AUC_train_tab$X <- NULL

#### Make a dataframe for plotting overlaying histograms in R
### cv AUC
Alobs_AUC_cv_tab = rbind(Alobs_DefaultDispDist_AUC_cv_tab, Alobs_300mDispDist_AUC_cv_tab, Alobs_1KmDispDist_AUC_cv_tab, Alobs_2KmDispDist_AUC_cv_tab, 
                         Alobs_4KmDispDist_AUC_cv_tab,
                         Alobs_6KmDispDist_AUC_cv_tab, Alobs_8KmDispDist_AUC_cv_tab, Alobs_10KmDispDist_AUC_cv_tab)
#Change order of factors to display noTopo at the edge
Alobs_AUC_cv_tab$model<- as.factor(Alobs_AUC_cv_tab$model)
levels(Alobs_AUC_cv_tab$model)
Alobs_AUC_cv_tab$model<-factor(Alobs_AUC_cv_tab$model, levels=c("Alobs_300mDispDist", "Alobs_1KmDispDist", "Alobs", 
                                                                "Alobs_2KmDispDist", "Alobs_4KmDispDist", 
                                                                "Alobs_6KmDispDist",  "Alobs_8KmDispDist", "Alobs_10KmDispDist"))
#Plot
ggplot(Alobs_AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Alobs_AUC_cv_tab$value~Alobs_AUC_cv_tab$model, ylab= "Cross-validated AUC")


### train AUC
Alobs_AUC_train_tab = rbind(Alobs_DefaultDispDist_AUC_train_tab, Alobs_300mDispDist_AUC_train_tab, Alobs_1KmDispDist_AUC_train_tab, Alobs_2KmDispDist_AUC_train_tab, 
                            Alobs_4KmDispDist_AUC_train_tab,
                            Alobs_6KmDispDist_AUC_train_tab, Alobs_8KmDispDist_AUC_train_tab, Alobs_10KmDispDist_AUC_train_tab)

levels(Alobs_AUC_train_tab$model)
Alobs_AUC_train_tab$model<-factor(Alobs_AUC_train_tab$model, levels=c("Alobs_300mDispDist", "Alobs_1KmDispDist", "Alobs", 
                                                                      "Alobs_2KmDispDist", "Alobs_4KmDispDist", 
                                                                      "Alobs_6KmDispDist",  "Alobs_8KmDispDist", "Alobs_10KmDispDist"))

ggplot(Alobs_AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Alobs_AUC_train_tab$value~Alobs_AUC_train_tab$model, ylab= "Training AUC")

summary(Alobs_DefaultDispDist_output_tab$AUC_cv)
summary(Alobs_300mDispDist_output_tab$AUC_cv)
summary(Alobs_1KmDispDist_output_tab$AUC_cv)
summary(Alobs_2KmDispDist_output_tab$AUC_cv)
summary(Alobs_4KmDispDist_output_tab$AUC_cv)
summary(Alobs_6KmDispDist_output_tab$AUC_cv)
summary(Alobs_8KmDispDist_output_tab$AUC_cv)
summary(Alobs_10KmDispDist_output_tab$AUC_cv)

summary(Alobs_DefaultDispDist_output_tab$AUC_train)
summary(Alobs_300mDispDist_output_tab$AUC_train)
summary(Alobs_1KmDispDist_output_tab$AUC_train)
summary(Alobs_2KmDispDist_output_tab$AUC_train)
summary(Alobs_4KmDispDist_output_tab$AUC_train)
summary(Alobs_6KmDispDist_output_tab$AUC_train)
summary(Alobs_8KmDispDist_output_tab$AUC_train)
summary(Alobs_10KmDispDist_output_tab$AUC_train)

summary(Alobs_DefaultDispDist_output_tab$ntrees)
summary(Alobs_300mDispDist_output_tab$ntrees)
summary(Alobs_1KmDispDist_output_tab$ntrees)
summary(Alobs_2KmDispDist_output_tab$ntrees)
summary(Alobs_4KmDispDist_output_tab$ntrees)
summary(Alobs_6KmDispDist_output_tab$ntrees)
summary(Alobs_8KmDispDist_output_tab$ntrees)
summary(Alobs_10KmDispDist_output_tab$ntrees)



################################################################################
################ Epcal #########################################################

#### 2Km #######################################################################

# Create the output table that contains all the values of each of the runs
Epcal_2KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal

# Create vectors to record performance measures
Epcal_2KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Epcal_2KmDispDist_AUC_vec = vector() #

Epcal_2KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Epcal_2KmDispDist_HSI_vec = vector()
Epcal_2KmDispDist_EgoSize_vec = vector()
Epcal_2KmDispDist_strength_vec  = vector()
Epcal_2KmDispDist_deg_vec  = vector()
Epcal_2KmDispDist_habAv_vec  = vector()
Epcal_2KmDispDist_unw_b_c_vec  = vector()
Epcal_2KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Epcal_2KmDispDist_predict_mat = matrix(nrow = length(Epcal_2KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Epcal_2KmDispDist = gbm.step(data=Epcal_2KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Epcal_2KmDispDist)$var, imp = summary(gbm_mod_Epcal_2KmDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Epcal_2KmDispDist_HSI_vec = append(Epcal_2KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Epcal_2KmDispDist_EgoSize_vec = append(Epcal_2KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Epcal_2KmDispDist_strength_vec  = append(Epcal_2KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Epcal_2KmDispDist_deg_vec  = append(Epcal_2KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Epcal_2KmDispDist_habAv_vec  = append(Epcal_2KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Epcal_2KmDispDist_unw_b_c_vec  = append(Epcal_2KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Epcal_2KmDispDist_Patch_Area_vec  = append(Epcal_2KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Epcal_2KmDispDist_AUC_cv_vec = append(Epcal_2KmDispDist_AUC_cv_vec, gbm_mod_Epcal_2KmDispDist$cv.statistics$discrimination.mean)
  Epcal_2KmDispDist_AUC_vec = append(Epcal_2KmDispDist_AUC_vec, gbm_mod_Epcal_2KmDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Epcal_2KmDispDist$n.trees
  Epcal_2KmDispDist_nt_vec = append(Epcal_2KmDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Epcal_2KmDispDist_predict_mat[,i] = predict(gbm_mod_Epcal_2KmDispDist, Epcal_2KmDispDist_stattest, gbm_mod_Epcal_2KmDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Epcal_2KmDispDist_output_tab[,"AUC_train"] = Epcal_2KmDispDist_AUC_vec
Epcal_2KmDispDist_output_tab[,"AUC_cv"] = Epcal_2KmDispDist_AUC_cv_vec
Epcal_2KmDispDist_output_tab[,"ntrees"] = Epcal_2KmDispDist_nt_vec

#Var. importance columns
Epcal_2KmDispDist_output_tab[,"HSI_imp"] = Epcal_2KmDispDist_HSI_vec
Epcal_2KmDispDist_output_tab[,"EgoSize_imp"] = Epcal_2KmDispDist_EgoSize_vec
Epcal_2KmDispDist_output_tab[,"strength_imp"] = Epcal_2KmDispDist_strength_vec
Epcal_2KmDispDist_output_tab[,"deg_imp"] = Epcal_2KmDispDist_deg_vec
Epcal_2KmDispDist_output_tab[,"habAv_imp"] = Epcal_2KmDispDist_habAv_vec
Epcal_2KmDispDist_output_tab[,"unw_b_c_imp"] = Epcal_2KmDispDist_unw_b_c_vec
Epcal_2KmDispDist_output_tab[,"Patch_Area_imp"] = Epcal_2KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Epcal_2KmDispDist_HSI_tab = data.frame(imp = Epcal_2KmDispDist_HSI_vec)
Epcal_2KmDispDist_EgoSize_tab = data.frame(imp = Epcal_2KmDispDist_EgoSize_vec)
Epcal_2KmDispDist_strength_tab = data.frame(imp = Epcal_2KmDispDist_strength_vec)
Epcal_2KmDispDist_deg_tab = data.frame(imp = Epcal_2KmDispDist_deg_vec)
Epcal_2KmDispDist_habAv_tab = data.frame(imp = Epcal_2KmDispDist_habAv_vec)
Epcal_2KmDispDist_unw_b_c_tab = data.frame(imp = Epcal_2KmDispDist_unw_b_c_vec)
Epcal_2KmDispDist_Patch_Area_tab = data.frame(imp = Epcal_2KmDispDist_Patch_Area_vec)

Epcal_2KmDispDist_HSI_tab$variable = "HSI"
Epcal_2KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Epcal_2KmDispDist_strength_tab$variable = "Strength"
Epcal_2KmDispDist_deg_tab$variable = "Degree"
Epcal_2KmDispDist_habAv_tab$variable = "Hab. Av."
Epcal_2KmDispDist_unw_b_c_tab$variable = "B.C."
Epcal_2KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Epcal_2KmDispDist_AUC_cv_tab = data.frame(value = Epcal_2KmDispDist_AUC_cv_vec)
Epcal_2KmDispDist_AUC_train_tab = data.frame(value = Epcal_2KmDispDist_AUC_vec)
#Make label of model for plot
Epcal_2KmDispDist_AUC_cv_tab$model = "Epcal_2KmDispDist"
Epcal_2KmDispDist_AUC_train_tab$model = "Epcal_2KmDispDist"

#combine pred. vars. into new data frame
Epcal_2KmDispDist_var_imp_tab = rbind(Epcal_2KmDispDist_HSI_tab,Epcal_2KmDispDist_EgoSize_tab,Epcal_2KmDispDist_strength_tab,Epcal_2KmDispDist_habAv_tab,Epcal_2KmDispDist_deg_tab,Epcal_2KmDispDist_unw_b_c_tab,Epcal_2KmDispDist_Patch_Area_tab)

ggplot(Epcal_2KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Epcal_2KmDispDist_var_imp_tab$imp~Epcal_2KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "2 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)


#Get mean var. importance of all of the vars.
mean(Epcal_2KmDispDist_HSI_vec)
mean(Epcal_2KmDispDist_EgoSize_vec)
mean(Epcal_2KmDispDist_strength_vec)
mean(Epcal_2KmDispDist_deg_vec)
mean(Epcal_2KmDispDist_habAv_vec)
mean(Epcal_2KmDispDist_unw_b_c_vec)
mean(Epcal_2KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Epcal_2KmDispDist_output_tab$AUC_cv)
summary(Epcal_2KmDispDist_output_tab$AUC_train)


###### 300m ####################################################################
 
# Create the output table that contains all the values of each of the runs
Epcal_300mDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal

# Create vectors to record performance measures
Epcal_300mDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Epcal_300mDispDist_AUC_vec = vector() #

Epcal_300mDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Epcal_300mDispDist_HSI_vec = vector()
Epcal_300mDispDist_EgoSize_vec = vector()
Epcal_300mDispDist_strength_vec  = vector()
Epcal_300mDispDist_deg_vec  = vector()
Epcal_300mDispDist_habAv_vec  = vector()
Epcal_300mDispDist_unw_b_c_vec  = vector()
Epcal_300mDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Epcal_300mDispDist_predict_mat = matrix(nrow = length(Epcal_300mDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Epcal_300mDispDist = gbm.step(data=Epcal_300mDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Epcal_300mDispDist)$var, imp = summary(gbm_mod_Epcal_300mDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Epcal_300mDispDist_HSI_vec = append(Epcal_300mDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Epcal_300mDispDist_EgoSize_vec = append(Epcal_300mDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Epcal_300mDispDist_strength_vec  = append(Epcal_300mDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Epcal_300mDispDist_deg_vec  = append(Epcal_300mDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Epcal_300mDispDist_habAv_vec  = append(Epcal_300mDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Epcal_300mDispDist_unw_b_c_vec  = append(Epcal_300mDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Epcal_300mDispDist_Patch_Area_vec  = append(Epcal_300mDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Epcal_300mDispDist_AUC_cv_vec = append(Epcal_300mDispDist_AUC_cv_vec, gbm_mod_Epcal_300mDispDist$cv.statistics$discrimination.mean)
  Epcal_300mDispDist_AUC_vec = append(Epcal_300mDispDist_AUC_vec, gbm_mod_Epcal_300mDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Epcal_300mDispDist$n.trees
  Epcal_300mDispDist_nt_vec = append(Epcal_300mDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Epcal_300mDispDist_predict_mat[,i] = predict(gbm_mod_Epcal_300mDispDist, Epcal_300mDispDist_stattest, gbm_mod_Epcal_300mDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Epcal_300mDispDist_output_tab[,"AUC_train"] = Epcal_300mDispDist_AUC_vec
Epcal_300mDispDist_output_tab[,"AUC_cv"] = Epcal_300mDispDist_AUC_cv_vec
Epcal_300mDispDist_output_tab[,"ntrees"] = Epcal_300mDispDist_nt_vec

#Var. importance columns
Epcal_300mDispDist_output_tab[,"HSI_imp"] = Epcal_300mDispDist_HSI_vec
Epcal_300mDispDist_output_tab[,"EgoSize_imp"] = Epcal_300mDispDist_EgoSize_vec
Epcal_300mDispDist_output_tab[,"strength_imp"] = Epcal_300mDispDist_strength_vec
Epcal_300mDispDist_output_tab[,"deg_imp"] = Epcal_300mDispDist_deg_vec
Epcal_300mDispDist_output_tab[,"habAv_imp"] = Epcal_300mDispDist_habAv_vec
Epcal_300mDispDist_output_tab[,"unw_b_c_imp"] = Epcal_300mDispDist_unw_b_c_vec
Epcal_300mDispDist_output_tab[,"Patch_Area_imp"] = Epcal_300mDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Epcal_300mDispDist_HSI_tab = data.frame(imp = Epcal_300mDispDist_HSI_vec)
Epcal_300mDispDist_EgoSize_tab = data.frame(imp = Epcal_300mDispDist_EgoSize_vec)
Epcal_300mDispDist_strength_tab = data.frame(imp = Epcal_300mDispDist_strength_vec)
Epcal_300mDispDist_deg_tab = data.frame(imp = Epcal_300mDispDist_deg_vec)
Epcal_300mDispDist_habAv_tab = data.frame(imp = Epcal_300mDispDist_habAv_vec)
Epcal_300mDispDist_unw_b_c_tab = data.frame(imp = Epcal_300mDispDist_unw_b_c_vec)
Epcal_300mDispDist_Patch_Area_tab = data.frame(imp = Epcal_300mDispDist_Patch_Area_vec)

Epcal_300mDispDist_HSI_tab$variable = "HSI"
Epcal_300mDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Epcal_300mDispDist_strength_tab$variable = "Strength"
Epcal_300mDispDist_deg_tab$variable = "Degree"
Epcal_300mDispDist_habAv_tab$variable = "Hab. Av."
Epcal_300mDispDist_unw_b_c_tab$variable = "B.C."
Epcal_300mDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Epcal_300mDispDist_AUC_cv_tab = data.frame(value = Epcal_300mDispDist_AUC_cv_vec)
Epcal_300mDispDist_AUC_train_tab = data.frame(value = Epcal_300mDispDist_AUC_vec)
#Make label of model for plot
Epcal_300mDispDist_AUC_cv_tab$model = "Epcal_300mDispDist"
Epcal_300mDispDist_AUC_train_tab$model = "Epcal_300mDispDist"

#combine pred. vars. into new data frame
Epcal_300mDispDist_var_imp_tab = rbind(Epcal_300mDispDist_HSI_tab,Epcal_300mDispDist_EgoSize_tab,Epcal_300mDispDist_strength_tab,Epcal_300mDispDist_habAv_tab,Epcal_300mDispDist_deg_tab,Epcal_300mDispDist_unw_b_c_tab,Epcal_300mDispDist_Patch_Area_tab)

ggplot(Epcal_300mDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Epcal_300mDispDist_var_imp_tab$imp~Epcal_300mDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "300 m maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Epcal_300mDispDist_HSI_vec)
mean(Epcal_300mDispDist_EgoSize_vec)
mean(Epcal_300mDispDist_strength_vec)
mean(Epcal_300mDispDist_deg_vec)
mean(Epcal_300mDispDist_habAv_vec)
mean(Epcal_300mDispDist_unw_b_c_vec)
mean(Epcal_300mDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Epcal_300mDispDist_output_tab$AUC_cv)
summary(Epcal_300mDispDist_output_tab$AUC_train)


#### 10Km ######################################################################

# Create the output table that contains all the values of each of the runs
Epcal_10KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal

# Create vectors to record performance measures
Epcal_10KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Epcal_10KmDispDist_AUC_vec = vector() #

Epcal_10KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Epcal_10KmDispDist_HSI_vec = vector()
Epcal_10KmDispDist_EgoSize_vec = vector()
Epcal_10KmDispDist_strength_vec  = vector()
Epcal_10KmDispDist_deg_vec  = vector()
Epcal_10KmDispDist_habAv_vec  = vector()
Epcal_10KmDispDist_unw_b_c_vec  = vector()
Epcal_10KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Epcal_10KmDispDist_predict_mat = matrix(nrow = length(Epcal_10KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Epcal_10KmDispDist = gbm.step(data=Epcal_10KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Epcal_10KmDispDist)$var, imp = summary(gbm_mod_Epcal_10KmDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Epcal_10KmDispDist_HSI_vec = append(Epcal_10KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Epcal_10KmDispDist_EgoSize_vec = append(Epcal_10KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Epcal_10KmDispDist_strength_vec  = append(Epcal_10KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Epcal_10KmDispDist_deg_vec  = append(Epcal_10KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Epcal_10KmDispDist_habAv_vec  = append(Epcal_10KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Epcal_10KmDispDist_unw_b_c_vec  = append(Epcal_10KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Epcal_10KmDispDist_Patch_Area_vec  = append(Epcal_10KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Epcal_10KmDispDist_AUC_cv_vec = append(Epcal_10KmDispDist_AUC_cv_vec, gbm_mod_Epcal_10KmDispDist$cv.statistics$discrimination.mean)
  Epcal_10KmDispDist_AUC_vec = append(Epcal_10KmDispDist_AUC_vec, gbm_mod_Epcal_10KmDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Epcal_10KmDispDist$n.trees
  Epcal_10KmDispDist_nt_vec = append(Epcal_10KmDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Epcal_10KmDispDist_predict_mat[,i] = predict(gbm_mod_Epcal_10KmDispDist, Epcal_10KmDispDist_stattest, gbm_mod_Epcal_10KmDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Epcal_10KmDispDist_output_tab[,"AUC_train"] = Epcal_10KmDispDist_AUC_vec
Epcal_10KmDispDist_output_tab[,"AUC_cv"] = Epcal_10KmDispDist_AUC_cv_vec
Epcal_10KmDispDist_output_tab[,"ntrees"] = Epcal_10KmDispDist_nt_vec

#Var. importance columns
Epcal_10KmDispDist_output_tab[,"HSI_imp"] = Epcal_10KmDispDist_HSI_vec
Epcal_10KmDispDist_output_tab[,"EgoSize_imp"] = Epcal_10KmDispDist_EgoSize_vec
Epcal_10KmDispDist_output_tab[,"strength_imp"] = Epcal_10KmDispDist_strength_vec
Epcal_10KmDispDist_output_tab[,"deg_imp"] = Epcal_10KmDispDist_deg_vec
Epcal_10KmDispDist_output_tab[,"habAv_imp"] = Epcal_10KmDispDist_habAv_vec
Epcal_10KmDispDist_output_tab[,"unw_b_c_imp"] = Epcal_10KmDispDist_unw_b_c_vec
Epcal_10KmDispDist_output_tab[,"Patch_Area_imp"] = Epcal_10KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Epcal_10KmDispDist_HSI_tab = data.frame(imp = Epcal_10KmDispDist_HSI_vec)
Epcal_10KmDispDist_EgoSize_tab = data.frame(imp = Epcal_10KmDispDist_EgoSize_vec)
Epcal_10KmDispDist_strength_tab = data.frame(imp = Epcal_10KmDispDist_strength_vec)
Epcal_10KmDispDist_deg_tab = data.frame(imp = Epcal_10KmDispDist_deg_vec)
Epcal_10KmDispDist_habAv_tab = data.frame(imp = Epcal_10KmDispDist_habAv_vec)
Epcal_10KmDispDist_unw_b_c_tab = data.frame(imp = Epcal_10KmDispDist_unw_b_c_vec)
Epcal_10KmDispDist_Patch_Area_tab = data.frame(imp = Epcal_10KmDispDist_Patch_Area_vec)

Epcal_10KmDispDist_HSI_tab$variable = "HSI"
Epcal_10KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Epcal_10KmDispDist_strength_tab$variable = "Strength"
Epcal_10KmDispDist_deg_tab$variable = "Degree"
Epcal_10KmDispDist_habAv_tab$variable = "Hab. Av."
Epcal_10KmDispDist_unw_b_c_tab$variable = "B.C."
Epcal_10KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Epcal_10KmDispDist_AUC_cv_tab = data.frame(value = Epcal_10KmDispDist_AUC_cv_vec)
Epcal_10KmDispDist_AUC_train_tab = data.frame(value = Epcal_10KmDispDist_AUC_vec)
#Make label of model for plot
Epcal_10KmDispDist_AUC_cv_tab$model = "Epcal_10KmDispDist"
Epcal_10KmDispDist_AUC_train_tab$model = "Epcal_10KmDispDist"

#combine pred. vars. into new data frame
Epcal_10KmDispDist_var_imp_tab = rbind(Epcal_10KmDispDist_HSI_tab,Epcal_10KmDispDist_EgoSize_tab,Epcal_10KmDispDist_strength_tab,Epcal_10KmDispDist_habAv_tab,Epcal_10KmDispDist_deg_tab,Epcal_10KmDispDist_unw_b_c_tab,Epcal_10KmDispDist_Patch_Area_tab)

ggplot(Epcal_10KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Epcal_10KmDispDist_var_imp_tab$imp~Epcal_10KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "10 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Epcal_10KmDispDist_HSI_vec)
mean(Epcal_10KmDispDist_EgoSize_vec)
mean(Epcal_10KmDispDist_strength_vec)
mean(Epcal_10KmDispDist_deg_vec)
mean(Epcal_10KmDispDist_habAv_vec)
mean(Epcal_10KmDispDist_unw_b_c_vec)
mean(Epcal_10KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Epcal_10KmDispDist_output_tab$AUC_cv)
summary(Epcal_10KmDispDist_output_tab$AUC_train)


####################################### 1Km ####################################

# Create the output table that contains all the values of each of the runs
Epcal_1KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal

# Create vectors to record performance measures
Epcal_1KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Epcal_1KmDispDist_AUC_vec = vector() #

Epcal_1KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Epcal_1KmDispDist_HSI_vec = vector() 
Epcal_1KmDispDist_EgoSize_vec = vector() 
Epcal_1KmDispDist_strength_vec  = vector()
Epcal_1KmDispDist_deg_vec  = vector()
Epcal_1KmDispDist_habAv_vec  = vector()
Epcal_1KmDispDist_unw_b_c_vec  = vector()
Epcal_1KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Epcal_1KmDispDist_predict_mat = matrix(nrow = length(Epcal_1KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Epcal_1KmDispDist = gbm.step(data=Epcal_1KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Epcal_1KmDispDist)$var, imp = summary(gbm_mod_Epcal_1KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Epcal_1KmDispDist_HSI_vec = append(Epcal_1KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Epcal_1KmDispDist_EgoSize_vec = append(Epcal_1KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Epcal_1KmDispDist_strength_vec  = append(Epcal_1KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Epcal_1KmDispDist_deg_vec  = append(Epcal_1KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Epcal_1KmDispDist_habAv_vec  = append(Epcal_1KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Epcal_1KmDispDist_unw_b_c_vec  = append(Epcal_1KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Epcal_1KmDispDist_Patch_Area_vec  = append(Epcal_1KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Epcal_1KmDispDist_AUC_cv_vec = append(Epcal_1KmDispDist_AUC_cv_vec, gbm_mod_Epcal_1KmDispDist$cv.statistics$discrimination.mean)
  Epcal_1KmDispDist_AUC_vec = append(Epcal_1KmDispDist_AUC_vec, gbm_mod_Epcal_1KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Epcal_1KmDispDist$n.trees
  Epcal_1KmDispDist_nt_vec = append(Epcal_1KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Epcal_1KmDispDist_predict_mat[,i] = predict(gbm_mod_Epcal_1KmDispDist, Epcal_1KmDispDist_stattest, gbm_mod_Epcal_1KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Epcal_1KmDispDist_output_tab[,"AUC_train"] = Epcal_1KmDispDist_AUC_vec
Epcal_1KmDispDist_output_tab[,"AUC_cv"] = Epcal_1KmDispDist_AUC_cv_vec
Epcal_1KmDispDist_output_tab[,"ntrees"] = Epcal_1KmDispDist_nt_vec

#Var. importance columns
Epcal_1KmDispDist_output_tab[,"HSI_imp"] = Epcal_1KmDispDist_HSI_vec
Epcal_1KmDispDist_output_tab[,"EgoSize_imp"] = Epcal_1KmDispDist_EgoSize_vec
Epcal_1KmDispDist_output_tab[,"strength_imp"] = Epcal_1KmDispDist_strength_vec
Epcal_1KmDispDist_output_tab[,"deg_imp"] = Epcal_1KmDispDist_deg_vec
Epcal_1KmDispDist_output_tab[,"habAv_imp"] = Epcal_1KmDispDist_habAv_vec
Epcal_1KmDispDist_output_tab[,"unw_b_c_imp"] = Epcal_1KmDispDist_unw_b_c_vec
Epcal_1KmDispDist_output_tab[,"Patch_Area_imp"] = Epcal_1KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Epcal_1KmDispDist_HSI_tab = data.frame(imp = Epcal_1KmDispDist_HSI_vec)
Epcal_1KmDispDist_EgoSize_tab = data.frame(imp = Epcal_1KmDispDist_EgoSize_vec)
Epcal_1KmDispDist_strength_tab = data.frame(imp = Epcal_1KmDispDist_strength_vec)
Epcal_1KmDispDist_deg_tab = data.frame(imp = Epcal_1KmDispDist_deg_vec)
Epcal_1KmDispDist_habAv_tab = data.frame(imp = Epcal_1KmDispDist_habAv_vec)
Epcal_1KmDispDist_unw_b_c_tab = data.frame(imp = Epcal_1KmDispDist_unw_b_c_vec)
Epcal_1KmDispDist_Patch_Area_tab = data.frame(imp = Epcal_1KmDispDist_Patch_Area_vec)

Epcal_1KmDispDist_HSI_tab$variable = "HSI"
Epcal_1KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Epcal_1KmDispDist_strength_tab$variable = "Strength"
Epcal_1KmDispDist_deg_tab$variable = "Degree"
Epcal_1KmDispDist_habAv_tab$variable = "Hab. Av."
Epcal_1KmDispDist_unw_b_c_tab$variable = "B.C."
Epcal_1KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Epcal_1KmDispDist_AUC_cv_tab = data.frame(value = Epcal_1KmDispDist_AUC_cv_vec)
Epcal_1KmDispDist_AUC_train_tab = data.frame(value = Epcal_1KmDispDist_AUC_vec)
#Make label of model for plot
Epcal_1KmDispDist_AUC_cv_tab$model = "Epcal_1KmDispDist"
Epcal_1KmDispDist_AUC_train_tab$model = "Epcal_1KmDispDist"

#combine pred. vars. into new data frame 
Epcal_1KmDispDist_var_imp_tab = rbind(Epcal_1KmDispDist_HSI_tab,Epcal_1KmDispDist_EgoSize_tab,Epcal_1KmDispDist_strength_tab,Epcal_1KmDispDist_habAv_tab,Epcal_1KmDispDist_deg_tab,Epcal_1KmDispDist_unw_b_c_tab,Epcal_1KmDispDist_Patch_Area_tab)

ggplot(Epcal_1KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Epcal_1KmDispDist_var_imp_tab$imp~Epcal_1KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "1 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Epcal_1KmDispDist_HSI_vec)
mean(Epcal_1KmDispDist_EgoSize_vec)
mean(Epcal_1KmDispDist_strength_vec)
mean(Epcal_1KmDispDist_deg_vec)
mean(Epcal_1KmDispDist_habAv_vec)
mean(Epcal_1KmDispDist_unw_b_c_vec)
mean(Epcal_1KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Epcal_1KmDispDist_output_tab$AUC_cv)
summary(Epcal_1KmDispDist_output_tab$AUC_train)


####################################### 6Km ####################################

# Create the output table that contains all the values of each of the runs
Epcal_6KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal

# Create vectors to record performance measures
Epcal_6KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Epcal_6KmDispDist_AUC_vec = vector() #

Epcal_6KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Epcal_6KmDispDist_HSI_vec = vector() 
Epcal_6KmDispDist_EgoSize_vec = vector() 
Epcal_6KmDispDist_strength_vec  = vector()
Epcal_6KmDispDist_deg_vec  = vector()
Epcal_6KmDispDist_habAv_vec  = vector()
Epcal_6KmDispDist_unw_b_c_vec  = vector()
Epcal_6KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Epcal_6KmDispDist_predict_mat = matrix(nrow = length(Epcal_6KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Epcal_6KmDispDist = gbm.step(data=Epcal_6KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Epcal_6KmDispDist)$var, imp = summary(gbm_mod_Epcal_6KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Epcal_6KmDispDist_HSI_vec = append(Epcal_6KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Epcal_6KmDispDist_EgoSize_vec = append(Epcal_6KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Epcal_6KmDispDist_strength_vec  = append(Epcal_6KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Epcal_6KmDispDist_deg_vec  = append(Epcal_6KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Epcal_6KmDispDist_habAv_vec  = append(Epcal_6KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Epcal_6KmDispDist_unw_b_c_vec  = append(Epcal_6KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Epcal_6KmDispDist_Patch_Area_vec  = append(Epcal_6KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Epcal_6KmDispDist_AUC_cv_vec = append(Epcal_6KmDispDist_AUC_cv_vec, gbm_mod_Epcal_6KmDispDist$cv.statistics$discrimination.mean)
  Epcal_6KmDispDist_AUC_vec = append(Epcal_6KmDispDist_AUC_vec, gbm_mod_Epcal_6KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Epcal_6KmDispDist$n.trees
  Epcal_6KmDispDist_nt_vec = append(Epcal_6KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Epcal_6KmDispDist_predict_mat[,i] = predict(gbm_mod_Epcal_6KmDispDist, Epcal_6KmDispDist_stattest, gbm_mod_Epcal_6KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Epcal_6KmDispDist_output_tab[,"AUC_train"] = Epcal_6KmDispDist_AUC_vec
Epcal_6KmDispDist_output_tab[,"AUC_cv"] = Epcal_6KmDispDist_AUC_cv_vec
Epcal_6KmDispDist_output_tab[,"ntrees"] = Epcal_6KmDispDist_nt_vec

#Var. importance columns
Epcal_6KmDispDist_output_tab[,"HSI_imp"] = Epcal_6KmDispDist_HSI_vec
Epcal_6KmDispDist_output_tab[,"EgoSize_imp"] = Epcal_6KmDispDist_EgoSize_vec
Epcal_6KmDispDist_output_tab[,"strength_imp"] = Epcal_6KmDispDist_strength_vec
Epcal_6KmDispDist_output_tab[,"deg_imp"] = Epcal_6KmDispDist_deg_vec
Epcal_6KmDispDist_output_tab[,"habAv_imp"] = Epcal_6KmDispDist_habAv_vec
Epcal_6KmDispDist_output_tab[,"unw_b_c_imp"] = Epcal_6KmDispDist_unw_b_c_vec
Epcal_6KmDispDist_output_tab[,"Patch_Area_imp"] = Epcal_6KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Epcal_6KmDispDist_HSI_tab = data.frame(imp = Epcal_6KmDispDist_HSI_vec)
Epcal_6KmDispDist_EgoSize_tab = data.frame(imp = Epcal_6KmDispDist_EgoSize_vec)
Epcal_6KmDispDist_strength_tab = data.frame(imp = Epcal_6KmDispDist_strength_vec)
Epcal_6KmDispDist_deg_tab = data.frame(imp = Epcal_6KmDispDist_deg_vec)
Epcal_6KmDispDist_habAv_tab = data.frame(imp = Epcal_6KmDispDist_habAv_vec)
Epcal_6KmDispDist_unw_b_c_tab = data.frame(imp = Epcal_6KmDispDist_unw_b_c_vec)
Epcal_6KmDispDist_Patch_Area_tab = data.frame(imp = Epcal_6KmDispDist_Patch_Area_vec)

Epcal_6KmDispDist_HSI_tab$variable = "HSI"
Epcal_6KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Epcal_6KmDispDist_strength_tab$variable = "Strength"
Epcal_6KmDispDist_deg_tab$variable = "Degree"
Epcal_6KmDispDist_habAv_tab$variable = "Hab. Av."
Epcal_6KmDispDist_unw_b_c_tab$variable = "B.C."
Epcal_6KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Epcal_6KmDispDist_AUC_cv_tab = data.frame(value = Epcal_6KmDispDist_AUC_cv_vec)
Epcal_6KmDispDist_AUC_train_tab = data.frame(value = Epcal_6KmDispDist_AUC_vec)
#Make label of model for plot
Epcal_6KmDispDist_AUC_cv_tab$model = "Epcal_6KmDispDist"
Epcal_6KmDispDist_AUC_train_tab$model = "Epcal_6KmDispDist"

#combine pred. vars. into new data frame 
Epcal_6KmDispDist_var_imp_tab = rbind(Epcal_6KmDispDist_HSI_tab,Epcal_6KmDispDist_EgoSize_tab,Epcal_6KmDispDist_strength_tab,Epcal_6KmDispDist_habAv_tab,Epcal_6KmDispDist_deg_tab,Epcal_6KmDispDist_unw_b_c_tab,Epcal_6KmDispDist_Patch_Area_tab)

ggplot(Epcal_6KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Epcal_6KmDispDist_var_imp_tab$imp~Epcal_6KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "6 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Epcal_6KmDispDist_HSI_vec)
mean(Epcal_6KmDispDist_EgoSize_vec)
mean(Epcal_6KmDispDist_strength_vec)
mean(Epcal_6KmDispDist_deg_vec)
mean(Epcal_6KmDispDist_habAv_vec)
mean(Epcal_6KmDispDist_unw_b_c_vec)
mean(Epcal_6KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Epcal_6KmDispDist_output_tab$AUC_cv)
summary(Epcal_6KmDispDist_output_tab$AUC_train)


####################################### 4Km ####################################

# Create the output table that contains all the values of each of the runs
Epcal_4KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal

# Create vectors to record performance measures
Epcal_4KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Epcal_4KmDispDist_AUC_vec = vector() #

Epcal_4KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Epcal_4KmDispDist_HSI_vec = vector()
Epcal_4KmDispDist_EgoSize_vec = vector()
Epcal_4KmDispDist_strength_vec  = vector()
Epcal_4KmDispDist_deg_vec  = vector()
Epcal_4KmDispDist_habAv_vec  = vector()
Epcal_4KmDispDist_unw_b_c_vec  = vector()
Epcal_4KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Epcal_4KmDispDist_predict_mat = matrix(nrow = length(Epcal_4KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Epcal_4KmDispDist = gbm.step(data=Epcal_4KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Epcal_4KmDispDist)$var, imp = summary(gbm_mod_Epcal_4KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Epcal_4KmDispDist_HSI_vec = append(Epcal_4KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Epcal_4KmDispDist_EgoSize_vec = append(Epcal_4KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Epcal_4KmDispDist_strength_vec  = append(Epcal_4KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Epcal_4KmDispDist_deg_vec  = append(Epcal_4KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Epcal_4KmDispDist_habAv_vec  = append(Epcal_4KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Epcal_4KmDispDist_unw_b_c_vec  = append(Epcal_4KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Epcal_4KmDispDist_Patch_Area_vec  = append(Epcal_4KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Epcal_4KmDispDist_AUC_cv_vec = append(Epcal_4KmDispDist_AUC_cv_vec, gbm_mod_Epcal_4KmDispDist$cv.statistics$discrimination.mean)
  Epcal_4KmDispDist_AUC_vec = append(Epcal_4KmDispDist_AUC_vec, gbm_mod_Epcal_4KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Epcal_4KmDispDist$n.trees
  Epcal_4KmDispDist_nt_vec = append(Epcal_4KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Epcal_4KmDispDist_predict_mat[,i] = predict(gbm_mod_Epcal_4KmDispDist, Epcal_4KmDispDist_stattest, gbm_mod_Epcal_4KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Epcal_4KmDispDist_output_tab[,"AUC_train"] = Epcal_4KmDispDist_AUC_vec
Epcal_4KmDispDist_output_tab[,"AUC_cv"] = Epcal_4KmDispDist_AUC_cv_vec
Epcal_4KmDispDist_output_tab[,"ntrees"] = Epcal_4KmDispDist_nt_vec

#Var. importance columns
Epcal_4KmDispDist_output_tab[,"HSI_imp"] = Epcal_4KmDispDist_HSI_vec
Epcal_4KmDispDist_output_tab[,"EgoSize_imp"] = Epcal_4KmDispDist_EgoSize_vec
Epcal_4KmDispDist_output_tab[,"strength_imp"] = Epcal_4KmDispDist_strength_vec
Epcal_4KmDispDist_output_tab[,"deg_imp"] = Epcal_4KmDispDist_deg_vec
Epcal_4KmDispDist_output_tab[,"habAv_imp"] = Epcal_4KmDispDist_habAv_vec
Epcal_4KmDispDist_output_tab[,"unw_b_c_imp"] = Epcal_4KmDispDist_unw_b_c_vec
Epcal_4KmDispDist_output_tab[,"Patch_Area_imp"] = Epcal_4KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Epcal_4KmDispDist_HSI_tab = data.frame(imp = Epcal_4KmDispDist_HSI_vec)
Epcal_4KmDispDist_EgoSize_tab = data.frame(imp = Epcal_4KmDispDist_EgoSize_vec)
Epcal_4KmDispDist_strength_tab = data.frame(imp = Epcal_4KmDispDist_strength_vec)
Epcal_4KmDispDist_deg_tab = data.frame(imp = Epcal_4KmDispDist_deg_vec)
Epcal_4KmDispDist_habAv_tab = data.frame(imp = Epcal_4KmDispDist_habAv_vec)
Epcal_4KmDispDist_unw_b_c_tab = data.frame(imp = Epcal_4KmDispDist_unw_b_c_vec)
Epcal_4KmDispDist_Patch_Area_tab = data.frame(imp = Epcal_4KmDispDist_Patch_Area_vec)

Epcal_4KmDispDist_HSI_tab$variable = "HSI"
Epcal_4KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Epcal_4KmDispDist_strength_tab$variable = "Strength"
Epcal_4KmDispDist_deg_tab$variable = "Degree"
Epcal_4KmDispDist_habAv_tab$variable = "Hab. Av."
Epcal_4KmDispDist_unw_b_c_tab$variable = "B.C."
Epcal_4KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Epcal_4KmDispDist_AUC_cv_tab = data.frame(value = Epcal_4KmDispDist_AUC_cv_vec)
Epcal_4KmDispDist_AUC_train_tab = data.frame(value = Epcal_4KmDispDist_AUC_vec)
#Make label of model for plot
Epcal_4KmDispDist_AUC_cv_tab$model = "Epcal_4KmDispDist"
Epcal_4KmDispDist_AUC_train_tab$model = "Epcal_4KmDispDist"

#combine pred. vars. into new data frame
Epcal_4KmDispDist_var_imp_tab = rbind(Epcal_4KmDispDist_HSI_tab,Epcal_4KmDispDist_EgoSize_tab,Epcal_4KmDispDist_strength_tab,Epcal_4KmDispDist_habAv_tab,Epcal_4KmDispDist_deg_tab,Epcal_4KmDispDist_unw_b_c_tab,Epcal_4KmDispDist_Patch_Area_tab)

ggplot(Epcal_4KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Epcal_4KmDispDist_var_imp_tab$imp~Epcal_4KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "4 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Epcal_4KmDispDist_HSI_vec)
mean(Epcal_4KmDispDist_EgoSize_vec)
mean(Epcal_4KmDispDist_strength_vec)
mean(Epcal_4KmDispDist_deg_vec)
mean(Epcal_4KmDispDist_habAv_vec)
mean(Epcal_4KmDispDist_unw_b_c_vec)
mean(Epcal_4KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Epcal_4KmDispDist_output_tab$AUC_cv)
summary(Epcal_4KmDispDist_output_tab$AUC_train)


####################################### 8Km ####################################

# Create the output table that contains all the values of each of the runs
Epcal_8KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal

# Create vectors to record performance measures
Epcal_8KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Epcal_8KmDispDist_AUC_vec = vector() #

Epcal_8KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Epcal_8KmDispDist_HSI_vec = vector() 
Epcal_8KmDispDist_EgoSize_vec = vector() 
Epcal_8KmDispDist_strength_vec  = vector()
Epcal_8KmDispDist_deg_vec  = vector()
Epcal_8KmDispDist_habAv_vec  = vector()
Epcal_8KmDispDist_unw_b_c_vec  = vector()
Epcal_8KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Epcal_8KmDispDist_predict_mat = matrix(nrow = length(Epcal_8KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Epcal_8KmDispDist = gbm.step(data=Epcal_8KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Epcal_8KmDispDist)$var, imp = summary(gbm_mod_Epcal_8KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Epcal_8KmDispDist_HSI_vec = append(Epcal_8KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Epcal_8KmDispDist_EgoSize_vec = append(Epcal_8KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Epcal_8KmDispDist_strength_vec  = append(Epcal_8KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Epcal_8KmDispDist_deg_vec  = append(Epcal_8KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Epcal_8KmDispDist_habAv_vec  = append(Epcal_8KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Epcal_8KmDispDist_unw_b_c_vec  = append(Epcal_8KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Epcal_8KmDispDist_Patch_Area_vec  = append(Epcal_8KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Epcal_8KmDispDist_AUC_cv_vec = append(Epcal_8KmDispDist_AUC_cv_vec, gbm_mod_Epcal_8KmDispDist$cv.statistics$discrimination.mean)
  Epcal_8KmDispDist_AUC_vec = append(Epcal_8KmDispDist_AUC_vec, gbm_mod_Epcal_8KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Epcal_8KmDispDist$n.trees
  Epcal_8KmDispDist_nt_vec = append(Epcal_8KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Epcal_8KmDispDist_predict_mat[,i] = predict(gbm_mod_Epcal_8KmDispDist, Epcal_8KmDispDist_stattest, gbm_mod_Epcal_8KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Epcal_8KmDispDist_output_tab[,"AUC_train"] = Epcal_8KmDispDist_AUC_vec
Epcal_8KmDispDist_output_tab[,"AUC_cv"] = Epcal_8KmDispDist_AUC_cv_vec
Epcal_8KmDispDist_output_tab[,"ntrees"] = Epcal_8KmDispDist_nt_vec

#Var. importance columns
Epcal_8KmDispDist_output_tab[,"HSI_imp"] = Epcal_8KmDispDist_HSI_vec
Epcal_8KmDispDist_output_tab[,"EgoSize_imp"] = Epcal_8KmDispDist_EgoSize_vec
Epcal_8KmDispDist_output_tab[,"strength_imp"] = Epcal_8KmDispDist_strength_vec
Epcal_8KmDispDist_output_tab[,"deg_imp"] = Epcal_8KmDispDist_deg_vec
Epcal_8KmDispDist_output_tab[,"habAv_imp"] = Epcal_8KmDispDist_habAv_vec
Epcal_8KmDispDist_output_tab[,"unw_b_c_imp"] = Epcal_8KmDispDist_unw_b_c_vec
Epcal_8KmDispDist_output_tab[,"Patch_Area_imp"] = Epcal_8KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Epcal_8KmDispDist_HSI_tab = data.frame(imp = Epcal_8KmDispDist_HSI_vec)
Epcal_8KmDispDist_EgoSize_tab = data.frame(imp = Epcal_8KmDispDist_EgoSize_vec)
Epcal_8KmDispDist_strength_tab = data.frame(imp = Epcal_8KmDispDist_strength_vec)
Epcal_8KmDispDist_deg_tab = data.frame(imp = Epcal_8KmDispDist_deg_vec)
Epcal_8KmDispDist_habAv_tab = data.frame(imp = Epcal_8KmDispDist_habAv_vec)
Epcal_8KmDispDist_unw_b_c_tab = data.frame(imp = Epcal_8KmDispDist_unw_b_c_vec)
Epcal_8KmDispDist_Patch_Area_tab = data.frame(imp = Epcal_8KmDispDist_Patch_Area_vec)

Epcal_8KmDispDist_HSI_tab$variable = "HSI"
Epcal_8KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Epcal_8KmDispDist_strength_tab$variable = "Strength"
Epcal_8KmDispDist_deg_tab$variable = "Degree"
Epcal_8KmDispDist_habAv_tab$variable = "Hab. Av."
Epcal_8KmDispDist_unw_b_c_tab$variable = "B.C."
Epcal_8KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Epcal_8KmDispDist_AUC_cv_tab = data.frame(value = Epcal_8KmDispDist_AUC_cv_vec)
Epcal_8KmDispDist_AUC_train_tab = data.frame(value = Epcal_8KmDispDist_AUC_vec)
#Make label of model for plot
Epcal_8KmDispDist_AUC_cv_tab$model = "Epcal_8KmDispDist"
Epcal_8KmDispDist_AUC_train_tab$model = "Epcal_8KmDispDist"

#combine pred. vars. into new data frame 
Epcal_8KmDispDist_var_imp_tab = rbind(Epcal_8KmDispDist_HSI_tab,Epcal_8KmDispDist_EgoSize_tab,Epcal_8KmDispDist_strength_tab,Epcal_8KmDispDist_habAv_tab,Epcal_8KmDispDist_deg_tab,Epcal_8KmDispDist_unw_b_c_tab,Epcal_8KmDispDist_Patch_Area_tab)

ggplot(Epcal_8KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Epcal_8KmDispDist_var_imp_tab$imp~Epcal_8KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "8 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Epcal_8KmDispDist_HSI_vec)
mean(Epcal_8KmDispDist_EgoSize_vec)
mean(Epcal_8KmDispDist_strength_vec)
mean(Epcal_8KmDispDist_deg_vec)
mean(Epcal_8KmDispDist_habAv_vec)
mean(Epcal_8KmDispDist_unw_b_c_vec)
mean(Epcal_8KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Epcal_8KmDispDist_output_tab$AUC_cv)
summary(Epcal_8KmDispDist_output_tab$AUC_train)


####################################################################################################
#### Compare scores between networks w/d0 variations of the same species ###########################
### Epcal ##################

#Import evaluation df's of original run with species-specific dispersal distance
Epcal_DefaultDispDist_AUC_cv_tab <- read.csv("Epcal_DefaultDispDist_AUC_cv_tab.csv")
Epcal_DefaultDispDist_AUC_train_tab <- read.csv("Epcal_DefaultDispDist_AUC_train_tab.csv")
Epcal_DefaultDispDist_var_imp_tab <- read.csv("Epcal_DefaultDispDist_var_imp_tab.csv")
Epcal_DefaultDispDist_output_tab <- read.csv("Epcal_DefaultDispDist_BRToutput_tab.csv")
Epcal_DefaultDispDist_predict_df <- read.csv("Epcal_DefaultDispDist_predict_df.csv")

head(Epcal_DefaultDispDist_AUC_cv_tab)
Epcal_DefaultDispDist_AUC_cv_tab$X <- NULL
Epcal_DefaultDispDist_AUC_train_tab$X <- NULL

#### Make a dataframe for plotting overlaying histograms in R
### cv AUC
Epcal_AUC_cv_tab = rbind(Epcal_DefaultDispDist_AUC_cv_tab, Epcal_300mDispDist_AUC_cv_tab, Epcal_1KmDispDist_AUC_cv_tab, Epcal_2KmDispDist_AUC_cv_tab, 
                         Epcal_4KmDispDist_AUC_cv_tab,
                         Epcal_6KmDispDist_AUC_cv_tab, Epcal_8KmDispDist_AUC_cv_tab, Epcal_10KmDispDist_AUC_cv_tab)

#Change order of factors to display noTopo at the edge
Epcal_AUC_cv_tab$model<- as.factor(Epcal_AUC_cv_tab$model)
levels(Epcal_AUC_cv_tab$model)
Epcal_AUC_cv_tab$model<-factor(Epcal_AUC_cv_tab$model, levels=c("Epcal_300mDispDist", "Epcal_1KmDispDist", 
                                                                "Epcal_2KmDispDist", "Epcal_4KmDispDist", "Epcal", "Epcal_6KmDispDist",
                                                                "Epcal_8KmDispDist", "Epcal_10KmDispDist"))
#Plot
ggplot(Epcal_AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Epcal_AUC_cv_tab$value~Epcal_AUC_cv_tab$model, ylab= "Cross-validated AUC")

### train AUC
Epcal_AUC_train_tab = rbind(Epcal_DefaultDispDist_AUC_train_tab, Epcal_300mDispDist_AUC_train_tab, Epcal_1KmDispDist_AUC_train_tab, Epcal_2KmDispDist_AUC_train_tab, 
                            Epcal_4KmDispDist_AUC_train_tab,
                            Epcal_6KmDispDist_AUC_train_tab, Epcal_8KmDispDist_AUC_train_tab, Epcal_10KmDispDist_AUC_train_tab)
levels(Epcal_AUC_train_tab$model)
Epcal_AUC_train_tab$model<-factor(Epcal_AUC_train_tab$model, levels=c("Epcal_300mDispDist", "Epcal_1KmDispDist", 
                                                                      "Epcal_2KmDispDist", "Epcal_4KmDispDist", "Epcal", "Epcal_6KmDispDist",
                                                                      "Epcal_8KmDispDist", "Epcal_10KmDispDist"))

ggplot(Epcal_AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Epcal_AUC_train_tab$value~Epcal_AUC_train_tab$model, ylab= "Training AUC")

summary(Epcal_DefaultDispDist_output_tab$AUC_cv)
summary(Epcal_300mDispDist_output_tab$AUC_cv)
summary(Epcal_1KmDispDist_output_tab$AUC_cv)
summary(Epcal_2KmDispDist_output_tab$AUC_cv)
summary(Epcal_4KmDispDist_output_tab$AUC_cv)
summary(Epcal_6KmDispDist_output_tab$AUC_cv)
summary(Epcal_8KmDispDist_output_tab$AUC_cv)
summary(Epcal_10KmDispDist_output_tab$AUC_cv)

summary(Epcal_DefaultDispDist_output_tab$AUC_train)
summary(Epcal_300mDispDist_output_tab$AUC_train)
summary(Epcal_1KmDispDist_output_tab$AUC_train)
summary(Epcal_2KmDispDist_output_tab$AUC_train)
summary(Epcal_4KmDispDist_output_tab$AUC_train)
summary(Epcal_6KmDispDist_output_tab$AUC_train)
summary(Epcal_8KmDispDist_output_tab$AUC_train)
summary(Epcal_10KmDispDist_output_tab$AUC_train)

summary(Epcal_DefaultDispDist_output_tab$ntrees)
summary(Epcal_300mDispDist_output_tab$ntrees)
summary(Epcal_1KmDispDist_output_tab$ntrees)
summary(Epcal_2KmDispDist_output_tab$ntrees)
summary(Epcal_4KmDispDist_output_tab$ntrees)
summary(Epcal_6KmDispDist_output_tab$ntrees)
summary(Epcal_8KmDispDist_output_tab$ntrees)
summary(Epcal_10KmDispDist_output_tab$ntrees)



################################################################################
################ Peagg #########################################################

#### 2Km #######################################################################

# Create the output table that contains all the values of each of the runs
Peagg_2KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg

# Create vectors to record performance measures
Peagg_2KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Peagg_2KmDispDist_AUC_vec = vector() #

Peagg_2KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Peagg_2KmDispDist_HSI_vec = vector()
Peagg_2KmDispDist_EgoSize_vec = vector()
Peagg_2KmDispDist_strength_vec  = vector()
Peagg_2KmDispDist_deg_vec  = vector()
Peagg_2KmDispDist_habAv_vec  = vector()
Peagg_2KmDispDist_unw_b_c_vec  = vector()
Peagg_2KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Peagg_2KmDispDist_predict_mat = matrix(nrow = length(Peagg_2KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Peagg_2KmDispDist = gbm.step(data=Peagg_2KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Peagg_2KmDispDist)$var, imp = summary(gbm_mod_Peagg_2KmDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Peagg_2KmDispDist_HSI_vec = append(Peagg_2KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Peagg_2KmDispDist_EgoSize_vec = append(Peagg_2KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Peagg_2KmDispDist_strength_vec  = append(Peagg_2KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Peagg_2KmDispDist_deg_vec  = append(Peagg_2KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Peagg_2KmDispDist_habAv_vec  = append(Peagg_2KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Peagg_2KmDispDist_unw_b_c_vec  = append(Peagg_2KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Peagg_2KmDispDist_Patch_Area_vec  = append(Peagg_2KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Peagg_2KmDispDist_AUC_cv_vec = append(Peagg_2KmDispDist_AUC_cv_vec, gbm_mod_Peagg_2KmDispDist$cv.statistics$discrimination.mean)
  Peagg_2KmDispDist_AUC_vec = append(Peagg_2KmDispDist_AUC_vec, gbm_mod_Peagg_2KmDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Peagg_2KmDispDist$n.trees
  Peagg_2KmDispDist_nt_vec = append(Peagg_2KmDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Peagg_2KmDispDist_predict_mat[,i] = predict(gbm_mod_Peagg_2KmDispDist, Peagg_2KmDispDist_stattest, gbm_mod_Peagg_2KmDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Peagg_2KmDispDist_output_tab[,"AUC_train"] = Peagg_2KmDispDist_AUC_vec
Peagg_2KmDispDist_output_tab[,"AUC_cv"] = Peagg_2KmDispDist_AUC_cv_vec
Peagg_2KmDispDist_output_tab[,"ntrees"] = Peagg_2KmDispDist_nt_vec

#Var. importance columns
Peagg_2KmDispDist_output_tab[,"HSI_imp"] = Peagg_2KmDispDist_HSI_vec
Peagg_2KmDispDist_output_tab[,"EgoSize_imp"] = Peagg_2KmDispDist_EgoSize_vec
Peagg_2KmDispDist_output_tab[,"strength_imp"] = Peagg_2KmDispDist_strength_vec
Peagg_2KmDispDist_output_tab[,"deg_imp"] = Peagg_2KmDispDist_deg_vec
Peagg_2KmDispDist_output_tab[,"habAv_imp"] = Peagg_2KmDispDist_habAv_vec
Peagg_2KmDispDist_output_tab[,"unw_b_c_imp"] = Peagg_2KmDispDist_unw_b_c_vec
Peagg_2KmDispDist_output_tab[,"Patch_Area_imp"] = Peagg_2KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Peagg_2KmDispDist_HSI_tab = data.frame(imp = Peagg_2KmDispDist_HSI_vec)
Peagg_2KmDispDist_EgoSize_tab = data.frame(imp = Peagg_2KmDispDist_EgoSize_vec)
Peagg_2KmDispDist_strength_tab = data.frame(imp = Peagg_2KmDispDist_strength_vec)
Peagg_2KmDispDist_deg_tab = data.frame(imp = Peagg_2KmDispDist_deg_vec)
Peagg_2KmDispDist_habAv_tab = data.frame(imp = Peagg_2KmDispDist_habAv_vec)
Peagg_2KmDispDist_unw_b_c_tab = data.frame(imp = Peagg_2KmDispDist_unw_b_c_vec)
Peagg_2KmDispDist_Patch_Area_tab = data.frame(imp = Peagg_2KmDispDist_Patch_Area_vec)

Peagg_2KmDispDist_HSI_tab$variable = "HSI"
Peagg_2KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Peagg_2KmDispDist_strength_tab$variable = "Strength"
Peagg_2KmDispDist_deg_tab$variable = "Degree"
Peagg_2KmDispDist_habAv_tab$variable = "Hab. Av."
Peagg_2KmDispDist_unw_b_c_tab$variable = "B.C."
Peagg_2KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Peagg_2KmDispDist_AUC_cv_tab = data.frame(value = Peagg_2KmDispDist_AUC_cv_vec)
Peagg_2KmDispDist_AUC_train_tab = data.frame(value = Peagg_2KmDispDist_AUC_vec)
#Make label of model for plot
Peagg_2KmDispDist_AUC_cv_tab$model = "Peagg_2KmDispDist"
Peagg_2KmDispDist_AUC_train_tab$model = "Peagg_2KmDispDist"

#combine pred. vars. into new data frame
Peagg_2KmDispDist_var_imp_tab = rbind(Peagg_2KmDispDist_HSI_tab,Peagg_2KmDispDist_EgoSize_tab,Peagg_2KmDispDist_strength_tab,Peagg_2KmDispDist_habAv_tab,Peagg_2KmDispDist_deg_tab,Peagg_2KmDispDist_unw_b_c_tab,Peagg_2KmDispDist_Patch_Area_tab)

ggplot(Peagg_2KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Peagg_2KmDispDist_var_imp_tab$imp~Peagg_2KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "2 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Peagg_2KmDispDist_HSI_vec)
mean(Peagg_2KmDispDist_EgoSize_vec)
mean(Peagg_2KmDispDist_strength_vec)
mean(Peagg_2KmDispDist_deg_vec)
mean(Peagg_2KmDispDist_habAv_vec)
mean(Peagg_2KmDispDist_unw_b_c_vec)
mean(Peagg_2KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Peagg_2KmDispDist_output_tab$AUC_cv)
summary(Peagg_2KmDispDist_output_tab$AUC_train)


###### 300m ####################################################################
 
# Create the output table that contains all the values of each of the runs
Peagg_300mDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg

# Create vectors to record performance measures
Peagg_300mDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Peagg_300mDispDist_AUC_vec = vector() #

Peagg_300mDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Peagg_300mDispDist_HSI_vec = vector()
Peagg_300mDispDist_EgoSize_vec = vector()
Peagg_300mDispDist_strength_vec  = vector()
Peagg_300mDispDist_deg_vec  = vector()
Peagg_300mDispDist_habAv_vec  = vector()
Peagg_300mDispDist_unw_b_c_vec  = vector()
Peagg_300mDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Peagg_300mDispDist_predict_mat = matrix(nrow = length(Peagg_300mDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Peagg_300mDispDist = gbm.step(data=Peagg_300mDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Peagg_300mDispDist)$var, imp = summary(gbm_mod_Peagg_300mDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Peagg_300mDispDist_HSI_vec = append(Peagg_300mDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Peagg_300mDispDist_EgoSize_vec = append(Peagg_300mDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Peagg_300mDispDist_strength_vec  = append(Peagg_300mDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Peagg_300mDispDist_deg_vec  = append(Peagg_300mDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Peagg_300mDispDist_habAv_vec  = append(Peagg_300mDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Peagg_300mDispDist_unw_b_c_vec  = append(Peagg_300mDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Peagg_300mDispDist_Patch_Area_vec  = append(Peagg_300mDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Peagg_300mDispDist_AUC_cv_vec = append(Peagg_300mDispDist_AUC_cv_vec, gbm_mod_Peagg_300mDispDist$cv.statistics$discrimination.mean)
  Peagg_300mDispDist_AUC_vec = append(Peagg_300mDispDist_AUC_vec, gbm_mod_Peagg_300mDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Peagg_300mDispDist$n.trees
  Peagg_300mDispDist_nt_vec = append(Peagg_300mDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Peagg_300mDispDist_predict_mat[,i] = predict(gbm_mod_Peagg_300mDispDist, Peagg_300mDispDist_stattest, gbm_mod_Peagg_300mDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Peagg_300mDispDist_output_tab[,"AUC_train"] = Peagg_300mDispDist_AUC_vec
Peagg_300mDispDist_output_tab[,"AUC_cv"] = Peagg_300mDispDist_AUC_cv_vec
Peagg_300mDispDist_output_tab[,"ntrees"] = Peagg_300mDispDist_nt_vec

#Var. importance columns
Peagg_300mDispDist_output_tab[,"HSI_imp"] = Peagg_300mDispDist_HSI_vec
Peagg_300mDispDist_output_tab[,"EgoSize_imp"] = Peagg_300mDispDist_EgoSize_vec
Peagg_300mDispDist_output_tab[,"strength_imp"] = Peagg_300mDispDist_strength_vec
Peagg_300mDispDist_output_tab[,"deg_imp"] = Peagg_300mDispDist_deg_vec
Peagg_300mDispDist_output_tab[,"habAv_imp"] = Peagg_300mDispDist_habAv_vec
Peagg_300mDispDist_output_tab[,"unw_b_c_imp"] = Peagg_300mDispDist_unw_b_c_vec
Peagg_300mDispDist_output_tab[,"Patch_Area_imp"] = Peagg_300mDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Peagg_300mDispDist_HSI_tab = data.frame(imp = Peagg_300mDispDist_HSI_vec)
Peagg_300mDispDist_EgoSize_tab = data.frame(imp = Peagg_300mDispDist_EgoSize_vec)
Peagg_300mDispDist_strength_tab = data.frame(imp = Peagg_300mDispDist_strength_vec)
Peagg_300mDispDist_deg_tab = data.frame(imp = Peagg_300mDispDist_deg_vec)
Peagg_300mDispDist_habAv_tab = data.frame(imp = Peagg_300mDispDist_habAv_vec)
Peagg_300mDispDist_unw_b_c_tab = data.frame(imp = Peagg_300mDispDist_unw_b_c_vec)
Peagg_300mDispDist_Patch_Area_tab = data.frame(imp = Peagg_300mDispDist_Patch_Area_vec)

Peagg_300mDispDist_HSI_tab$variable = "HSI"
Peagg_300mDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Peagg_300mDispDist_strength_tab$variable = "Strength"
Peagg_300mDispDist_deg_tab$variable = "Degree"
Peagg_300mDispDist_habAv_tab$variable = "Hab. Av."
Peagg_300mDispDist_unw_b_c_tab$variable = "B.C."
Peagg_300mDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Peagg_300mDispDist_AUC_cv_tab = data.frame(value = Peagg_300mDispDist_AUC_cv_vec)
Peagg_300mDispDist_AUC_train_tab = data.frame(value = Peagg_300mDispDist_AUC_vec)
#Make label of model for plot
Peagg_300mDispDist_AUC_cv_tab$model = "Peagg_300mDispDist"
Peagg_300mDispDist_AUC_train_tab$model = "Peagg_300mDispDist"

#combine pred. vars. into new data frame
Peagg_300mDispDist_var_imp_tab = rbind(Peagg_300mDispDist_HSI_tab,Peagg_300mDispDist_EgoSize_tab,Peagg_300mDispDist_strength_tab,Peagg_300mDispDist_habAv_tab,Peagg_300mDispDist_deg_tab,Peagg_300mDispDist_unw_b_c_tab,Peagg_300mDispDist_Patch_Area_tab)

ggplot(Peagg_300mDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Peagg_300mDispDist_var_imp_tab$imp~Peagg_300mDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "300 m maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Peagg_300mDispDist_HSI_vec)
mean(Peagg_300mDispDist_EgoSize_vec)
mean(Peagg_300mDispDist_strength_vec)
mean(Peagg_300mDispDist_deg_vec)
mean(Peagg_300mDispDist_habAv_vec)
mean(Peagg_300mDispDist_unw_b_c_vec)
mean(Peagg_300mDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Peagg_300mDispDist_output_tab$AUC_cv)
summary(Peagg_300mDispDist_output_tab$AUC_train)


#### 10Km ######################################################################

# Create the output table that contains all the values of each of the runs
Peagg_10KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg

# Create vectors to record performance measures
Peagg_10KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Peagg_10KmDispDist_AUC_vec = vector() #

Peagg_10KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Peagg_10KmDispDist_HSI_vec = vector()
Peagg_10KmDispDist_EgoSize_vec = vector()
Peagg_10KmDispDist_strength_vec  = vector()
Peagg_10KmDispDist_deg_vec  = vector()
Peagg_10KmDispDist_habAv_vec  = vector()
Peagg_10KmDispDist_unw_b_c_vec  = vector()
Peagg_10KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Peagg_10KmDispDist_predict_mat = matrix(nrow = length(Peagg_10KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Peagg_10KmDispDist = gbm.step(data=Peagg_10KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Peagg_10KmDispDist)$var, imp = summary(gbm_mod_Peagg_10KmDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Peagg_10KmDispDist_HSI_vec = append(Peagg_10KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Peagg_10KmDispDist_EgoSize_vec = append(Peagg_10KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Peagg_10KmDispDist_strength_vec  = append(Peagg_10KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Peagg_10KmDispDist_deg_vec  = append(Peagg_10KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Peagg_10KmDispDist_habAv_vec  = append(Peagg_10KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Peagg_10KmDispDist_unw_b_c_vec  = append(Peagg_10KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Peagg_10KmDispDist_Patch_Area_vec  = append(Peagg_10KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Peagg_10KmDispDist_AUC_cv_vec = append(Peagg_10KmDispDist_AUC_cv_vec, gbm_mod_Peagg_10KmDispDist$cv.statistics$discrimination.mean)
  Peagg_10KmDispDist_AUC_vec = append(Peagg_10KmDispDist_AUC_vec, gbm_mod_Peagg_10KmDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Peagg_10KmDispDist$n.trees
  Peagg_10KmDispDist_nt_vec = append(Peagg_10KmDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Peagg_10KmDispDist_predict_mat[,i] = predict(gbm_mod_Peagg_10KmDispDist, Peagg_10KmDispDist_stattest, gbm_mod_Peagg_10KmDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Peagg_10KmDispDist_output_tab[,"AUC_train"] = Peagg_10KmDispDist_AUC_vec
Peagg_10KmDispDist_output_tab[,"AUC_cv"] = Peagg_10KmDispDist_AUC_cv_vec
Peagg_10KmDispDist_output_tab[,"ntrees"] = Peagg_10KmDispDist_nt_vec

#Var. importance columns
Peagg_10KmDispDist_output_tab[,"HSI_imp"] = Peagg_10KmDispDist_HSI_vec
Peagg_10KmDispDist_output_tab[,"EgoSize_imp"] = Peagg_10KmDispDist_EgoSize_vec
Peagg_10KmDispDist_output_tab[,"strength_imp"] = Peagg_10KmDispDist_strength_vec
Peagg_10KmDispDist_output_tab[,"deg_imp"] = Peagg_10KmDispDist_deg_vec
Peagg_10KmDispDist_output_tab[,"habAv_imp"] = Peagg_10KmDispDist_habAv_vec
Peagg_10KmDispDist_output_tab[,"unw_b_c_imp"] = Peagg_10KmDispDist_unw_b_c_vec
Peagg_10KmDispDist_output_tab[,"Patch_Area_imp"] = Peagg_10KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Peagg_10KmDispDist_HSI_tab = data.frame(imp = Peagg_10KmDispDist_HSI_vec)
Peagg_10KmDispDist_EgoSize_tab = data.frame(imp = Peagg_10KmDispDist_EgoSize_vec)
Peagg_10KmDispDist_strength_tab = data.frame(imp = Peagg_10KmDispDist_strength_vec)
Peagg_10KmDispDist_deg_tab = data.frame(imp = Peagg_10KmDispDist_deg_vec)
Peagg_10KmDispDist_habAv_tab = data.frame(imp = Peagg_10KmDispDist_habAv_vec)
Peagg_10KmDispDist_unw_b_c_tab = data.frame(imp = Peagg_10KmDispDist_unw_b_c_vec)
Peagg_10KmDispDist_Patch_Area_tab = data.frame(imp = Peagg_10KmDispDist_Patch_Area_vec)

Peagg_10KmDispDist_HSI_tab$variable = "HSI"
Peagg_10KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Peagg_10KmDispDist_strength_tab$variable = "Strength"
Peagg_10KmDispDist_deg_tab$variable = "Degree"
Peagg_10KmDispDist_habAv_tab$variable = "Hab. Av."
Peagg_10KmDispDist_unw_b_c_tab$variable = "B.C."
Peagg_10KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Peagg_10KmDispDist_AUC_cv_tab = data.frame(value = Peagg_10KmDispDist_AUC_cv_vec)
Peagg_10KmDispDist_AUC_train_tab = data.frame(value = Peagg_10KmDispDist_AUC_vec)
#Make label of model for plot
Peagg_10KmDispDist_AUC_cv_tab$model = "Peagg_10KmDispDist"
Peagg_10KmDispDist_AUC_train_tab$model = "Peagg_10KmDispDist"

#combine pred. vars. into new data frame
Peagg_10KmDispDist_var_imp_tab = rbind(Peagg_10KmDispDist_HSI_tab,Peagg_10KmDispDist_EgoSize_tab,Peagg_10KmDispDist_strength_tab,Peagg_10KmDispDist_habAv_tab,Peagg_10KmDispDist_deg_tab,Peagg_10KmDispDist_unw_b_c_tab,Peagg_10KmDispDist_Patch_Area_tab)

ggplot(Peagg_10KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Peagg_10KmDispDist_var_imp_tab$imp~Peagg_10KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "10 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Peagg_10KmDispDist_HSI_vec)
mean(Peagg_10KmDispDist_EgoSize_vec)
mean(Peagg_10KmDispDist_strength_vec)
mean(Peagg_10KmDispDist_deg_vec)
mean(Peagg_10KmDispDist_habAv_vec)
mean(Peagg_10KmDispDist_unw_b_c_vec)
mean(Peagg_10KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Peagg_10KmDispDist_output_tab$AUC_cv)
summary(Peagg_10KmDispDist_output_tab$AUC_train)


####################################### 1Km ####################################

# Create the output table that contains all the values of each of the runs
Peagg_1KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg

# Create vectors to record performance measures
Peagg_1KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Peagg_1KmDispDist_AUC_vec = vector() #

Peagg_1KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Peagg_1KmDispDist_HSI_vec = vector() 
Peagg_1KmDispDist_EgoSize_vec = vector() 
Peagg_1KmDispDist_strength_vec  = vector()
Peagg_1KmDispDist_deg_vec  = vector()
Peagg_1KmDispDist_habAv_vec  = vector()
Peagg_1KmDispDist_unw_b_c_vec  = vector()
Peagg_1KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Peagg_1KmDispDist_predict_mat = matrix(nrow = length(Peagg_1KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Peagg_1KmDispDist = gbm.step(data=Peagg_1KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Peagg_1KmDispDist)$var, imp = summary(gbm_mod_Peagg_1KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Peagg_1KmDispDist_HSI_vec = append(Peagg_1KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Peagg_1KmDispDist_EgoSize_vec = append(Peagg_1KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Peagg_1KmDispDist_strength_vec  = append(Peagg_1KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Peagg_1KmDispDist_deg_vec  = append(Peagg_1KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Peagg_1KmDispDist_habAv_vec  = append(Peagg_1KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Peagg_1KmDispDist_unw_b_c_vec  = append(Peagg_1KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Peagg_1KmDispDist_Patch_Area_vec  = append(Peagg_1KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Peagg_1KmDispDist_AUC_cv_vec = append(Peagg_1KmDispDist_AUC_cv_vec, gbm_mod_Peagg_1KmDispDist$cv.statistics$discrimination.mean)
  Peagg_1KmDispDist_AUC_vec = append(Peagg_1KmDispDist_AUC_vec, gbm_mod_Peagg_1KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Peagg_1KmDispDist$n.trees
  Peagg_1KmDispDist_nt_vec = append(Peagg_1KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Peagg_1KmDispDist_predict_mat[,i] = predict(gbm_mod_Peagg_1KmDispDist, Peagg_1KmDispDist_stattest, gbm_mod_Peagg_1KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Peagg_1KmDispDist_output_tab[,"AUC_train"] = Peagg_1KmDispDist_AUC_vec
Peagg_1KmDispDist_output_tab[,"AUC_cv"] = Peagg_1KmDispDist_AUC_cv_vec
Peagg_1KmDispDist_output_tab[,"ntrees"] = Peagg_1KmDispDist_nt_vec

#Var. importance columns
Peagg_1KmDispDist_output_tab[,"HSI_imp"] = Peagg_1KmDispDist_HSI_vec
Peagg_1KmDispDist_output_tab[,"EgoSize_imp"] = Peagg_1KmDispDist_EgoSize_vec
Peagg_1KmDispDist_output_tab[,"strength_imp"] = Peagg_1KmDispDist_strength_vec
Peagg_1KmDispDist_output_tab[,"deg_imp"] = Peagg_1KmDispDist_deg_vec
Peagg_1KmDispDist_output_tab[,"habAv_imp"] = Peagg_1KmDispDist_habAv_vec
Peagg_1KmDispDist_output_tab[,"unw_b_c_imp"] = Peagg_1KmDispDist_unw_b_c_vec
Peagg_1KmDispDist_output_tab[,"Patch_Area_imp"] = Peagg_1KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Peagg_1KmDispDist_HSI_tab = data.frame(imp = Peagg_1KmDispDist_HSI_vec)
Peagg_1KmDispDist_EgoSize_tab = data.frame(imp = Peagg_1KmDispDist_EgoSize_vec)
Peagg_1KmDispDist_strength_tab = data.frame(imp = Peagg_1KmDispDist_strength_vec)
Peagg_1KmDispDist_deg_tab = data.frame(imp = Peagg_1KmDispDist_deg_vec)
Peagg_1KmDispDist_habAv_tab = data.frame(imp = Peagg_1KmDispDist_habAv_vec)
Peagg_1KmDispDist_unw_b_c_tab = data.frame(imp = Peagg_1KmDispDist_unw_b_c_vec)
Peagg_1KmDispDist_Patch_Area_tab = data.frame(imp = Peagg_1KmDispDist_Patch_Area_vec)

Peagg_1KmDispDist_HSI_tab$variable = "HSI"
Peagg_1KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Peagg_1KmDispDist_strength_tab$variable = "Strength"
Peagg_1KmDispDist_deg_tab$variable = "Degree"
Peagg_1KmDispDist_habAv_tab$variable = "Hab. Av."
Peagg_1KmDispDist_unw_b_c_tab$variable = "B.C."
Peagg_1KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Peagg_1KmDispDist_AUC_cv_tab = data.frame(value = Peagg_1KmDispDist_AUC_cv_vec)
Peagg_1KmDispDist_AUC_train_tab = data.frame(value = Peagg_1KmDispDist_AUC_vec)
#Make label of model for plot
Peagg_1KmDispDist_AUC_cv_tab$model = "Peagg_1KmDispDist"
Peagg_1KmDispDist_AUC_train_tab$model = "Peagg_1KmDispDist"

#combine pred. vars. into new data frame 
Peagg_1KmDispDist_var_imp_tab = rbind(Peagg_1KmDispDist_HSI_tab,Peagg_1KmDispDist_EgoSize_tab,Peagg_1KmDispDist_strength_tab,Peagg_1KmDispDist_habAv_tab,Peagg_1KmDispDist_deg_tab,Peagg_1KmDispDist_unw_b_c_tab,Peagg_1KmDispDist_Patch_Area_tab)

ggplot(Peagg_1KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Peagg_1KmDispDist_var_imp_tab$imp~Peagg_1KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "1 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Peagg_1KmDispDist_HSI_vec)
mean(Peagg_1KmDispDist_EgoSize_vec)
mean(Peagg_1KmDispDist_strength_vec)
mean(Peagg_1KmDispDist_deg_vec)
mean(Peagg_1KmDispDist_habAv_vec)
mean(Peagg_1KmDispDist_unw_b_c_vec)
mean(Peagg_1KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Peagg_1KmDispDist_output_tab$AUC_cv)
summary(Peagg_1KmDispDist_output_tab$AUC_train)


####################################### 6Km ####################################

# Create the output table that contains all the values of each of the runs
Peagg_6KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg

# Create vectors to record performance measures
Peagg_6KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Peagg_6KmDispDist_AUC_vec = vector() #

Peagg_6KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Peagg_6KmDispDist_HSI_vec = vector() 
Peagg_6KmDispDist_EgoSize_vec = vector() 
Peagg_6KmDispDist_strength_vec  = vector()
Peagg_6KmDispDist_deg_vec  = vector()
Peagg_6KmDispDist_habAv_vec  = vector()
Peagg_6KmDispDist_unw_b_c_vec  = vector()
Peagg_6KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Peagg_6KmDispDist_predict_mat = matrix(nrow = length(Peagg_6KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Peagg_6KmDispDist = gbm.step(data=Peagg_6KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Peagg_6KmDispDist)$var, imp = summary(gbm_mod_Peagg_6KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Peagg_6KmDispDist_HSI_vec = append(Peagg_6KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Peagg_6KmDispDist_EgoSize_vec = append(Peagg_6KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Peagg_6KmDispDist_strength_vec  = append(Peagg_6KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Peagg_6KmDispDist_deg_vec  = append(Peagg_6KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Peagg_6KmDispDist_habAv_vec  = append(Peagg_6KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Peagg_6KmDispDist_unw_b_c_vec  = append(Peagg_6KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Peagg_6KmDispDist_Patch_Area_vec  = append(Peagg_6KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Peagg_6KmDispDist_AUC_cv_vec = append(Peagg_6KmDispDist_AUC_cv_vec, gbm_mod_Peagg_6KmDispDist$cv.statistics$discrimination.mean)
  Peagg_6KmDispDist_AUC_vec = append(Peagg_6KmDispDist_AUC_vec, gbm_mod_Peagg_6KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Peagg_6KmDispDist$n.trees
  Peagg_6KmDispDist_nt_vec = append(Peagg_6KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Peagg_6KmDispDist_predict_mat[,i] = predict(gbm_mod_Peagg_6KmDispDist, Peagg_6KmDispDist_stattest, gbm_mod_Peagg_6KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Peagg_6KmDispDist_output_tab[,"AUC_train"] = Peagg_6KmDispDist_AUC_vec
Peagg_6KmDispDist_output_tab[,"AUC_cv"] = Peagg_6KmDispDist_AUC_cv_vec
Peagg_6KmDispDist_output_tab[,"ntrees"] = Peagg_6KmDispDist_nt_vec

#Var. importance columns
Peagg_6KmDispDist_output_tab[,"HSI_imp"] = Peagg_6KmDispDist_HSI_vec
Peagg_6KmDispDist_output_tab[,"EgoSize_imp"] = Peagg_6KmDispDist_EgoSize_vec
Peagg_6KmDispDist_output_tab[,"strength_imp"] = Peagg_6KmDispDist_strength_vec
Peagg_6KmDispDist_output_tab[,"deg_imp"] = Peagg_6KmDispDist_deg_vec
Peagg_6KmDispDist_output_tab[,"habAv_imp"] = Peagg_6KmDispDist_habAv_vec
Peagg_6KmDispDist_output_tab[,"unw_b_c_imp"] = Peagg_6KmDispDist_unw_b_c_vec
Peagg_6KmDispDist_output_tab[,"Patch_Area_imp"] = Peagg_6KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Peagg_6KmDispDist_HSI_tab = data.frame(imp = Peagg_6KmDispDist_HSI_vec)
Peagg_6KmDispDist_EgoSize_tab = data.frame(imp = Peagg_6KmDispDist_EgoSize_vec)
Peagg_6KmDispDist_strength_tab = data.frame(imp = Peagg_6KmDispDist_strength_vec)
Peagg_6KmDispDist_deg_tab = data.frame(imp = Peagg_6KmDispDist_deg_vec)
Peagg_6KmDispDist_habAv_tab = data.frame(imp = Peagg_6KmDispDist_habAv_vec)
Peagg_6KmDispDist_unw_b_c_tab = data.frame(imp = Peagg_6KmDispDist_unw_b_c_vec)
Peagg_6KmDispDist_Patch_Area_tab = data.frame(imp = Peagg_6KmDispDist_Patch_Area_vec)

Peagg_6KmDispDist_HSI_tab$variable = "HSI"
Peagg_6KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Peagg_6KmDispDist_strength_tab$variable = "Strength"
Peagg_6KmDispDist_deg_tab$variable = "Degree"
Peagg_6KmDispDist_habAv_tab$variable = "Hab. Av."
Peagg_6KmDispDist_unw_b_c_tab$variable = "B.C."
Peagg_6KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Peagg_6KmDispDist_AUC_cv_tab = data.frame(value = Peagg_6KmDispDist_AUC_cv_vec)
Peagg_6KmDispDist_AUC_train_tab = data.frame(value = Peagg_6KmDispDist_AUC_vec)
#Make label of model for plot
Peagg_6KmDispDist_AUC_cv_tab$model = "Peagg_6KmDispDist"
Peagg_6KmDispDist_AUC_train_tab$model = "Peagg_6KmDispDist"

#combine pred. vars. into new data frame 
Peagg_6KmDispDist_var_imp_tab = rbind(Peagg_6KmDispDist_HSI_tab,Peagg_6KmDispDist_EgoSize_tab,Peagg_6KmDispDist_strength_tab,Peagg_6KmDispDist_habAv_tab,Peagg_6KmDispDist_deg_tab,Peagg_6KmDispDist_unw_b_c_tab,Peagg_6KmDispDist_Patch_Area_tab)

ggplot(Peagg_6KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Peagg_6KmDispDist_var_imp_tab$imp~Peagg_6KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "6 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Peagg_6KmDispDist_HSI_vec)
mean(Peagg_6KmDispDist_EgoSize_vec)
mean(Peagg_6KmDispDist_strength_vec)
mean(Peagg_6KmDispDist_deg_vec)
mean(Peagg_6KmDispDist_habAv_vec)
mean(Peagg_6KmDispDist_unw_b_c_vec)
mean(Peagg_6KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Peagg_6KmDispDist_output_tab$AUC_cv)
summary(Peagg_6KmDispDist_output_tab$AUC_train)


####################################### 4Km ####################################

# Create the output table that contains all the values of each of the runs
Peagg_4KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg

# Create vectors to record performance measures
Peagg_4KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Peagg_4KmDispDist_AUC_vec = vector() #

Peagg_4KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Peagg_4KmDispDist_HSI_vec = vector()
Peagg_4KmDispDist_EgoSize_vec = vector()
Peagg_4KmDispDist_strength_vec  = vector()
Peagg_4KmDispDist_deg_vec  = vector()
Peagg_4KmDispDist_habAv_vec  = vector()
Peagg_4KmDispDist_unw_b_c_vec  = vector()
Peagg_4KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Peagg_4KmDispDist_predict_mat = matrix(nrow = length(Peagg_4KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Peagg_4KmDispDist = gbm.step(data=Peagg_4KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Peagg_4KmDispDist)$var, imp = summary(gbm_mod_Peagg_4KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Peagg_4KmDispDist_HSI_vec = append(Peagg_4KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Peagg_4KmDispDist_EgoSize_vec = append(Peagg_4KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Peagg_4KmDispDist_strength_vec  = append(Peagg_4KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Peagg_4KmDispDist_deg_vec  = append(Peagg_4KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Peagg_4KmDispDist_habAv_vec  = append(Peagg_4KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Peagg_4KmDispDist_unw_b_c_vec  = append(Peagg_4KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Peagg_4KmDispDist_Patch_Area_vec  = append(Peagg_4KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Peagg_4KmDispDist_AUC_cv_vec = append(Peagg_4KmDispDist_AUC_cv_vec, gbm_mod_Peagg_4KmDispDist$cv.statistics$discrimination.mean)
  Peagg_4KmDispDist_AUC_vec = append(Peagg_4KmDispDist_AUC_vec, gbm_mod_Peagg_4KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Peagg_4KmDispDist$n.trees
  Peagg_4KmDispDist_nt_vec = append(Peagg_4KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Peagg_4KmDispDist_predict_mat[,i] = predict(gbm_mod_Peagg_4KmDispDist, Peagg_4KmDispDist_stattest, gbm_mod_Peagg_4KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Peagg_4KmDispDist_output_tab[,"AUC_train"] = Peagg_4KmDispDist_AUC_vec
Peagg_4KmDispDist_output_tab[,"AUC_cv"] = Peagg_4KmDispDist_AUC_cv_vec
Peagg_4KmDispDist_output_tab[,"ntrees"] = Peagg_4KmDispDist_nt_vec

#Var. importance columns
Peagg_4KmDispDist_output_tab[,"HSI_imp"] = Peagg_4KmDispDist_HSI_vec
Peagg_4KmDispDist_output_tab[,"EgoSize_imp"] = Peagg_4KmDispDist_EgoSize_vec
Peagg_4KmDispDist_output_tab[,"strength_imp"] = Peagg_4KmDispDist_strength_vec
Peagg_4KmDispDist_output_tab[,"deg_imp"] = Peagg_4KmDispDist_deg_vec
Peagg_4KmDispDist_output_tab[,"habAv_imp"] = Peagg_4KmDispDist_habAv_vec
Peagg_4KmDispDist_output_tab[,"unw_b_c_imp"] = Peagg_4KmDispDist_unw_b_c_vec
Peagg_4KmDispDist_output_tab[,"Patch_Area_imp"] = Peagg_4KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Peagg_4KmDispDist_HSI_tab = data.frame(imp = Peagg_4KmDispDist_HSI_vec)
Peagg_4KmDispDist_EgoSize_tab = data.frame(imp = Peagg_4KmDispDist_EgoSize_vec)
Peagg_4KmDispDist_strength_tab = data.frame(imp = Peagg_4KmDispDist_strength_vec)
Peagg_4KmDispDist_deg_tab = data.frame(imp = Peagg_4KmDispDist_deg_vec)
Peagg_4KmDispDist_habAv_tab = data.frame(imp = Peagg_4KmDispDist_habAv_vec)
Peagg_4KmDispDist_unw_b_c_tab = data.frame(imp = Peagg_4KmDispDist_unw_b_c_vec)
Peagg_4KmDispDist_Patch_Area_tab = data.frame(imp = Peagg_4KmDispDist_Patch_Area_vec)

Peagg_4KmDispDist_HSI_tab$variable = "HSI"
Peagg_4KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Peagg_4KmDispDist_strength_tab$variable = "Strength"
Peagg_4KmDispDist_deg_tab$variable = "Degree"
Peagg_4KmDispDist_habAv_tab$variable = "Hab. Av."
Peagg_4KmDispDist_unw_b_c_tab$variable = "B.C."
Peagg_4KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Peagg_4KmDispDist_AUC_cv_tab = data.frame(value = Peagg_4KmDispDist_AUC_cv_vec)
Peagg_4KmDispDist_AUC_train_tab = data.frame(value = Peagg_4KmDispDist_AUC_vec)
#Make label of model for plot
Peagg_4KmDispDist_AUC_cv_tab$model = "Peagg_4KmDispDist"
Peagg_4KmDispDist_AUC_train_tab$model = "Peagg_4KmDispDist"

#combine pred. vars. into new data frame
Peagg_4KmDispDist_var_imp_tab = rbind(Peagg_4KmDispDist_HSI_tab,Peagg_4KmDispDist_EgoSize_tab,Peagg_4KmDispDist_strength_tab,Peagg_4KmDispDist_habAv_tab,Peagg_4KmDispDist_deg_tab,Peagg_4KmDispDist_unw_b_c_tab,Peagg_4KmDispDist_Patch_Area_tab)

ggplot(Peagg_4KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Peagg_4KmDispDist_var_imp_tab$imp~Peagg_4KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "4 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Peagg_4KmDispDist_HSI_vec)
mean(Peagg_4KmDispDist_EgoSize_vec)
mean(Peagg_4KmDispDist_strength_vec)
mean(Peagg_4KmDispDist_deg_vec)
mean(Peagg_4KmDispDist_habAv_vec)
mean(Peagg_4KmDispDist_unw_b_c_vec)
mean(Peagg_4KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Peagg_4KmDispDist_output_tab$AUC_cv)
summary(Peagg_4KmDispDist_output_tab$AUC_train)


####################################### 8Km ####################################

# Create the output table that contains all the values of each of the runs
Peagg_8KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg

# Create vectors to record performance measures
Peagg_8KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Peagg_8KmDispDist_AUC_vec = vector() #

Peagg_8KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Peagg_8KmDispDist_HSI_vec = vector() 
Peagg_8KmDispDist_EgoSize_vec = vector() 
Peagg_8KmDispDist_strength_vec  = vector()
Peagg_8KmDispDist_deg_vec  = vector()
Peagg_8KmDispDist_habAv_vec  = vector()
Peagg_8KmDispDist_unw_b_c_vec  = vector()
Peagg_8KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Peagg_8KmDispDist_predict_mat = matrix(nrow = length(Peagg_8KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Peagg_8KmDispDist = gbm.step(data=Peagg_8KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Peagg_8KmDispDist)$var, imp = summary(gbm_mod_Peagg_8KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Peagg_8KmDispDist_HSI_vec = append(Peagg_8KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Peagg_8KmDispDist_EgoSize_vec = append(Peagg_8KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Peagg_8KmDispDist_strength_vec  = append(Peagg_8KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Peagg_8KmDispDist_deg_vec  = append(Peagg_8KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Peagg_8KmDispDist_habAv_vec  = append(Peagg_8KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Peagg_8KmDispDist_unw_b_c_vec  = append(Peagg_8KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Peagg_8KmDispDist_Patch_Area_vec  = append(Peagg_8KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Peagg_8KmDispDist_AUC_cv_vec = append(Peagg_8KmDispDist_AUC_cv_vec, gbm_mod_Peagg_8KmDispDist$cv.statistics$discrimination.mean)
  Peagg_8KmDispDist_AUC_vec = append(Peagg_8KmDispDist_AUC_vec, gbm_mod_Peagg_8KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Peagg_8KmDispDist$n.trees
  Peagg_8KmDispDist_nt_vec = append(Peagg_8KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Peagg_8KmDispDist_predict_mat[,i] = predict(gbm_mod_Peagg_8KmDispDist, Peagg_8KmDispDist_stattest, gbm_mod_Peagg_8KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Peagg_8KmDispDist_output_tab[,"AUC_train"] = Peagg_8KmDispDist_AUC_vec
Peagg_8KmDispDist_output_tab[,"AUC_cv"] = Peagg_8KmDispDist_AUC_cv_vec
Peagg_8KmDispDist_output_tab[,"ntrees"] = Peagg_8KmDispDist_nt_vec

#Var. importance columns
Peagg_8KmDispDist_output_tab[,"HSI_imp"] = Peagg_8KmDispDist_HSI_vec
Peagg_8KmDispDist_output_tab[,"EgoSize_imp"] = Peagg_8KmDispDist_EgoSize_vec
Peagg_8KmDispDist_output_tab[,"strength_imp"] = Peagg_8KmDispDist_strength_vec
Peagg_8KmDispDist_output_tab[,"deg_imp"] = Peagg_8KmDispDist_deg_vec
Peagg_8KmDispDist_output_tab[,"habAv_imp"] = Peagg_8KmDispDist_habAv_vec
Peagg_8KmDispDist_output_tab[,"unw_b_c_imp"] = Peagg_8KmDispDist_unw_b_c_vec
Peagg_8KmDispDist_output_tab[,"Patch_Area_imp"] = Peagg_8KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Peagg_8KmDispDist_HSI_tab = data.frame(imp = Peagg_8KmDispDist_HSI_vec)
Peagg_8KmDispDist_EgoSize_tab = data.frame(imp = Peagg_8KmDispDist_EgoSize_vec)
Peagg_8KmDispDist_strength_tab = data.frame(imp = Peagg_8KmDispDist_strength_vec)
Peagg_8KmDispDist_deg_tab = data.frame(imp = Peagg_8KmDispDist_deg_vec)
Peagg_8KmDispDist_habAv_tab = data.frame(imp = Peagg_8KmDispDist_habAv_vec)
Peagg_8KmDispDist_unw_b_c_tab = data.frame(imp = Peagg_8KmDispDist_unw_b_c_vec)
Peagg_8KmDispDist_Patch_Area_tab = data.frame(imp = Peagg_8KmDispDist_Patch_Area_vec)

Peagg_8KmDispDist_HSI_tab$variable = "HSI"
Peagg_8KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Peagg_8KmDispDist_strength_tab$variable = "Strength"
Peagg_8KmDispDist_deg_tab$variable = "Degree"
Peagg_8KmDispDist_habAv_tab$variable = "Hab. Av."
Peagg_8KmDispDist_unw_b_c_tab$variable = "B.C."
Peagg_8KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Peagg_8KmDispDist_AUC_cv_tab = data.frame(value = Peagg_8KmDispDist_AUC_cv_vec)
Peagg_8KmDispDist_AUC_train_tab = data.frame(value = Peagg_8KmDispDist_AUC_vec)
#Make label of model for plot
Peagg_8KmDispDist_AUC_cv_tab$model = "Peagg_8KmDispDist"
Peagg_8KmDispDist_AUC_train_tab$model = "Peagg_8KmDispDist"

#combine pred. vars. into new data frame 
Peagg_8KmDispDist_var_imp_tab = rbind(Peagg_8KmDispDist_HSI_tab,Peagg_8KmDispDist_EgoSize_tab,Peagg_8KmDispDist_strength_tab,Peagg_8KmDispDist_habAv_tab,Peagg_8KmDispDist_deg_tab,Peagg_8KmDispDist_unw_b_c_tab,Peagg_8KmDispDist_Patch_Area_tab)

ggplot(Peagg_8KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Peagg_8KmDispDist_var_imp_tab$imp~Peagg_8KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "8 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Peagg_8KmDispDist_HSI_vec)
mean(Peagg_8KmDispDist_EgoSize_vec)
mean(Peagg_8KmDispDist_strength_vec)
mean(Peagg_8KmDispDist_deg_vec)
mean(Peagg_8KmDispDist_habAv_vec)
mean(Peagg_8KmDispDist_unw_b_c_vec)
mean(Peagg_8KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Peagg_8KmDispDist_output_tab$AUC_cv)
summary(Peagg_8KmDispDist_output_tab$AUC_train)


####################################################################################################
#### Compare scores between networks w/d0 variations of the same species ###########################
### Peagg ##################

#Import evaluation df's of original run with species-specific dispersal distance
Peagg_DefaultDispDist_AUC_cv_tab <- read.csv("Peagg_DefaultDispDist_AUC_cv_tab.csv")
Peagg_DefaultDispDist_AUC_train_tab <- read.csv("Peagg_DefaultDispDist_AUC_train_tab.csv")
Peagg_DefaultDispDist_var_imp_tab <- read.csv("Peagg_DefaultDispDist_var_imp_tab.csv")
Peagg_DefaultDispDist_output_tab <- read.csv("Peagg_DefaultDispDist_BRToutput_tab.csv")
Peagg_DefaultDispDist_predict_df <- read.csv("Peagg_DefaultDispDist_predict_df.csv")

head(Peagg_DefaultDispDist_AUC_cv_tab)
Peagg_DefaultDispDist_AUC_cv_tab$X <- NULL
Peagg_DefaultDispDist_AUC_train_tab$X <- NULL

#### Make a dataframe for plotting overlaying histograms in R
### cv AUC
Peagg_AUC_cv_tab = rbind(Peagg_DefaultDispDist_AUC_cv_tab, Peagg_300mDispDist_AUC_cv_tab, Peagg_1KmDispDist_AUC_cv_tab, Peagg_2KmDispDist_AUC_cv_tab, 
                         Peagg_4KmDispDist_AUC_cv_tab,
                         Peagg_6KmDispDist_AUC_cv_tab, Peagg_8KmDispDist_AUC_cv_tab, Peagg_10KmDispDist_AUC_cv_tab)

#Change order of factors to display noTopo at the edge
Peagg_AUC_cv_tab$model<- as.factor(Peagg_AUC_cv_tab$model)
levels(Peagg_AUC_cv_tab$model) #The level for DefaultDispDist is called only with the species name
Peagg_AUC_cv_tab$model<-factor(Peagg_AUC_cv_tab$model, levels=c("Peagg_300mDispDist", "Peagg_1KmDispDist", "Peagg", "Peagg_2KmDispDist",
                                                                "Peagg_4KmDispDist", 
                                                                "Peagg_6KmDispDist", "Peagg_8KmDispDist", "Peagg_10KmDispDist"))
#Plot
ggplot(Peagg_AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Peagg_AUC_cv_tab$value~Peagg_AUC_cv_tab$model, ylab= "Cross-validated AUC")

### train AUC
Peagg_AUC_train_tab = rbind(Peagg_DefaultDispDist_AUC_train_tab, Peagg_300mDispDist_AUC_train_tab, Peagg_1KmDispDist_AUC_train_tab, Peagg_2KmDispDist_AUC_train_tab, 
                            Peagg_4KmDispDist_AUC_train_tab,
                            Peagg_6KmDispDist_AUC_train_tab, Peagg_8KmDispDist_AUC_train_tab, Peagg_10KmDispDist_AUC_train_tab)

levels(Peagg_AUC_train_tab$model)
Peagg_AUC_train_tab$model<-factor(Peagg_AUC_train_tab$model, levels=c("Peagg_300mDispDist", "Peagg_1KmDispDist", "Peagg", "Peagg_2KmDispDist",
                                                                      "Peagg_4KmDispDist", 
                                                                      "Peagg_6KmDispDist", "Peagg_8KmDispDist", "Peagg_10KmDispDist"))
#Plot
ggplot(Peagg_AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Peagg_AUC_train_tab$value~Peagg_AUC_train_tab$model, ylab= "Training AUC")

summary(Peagg_DefaultDispDist_output_tab$AUC_cv)
summary(Peagg_300mDispDist_output_tab$AUC_cv)
summary(Peagg_1KmDispDist_output_tab$AUC_cv)
summary(Peagg_2KmDispDist_output_tab$AUC_cv)
summary(Peagg_4KmDispDist_output_tab$AUC_cv)
summary(Peagg_6KmDispDist_output_tab$AUC_cv)
summary(Peagg_8KmDispDist_output_tab$AUC_cv)
summary(Peagg_10KmDispDist_output_tab$AUC_cv)

summary(Peagg_DefaultDispDist_output_tab$AUC_train)
summary(Peagg_300mDispDist_output_tab$AUC_train)
summary(Peagg_1KmDispDist_output_tab$AUC_train)
summary(Peagg_2KmDispDist_output_tab$AUC_train)
summary(Peagg_4KmDispDist_output_tab$AUC_train)
summary(Peagg_6KmDispDist_output_tab$AUC_train)
summary(Peagg_8KmDispDist_output_tab$AUC_train)
summary(Peagg_10KmDispDist_output_tab$AUC_train)

summary(Peagg_DefaultDispDist_output_tab$ntrees)
summary(Peagg_300mDispDist_output_tab$ntrees)
summary(Peagg_1KmDispDist_output_tab$ntrees)
summary(Peagg_2KmDispDist_output_tab$ntrees)
summary(Peagg_4KmDispDist_output_tab$ntrees)
summary(Peagg_6KmDispDist_output_tab$ntrees)
summary(Peagg_8KmDispDist_output_tab$ntrees)
summary(Peagg_10KmDispDist_output_tab$ntrees)



################################################################################
################ Perid #########################################################

#### 2Km #######################################################################

# Create the output table that contains all the values of each of the runs
Perid_2KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid

# Create vectors to record performance measures
Perid_2KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Perid_2KmDispDist_AUC_vec = vector() #

Perid_2KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Perid_2KmDispDist_HSI_vec = vector()
Perid_2KmDispDist_EgoSize_vec = vector()
Perid_2KmDispDist_strength_vec  = vector()
Perid_2KmDispDist_deg_vec  = vector()
Perid_2KmDispDist_habAv_vec  = vector()
Perid_2KmDispDist_unw_b_c_vec  = vector()
Perid_2KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Perid_2KmDispDist_predict_mat = matrix(nrow = length(Perid_2KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Perid_2KmDispDist = gbm.step(data=Perid_2KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Perid_2KmDispDist)$var, imp = summary(gbm_mod_Perid_2KmDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Perid_2KmDispDist_HSI_vec = append(Perid_2KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Perid_2KmDispDist_EgoSize_vec = append(Perid_2KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Perid_2KmDispDist_strength_vec  = append(Perid_2KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Perid_2KmDispDist_deg_vec  = append(Perid_2KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Perid_2KmDispDist_habAv_vec  = append(Perid_2KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Perid_2KmDispDist_unw_b_c_vec  = append(Perid_2KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Perid_2KmDispDist_Patch_Area_vec  = append(Perid_2KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Perid_2KmDispDist_AUC_cv_vec = append(Perid_2KmDispDist_AUC_cv_vec, gbm_mod_Perid_2KmDispDist$cv.statistics$discrimination.mean)
  Perid_2KmDispDist_AUC_vec = append(Perid_2KmDispDist_AUC_vec, gbm_mod_Perid_2KmDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Perid_2KmDispDist$n.trees
  Perid_2KmDispDist_nt_vec = append(Perid_2KmDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Perid_2KmDispDist_predict_mat[,i] = predict(gbm_mod_Perid_2KmDispDist, Perid_2KmDispDist_stattest, gbm_mod_Perid_2KmDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Perid_2KmDispDist_output_tab[,"AUC_train"] = Perid_2KmDispDist_AUC_vec
Perid_2KmDispDist_output_tab[,"AUC_cv"] = Perid_2KmDispDist_AUC_cv_vec
Perid_2KmDispDist_output_tab[,"ntrees"] = Perid_2KmDispDist_nt_vec

#Var. importance columns
Perid_2KmDispDist_output_tab[,"HSI_imp"] = Perid_2KmDispDist_HSI_vec
Perid_2KmDispDist_output_tab[,"EgoSize_imp"] = Perid_2KmDispDist_EgoSize_vec
Perid_2KmDispDist_output_tab[,"strength_imp"] = Perid_2KmDispDist_strength_vec
Perid_2KmDispDist_output_tab[,"deg_imp"] = Perid_2KmDispDist_deg_vec
Perid_2KmDispDist_output_tab[,"habAv_imp"] = Perid_2KmDispDist_habAv_vec
Perid_2KmDispDist_output_tab[,"unw_b_c_imp"] = Perid_2KmDispDist_unw_b_c_vec
Perid_2KmDispDist_output_tab[,"Patch_Area_imp"] = Perid_2KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Perid_2KmDispDist_HSI_tab = data.frame(imp = Perid_2KmDispDist_HSI_vec)
Perid_2KmDispDist_EgoSize_tab = data.frame(imp = Perid_2KmDispDist_EgoSize_vec)
Perid_2KmDispDist_strength_tab = data.frame(imp = Perid_2KmDispDist_strength_vec)
Perid_2KmDispDist_deg_tab = data.frame(imp = Perid_2KmDispDist_deg_vec)
Perid_2KmDispDist_habAv_tab = data.frame(imp = Perid_2KmDispDist_habAv_vec)
Perid_2KmDispDist_unw_b_c_tab = data.frame(imp = Perid_2KmDispDist_unw_b_c_vec)
Perid_2KmDispDist_Patch_Area_tab = data.frame(imp = Perid_2KmDispDist_Patch_Area_vec)

Perid_2KmDispDist_HSI_tab$variable = "HSI"
Perid_2KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Perid_2KmDispDist_strength_tab$variable = "Strength"
Perid_2KmDispDist_deg_tab$variable = "Degree"
Perid_2KmDispDist_habAv_tab$variable = "Hab. Av."
Perid_2KmDispDist_unw_b_c_tab$variable = "B.C."
Perid_2KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Perid_2KmDispDist_AUC_cv_tab = data.frame(value = Perid_2KmDispDist_AUC_cv_vec)
Perid_2KmDispDist_AUC_train_tab = data.frame(value = Perid_2KmDispDist_AUC_vec)
#Make label of model for plot
Perid_2KmDispDist_AUC_cv_tab$model = "Perid_2KmDispDist"
Perid_2KmDispDist_AUC_train_tab$model = "Perid_2KmDispDist"

#combine pred. vars. into new data frame
Perid_2KmDispDist_var_imp_tab = rbind(Perid_2KmDispDist_HSI_tab,Perid_2KmDispDist_EgoSize_tab,Perid_2KmDispDist_strength_tab,Perid_2KmDispDist_habAv_tab,Perid_2KmDispDist_deg_tab,Perid_2KmDispDist_unw_b_c_tab,Perid_2KmDispDist_Patch_Area_tab)

ggplot(Perid_2KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Perid_2KmDispDist_var_imp_tab$imp~Perid_2KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "2 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Perid_2KmDispDist_HSI_vec)
mean(Perid_2KmDispDist_EgoSize_vec)
mean(Perid_2KmDispDist_strength_vec)
mean(Perid_2KmDispDist_deg_vec)
mean(Perid_2KmDispDist_habAv_vec)
mean(Perid_2KmDispDist_unw_b_c_vec)
mean(Perid_2KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Perid_2KmDispDist_output_tab$AUC_cv)
summary(Perid_2KmDispDist_output_tab$AUC_train)


# #### 300m ####################################################################

# Create the output table that contains all the values of each of the runs
Perid_300mDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid

# Create vectors to record performance measures
Perid_300mDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Perid_300mDispDist_AUC_vec = vector() #

Perid_300mDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Perid_300mDispDist_HSI_vec = vector()
Perid_300mDispDist_EgoSize_vec = vector()
Perid_300mDispDist_strength_vec  = vector()
Perid_300mDispDist_deg_vec  = vector()
Perid_300mDispDist_habAv_vec  = vector()
Perid_300mDispDist_unw_b_c_vec  = vector()
Perid_300mDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Perid_300mDispDist_predict_mat = matrix(nrow = length(Perid_300mDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Perid_300mDispDist = gbm.step(data=Perid_300mDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Perid_300mDispDist)$var, imp = summary(gbm_mod_Perid_300mDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Perid_300mDispDist_HSI_vec = append(Perid_300mDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Perid_300mDispDist_EgoSize_vec = append(Perid_300mDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Perid_300mDispDist_strength_vec  = append(Perid_300mDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Perid_300mDispDist_deg_vec  = append(Perid_300mDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Perid_300mDispDist_habAv_vec  = append(Perid_300mDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Perid_300mDispDist_unw_b_c_vec  = append(Perid_300mDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Perid_300mDispDist_Patch_Area_vec  = append(Perid_300mDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Perid_300mDispDist_AUC_cv_vec = append(Perid_300mDispDist_AUC_cv_vec, gbm_mod_Perid_300mDispDist$cv.statistics$discrimination.mean)
  Perid_300mDispDist_AUC_vec = append(Perid_300mDispDist_AUC_vec, gbm_mod_Perid_300mDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Perid_300mDispDist$n.trees
  Perid_300mDispDist_nt_vec = append(Perid_300mDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Perid_300mDispDist_predict_mat[,i] = predict(gbm_mod_Perid_300mDispDist, Perid_300mDispDist_stattest, gbm_mod_Perid_300mDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Perid_300mDispDist_output_tab[,"AUC_train"] = Perid_300mDispDist_AUC_vec
Perid_300mDispDist_output_tab[,"AUC_cv"] = Perid_300mDispDist_AUC_cv_vec
Perid_300mDispDist_output_tab[,"ntrees"] = Perid_300mDispDist_nt_vec

#Var. importance columns
Perid_300mDispDist_output_tab[,"HSI_imp"] = Perid_300mDispDist_HSI_vec
Perid_300mDispDist_output_tab[,"EgoSize_imp"] = Perid_300mDispDist_EgoSize_vec
Perid_300mDispDist_output_tab[,"strength_imp"] = Perid_300mDispDist_strength_vec
Perid_300mDispDist_output_tab[,"deg_imp"] = Perid_300mDispDist_deg_vec
Perid_300mDispDist_output_tab[,"habAv_imp"] = Perid_300mDispDist_habAv_vec
Perid_300mDispDist_output_tab[,"unw_b_c_imp"] = Perid_300mDispDist_unw_b_c_vec
Perid_300mDispDist_output_tab[,"Patch_Area_imp"] = Perid_300mDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Perid_300mDispDist_HSI_tab = data.frame(imp = Perid_300mDispDist_HSI_vec)
Perid_300mDispDist_EgoSize_tab = data.frame(imp = Perid_300mDispDist_EgoSize_vec)
Perid_300mDispDist_strength_tab = data.frame(imp = Perid_300mDispDist_strength_vec)
Perid_300mDispDist_deg_tab = data.frame(imp = Perid_300mDispDist_deg_vec)
Perid_300mDispDist_habAv_tab = data.frame(imp = Perid_300mDispDist_habAv_vec)
Perid_300mDispDist_unw_b_c_tab = data.frame(imp = Perid_300mDispDist_unw_b_c_vec)
Perid_300mDispDist_Patch_Area_tab = data.frame(imp = Perid_300mDispDist_Patch_Area_vec)

Perid_300mDispDist_HSI_tab$variable = "HSI"
Perid_300mDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Perid_300mDispDist_strength_tab$variable = "Strength"
Perid_300mDispDist_deg_tab$variable = "Degree"
Perid_300mDispDist_habAv_tab$variable = "Hab. Av."
Perid_300mDispDist_unw_b_c_tab$variable = "B.C."
Perid_300mDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Perid_300mDispDist_AUC_cv_tab = data.frame(value = Perid_300mDispDist_AUC_cv_vec)
Perid_300mDispDist_AUC_train_tab = data.frame(value = Perid_300mDispDist_AUC_vec)
#Make label of model for plot
Perid_300mDispDist_AUC_cv_tab$model = "Perid_300mDispDist"
Perid_300mDispDist_AUC_train_tab$model = "Perid_300mDispDist"

#combine pred. vars. into new data frame
Perid_300mDispDist_var_imp_tab = rbind(Perid_300mDispDist_HSI_tab,Perid_300mDispDist_EgoSize_tab,Perid_300mDispDist_strength_tab,Perid_300mDispDist_habAv_tab,Perid_300mDispDist_deg_tab,Perid_300mDispDist_unw_b_c_tab,Perid_300mDispDist_Patch_Area_tab)

ggplot(Perid_300mDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Perid_300mDispDist_var_imp_tab$imp~Perid_300mDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "300 m maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Perid_300mDispDist_HSI_vec)
mean(Perid_300mDispDist_EgoSize_vec)
mean(Perid_300mDispDist_strength_vec)
mean(Perid_300mDispDist_deg_vec)
mean(Perid_300mDispDist_habAv_vec)
mean(Perid_300mDispDist_unw_b_c_vec)
mean(Perid_300mDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Perid_300mDispDist_output_tab$AUC_cv)
summary(Perid_300mDispDist_output_tab$AUC_train)


# #### 10Km ####################################################################

# Create the output table that contains all the values of each of the runs
Perid_10KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid

# Create vectors to record performance measures
Perid_10KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Perid_10KmDispDist_AUC_vec = vector() #

Perid_10KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Perid_10KmDispDist_HSI_vec = vector()
Perid_10KmDispDist_EgoSize_vec = vector()
Perid_10KmDispDist_strength_vec  = vector()
Perid_10KmDispDist_deg_vec  = vector()
Perid_10KmDispDist_habAv_vec  = vector()
Perid_10KmDispDist_unw_b_c_vec  = vector()
Perid_10KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Perid_10KmDispDist_predict_mat = matrix(nrow = length(Perid_10KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){

  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Perid_10KmDispDist = gbm.step(data=Perid_10KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Perid_10KmDispDist)$var, imp = summary(gbm_mod_Perid_10KmDispDist)$rel.inf)

  #Write var_imp results to the vectors
  Perid_10KmDispDist_HSI_vec = append(Perid_10KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Perid_10KmDispDist_EgoSize_vec = append(Perid_10KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Perid_10KmDispDist_strength_vec  = append(Perid_10KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Perid_10KmDispDist_deg_vec  = append(Perid_10KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Perid_10KmDispDist_habAv_vec  = append(Perid_10KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Perid_10KmDispDist_unw_b_c_vec  = append(Perid_10KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Perid_10KmDispDist_Patch_Area_vec  = append(Perid_10KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])

  #Write the AUC & CV-AUC of this model to a vector.
  Perid_10KmDispDist_AUC_cv_vec = append(Perid_10KmDispDist_AUC_cv_vec, gbm_mod_Perid_10KmDispDist$cv.statistics$discrimination.mean)
  Perid_10KmDispDist_AUC_vec = append(Perid_10KmDispDist_AUC_vec, gbm_mod_Perid_10KmDispDist$self.statistics$discrimination)

  #Write the number of trees to a vector
  nt = gbm_mod_Perid_10KmDispDist$n.trees
  Perid_10KmDispDist_nt_vec = append(Perid_10KmDispDist_nt_vec, nt)
  print(nt)

  #write the continuous prediction over all the patches
  Perid_10KmDispDist_predict_mat[,i] = predict(gbm_mod_Perid_10KmDispDist, Perid_10KmDispDist_stattest, gbm_mod_Perid_10KmDispDist$n.trees, type = "response", single.tree = FALSE)

  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Perid_10KmDispDist_output_tab[,"AUC_train"] = Perid_10KmDispDist_AUC_vec
Perid_10KmDispDist_output_tab[,"AUC_cv"] = Perid_10KmDispDist_AUC_cv_vec
Perid_10KmDispDist_output_tab[,"ntrees"] = Perid_10KmDispDist_nt_vec

#Var. importance columns
Perid_10KmDispDist_output_tab[,"HSI_imp"] = Perid_10KmDispDist_HSI_vec
Perid_10KmDispDist_output_tab[,"EgoSize_imp"] = Perid_10KmDispDist_EgoSize_vec
Perid_10KmDispDist_output_tab[,"strength_imp"] = Perid_10KmDispDist_strength_vec
Perid_10KmDispDist_output_tab[,"deg_imp"] = Perid_10KmDispDist_deg_vec
Perid_10KmDispDist_output_tab[,"habAv_imp"] = Perid_10KmDispDist_habAv_vec
Perid_10KmDispDist_output_tab[,"unw_b_c_imp"] = Perid_10KmDispDist_unw_b_c_vec
Perid_10KmDispDist_output_tab[,"Patch_Area_imp"] = Perid_10KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Perid_10KmDispDist_HSI_tab = data.frame(imp = Perid_10KmDispDist_HSI_vec)
Perid_10KmDispDist_EgoSize_tab = data.frame(imp = Perid_10KmDispDist_EgoSize_vec)
Perid_10KmDispDist_strength_tab = data.frame(imp = Perid_10KmDispDist_strength_vec)
Perid_10KmDispDist_deg_tab = data.frame(imp = Perid_10KmDispDist_deg_vec)
Perid_10KmDispDist_habAv_tab = data.frame(imp = Perid_10KmDispDist_habAv_vec)
Perid_10KmDispDist_unw_b_c_tab = data.frame(imp = Perid_10KmDispDist_unw_b_c_vec)
Perid_10KmDispDist_Patch_Area_tab = data.frame(imp = Perid_10KmDispDist_Patch_Area_vec)

Perid_10KmDispDist_HSI_tab$variable = "HSI"
Perid_10KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Perid_10KmDispDist_strength_tab$variable = "Strength"
Perid_10KmDispDist_deg_tab$variable = "Degree"
Perid_10KmDispDist_habAv_tab$variable = "Hab. Av."
Perid_10KmDispDist_unw_b_c_tab$variable = "B.C."
Perid_10KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Perid_10KmDispDist_AUC_cv_tab = data.frame(value = Perid_10KmDispDist_AUC_cv_vec)
Perid_10KmDispDist_AUC_train_tab = data.frame(value = Perid_10KmDispDist_AUC_vec)
#Make label of model for plot
Perid_10KmDispDist_AUC_cv_tab$model = "Perid_10KmDispDist"
Perid_10KmDispDist_AUC_train_tab$model = "Perid_10KmDispDist"

#combine pred. vars. into new data frame
Perid_10KmDispDist_var_imp_tab = rbind(Perid_10KmDispDist_HSI_tab,Perid_10KmDispDist_EgoSize_tab,Perid_10KmDispDist_strength_tab,Perid_10KmDispDist_habAv_tab,Perid_10KmDispDist_deg_tab,Perid_10KmDispDist_unw_b_c_tab,Perid_10KmDispDist_Patch_Area_tab)

ggplot(Perid_10KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Perid_10KmDispDist_var_imp_tab$imp~Perid_10KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "10 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Perid_10KmDispDist_HSI_vec)
mean(Perid_10KmDispDist_EgoSize_vec)
mean(Perid_10KmDispDist_strength_vec)
mean(Perid_10KmDispDist_deg_vec)
mean(Perid_10KmDispDist_habAv_vec)
mean(Perid_10KmDispDist_unw_b_c_vec)
mean(Perid_10KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Perid_10KmDispDist_output_tab$AUC_cv)
summary(Perid_10KmDispDist_output_tab$AUC_train)


####################################### 1Km ####################################

# Create the output table that contains all the values of each of the runs
Perid_1KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid

# Create vectors to record performance measures
Perid_1KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Perid_1KmDispDist_AUC_vec = vector() #

Perid_1KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Perid_1KmDispDist_HSI_vec = vector() 
Perid_1KmDispDist_EgoSize_vec = vector() 
Perid_1KmDispDist_strength_vec  = vector()
Perid_1KmDispDist_deg_vec  = vector()
Perid_1KmDispDist_habAv_vec  = vector()
Perid_1KmDispDist_unw_b_c_vec  = vector()
Perid_1KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Perid_1KmDispDist_predict_mat = matrix(nrow = length(Perid_1KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Perid_1KmDispDist = gbm.step(data=Perid_1KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Perid_1KmDispDist)$var, imp = summary(gbm_mod_Perid_1KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Perid_1KmDispDist_HSI_vec = append(Perid_1KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Perid_1KmDispDist_EgoSize_vec = append(Perid_1KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Perid_1KmDispDist_strength_vec  = append(Perid_1KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Perid_1KmDispDist_deg_vec  = append(Perid_1KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Perid_1KmDispDist_habAv_vec  = append(Perid_1KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Perid_1KmDispDist_unw_b_c_vec  = append(Perid_1KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Perid_1KmDispDist_Patch_Area_vec  = append(Perid_1KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Perid_1KmDispDist_AUC_cv_vec = append(Perid_1KmDispDist_AUC_cv_vec, gbm_mod_Perid_1KmDispDist$cv.statistics$discrimination.mean)
  Perid_1KmDispDist_AUC_vec = append(Perid_1KmDispDist_AUC_vec, gbm_mod_Perid_1KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Perid_1KmDispDist$n.trees
  Perid_1KmDispDist_nt_vec = append(Perid_1KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Perid_1KmDispDist_predict_mat[,i] = predict(gbm_mod_Perid_1KmDispDist, Perid_1KmDispDist_stattest, gbm_mod_Perid_1KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Perid_1KmDispDist_output_tab[,"AUC_train"] = Perid_1KmDispDist_AUC_vec
Perid_1KmDispDist_output_tab[,"AUC_cv"] = Perid_1KmDispDist_AUC_cv_vec
Perid_1KmDispDist_output_tab[,"ntrees"] = Perid_1KmDispDist_nt_vec

#Var. importance columns
Perid_1KmDispDist_output_tab[,"HSI_imp"] = Perid_1KmDispDist_HSI_vec
Perid_1KmDispDist_output_tab[,"EgoSize_imp"] = Perid_1KmDispDist_EgoSize_vec
Perid_1KmDispDist_output_tab[,"strength_imp"] = Perid_1KmDispDist_strength_vec
Perid_1KmDispDist_output_tab[,"deg_imp"] = Perid_1KmDispDist_deg_vec
Perid_1KmDispDist_output_tab[,"habAv_imp"] = Perid_1KmDispDist_habAv_vec
Perid_1KmDispDist_output_tab[,"unw_b_c_imp"] = Perid_1KmDispDist_unw_b_c_vec
Perid_1KmDispDist_output_tab[,"Patch_Area_imp"] = Perid_1KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Perid_1KmDispDist_HSI_tab = data.frame(imp = Perid_1KmDispDist_HSI_vec)
Perid_1KmDispDist_EgoSize_tab = data.frame(imp = Perid_1KmDispDist_EgoSize_vec)
Perid_1KmDispDist_strength_tab = data.frame(imp = Perid_1KmDispDist_strength_vec)
Perid_1KmDispDist_deg_tab = data.frame(imp = Perid_1KmDispDist_deg_vec)
Perid_1KmDispDist_habAv_tab = data.frame(imp = Perid_1KmDispDist_habAv_vec)
Perid_1KmDispDist_unw_b_c_tab = data.frame(imp = Perid_1KmDispDist_unw_b_c_vec)
Perid_1KmDispDist_Patch_Area_tab = data.frame(imp = Perid_1KmDispDist_Patch_Area_vec)

Perid_1KmDispDist_HSI_tab$variable = "HSI"
Perid_1KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Perid_1KmDispDist_strength_tab$variable = "Strength"
Perid_1KmDispDist_deg_tab$variable = "Degree"
Perid_1KmDispDist_habAv_tab$variable = "Hab. Av."
Perid_1KmDispDist_unw_b_c_tab$variable = "B.C."
Perid_1KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Perid_1KmDispDist_AUC_cv_tab = data.frame(value = Perid_1KmDispDist_AUC_cv_vec)
Perid_1KmDispDist_AUC_train_tab = data.frame(value = Perid_1KmDispDist_AUC_vec)
#Make label of model for plot
Perid_1KmDispDist_AUC_cv_tab$model = "Perid_1KmDispDist"
Perid_1KmDispDist_AUC_train_tab$model = "Perid_1KmDispDist"

#combine pred. vars. into new data frame 
Perid_1KmDispDist_var_imp_tab = rbind(Perid_1KmDispDist_HSI_tab,Perid_1KmDispDist_EgoSize_tab,Perid_1KmDispDist_strength_tab,Perid_1KmDispDist_habAv_tab,Perid_1KmDispDist_deg_tab,Perid_1KmDispDist_unw_b_c_tab,Perid_1KmDispDist_Patch_Area_tab)

ggplot(Perid_1KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Perid_1KmDispDist_var_imp_tab$imp~Perid_1KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "1 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Perid_1KmDispDist_HSI_vec)
mean(Perid_1KmDispDist_EgoSize_vec)
mean(Perid_1KmDispDist_strength_vec)
mean(Perid_1KmDispDist_deg_vec)
mean(Perid_1KmDispDist_habAv_vec)
mean(Perid_1KmDispDist_unw_b_c_vec)
mean(Perid_1KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Perid_1KmDispDist_output_tab$AUC_cv)
summary(Perid_1KmDispDist_output_tab$AUC_train)


####################################### 6Km ####################################

# Create the output table that contains all the values of each of the runs
Perid_6KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid

# Create vectors to record performance measures
Perid_6KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Perid_6KmDispDist_AUC_vec = vector() #

Perid_6KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Perid_6KmDispDist_HSI_vec = vector() 
Perid_6KmDispDist_EgoSize_vec = vector() 
Perid_6KmDispDist_strength_vec  = vector()
Perid_6KmDispDist_deg_vec  = vector()
Perid_6KmDispDist_habAv_vec  = vector()
Perid_6KmDispDist_unw_b_c_vec  = vector()
Perid_6KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Perid_6KmDispDist_predict_mat = matrix(nrow = length(Perid_6KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Perid_6KmDispDist = gbm.step(data=Perid_6KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Perid_6KmDispDist)$var, imp = summary(gbm_mod_Perid_6KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Perid_6KmDispDist_HSI_vec = append(Perid_6KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Perid_6KmDispDist_EgoSize_vec = append(Perid_6KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Perid_6KmDispDist_strength_vec  = append(Perid_6KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Perid_6KmDispDist_deg_vec  = append(Perid_6KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Perid_6KmDispDist_habAv_vec  = append(Perid_6KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Perid_6KmDispDist_unw_b_c_vec  = append(Perid_6KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Perid_6KmDispDist_Patch_Area_vec  = append(Perid_6KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Perid_6KmDispDist_AUC_cv_vec = append(Perid_6KmDispDist_AUC_cv_vec, gbm_mod_Perid_6KmDispDist$cv.statistics$discrimination.mean)
  Perid_6KmDispDist_AUC_vec = append(Perid_6KmDispDist_AUC_vec, gbm_mod_Perid_6KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Perid_6KmDispDist$n.trees
  Perid_6KmDispDist_nt_vec = append(Perid_6KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Perid_6KmDispDist_predict_mat[,i] = predict(gbm_mod_Perid_6KmDispDist, Perid_6KmDispDist_stattest, gbm_mod_Perid_6KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Perid_6KmDispDist_output_tab[,"AUC_train"] = Perid_6KmDispDist_AUC_vec
Perid_6KmDispDist_output_tab[,"AUC_cv"] = Perid_6KmDispDist_AUC_cv_vec
Perid_6KmDispDist_output_tab[,"ntrees"] = Perid_6KmDispDist_nt_vec

#Var. importance columns
Perid_6KmDispDist_output_tab[,"HSI_imp"] = Perid_6KmDispDist_HSI_vec
Perid_6KmDispDist_output_tab[,"EgoSize_imp"] = Perid_6KmDispDist_EgoSize_vec
Perid_6KmDispDist_output_tab[,"strength_imp"] = Perid_6KmDispDist_strength_vec
Perid_6KmDispDist_output_tab[,"deg_imp"] = Perid_6KmDispDist_deg_vec
Perid_6KmDispDist_output_tab[,"habAv_imp"] = Perid_6KmDispDist_habAv_vec
Perid_6KmDispDist_output_tab[,"unw_b_c_imp"] = Perid_6KmDispDist_unw_b_c_vec
Perid_6KmDispDist_output_tab[,"Patch_Area_imp"] = Perid_6KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Perid_6KmDispDist_HSI_tab = data.frame(imp = Perid_6KmDispDist_HSI_vec)
Perid_6KmDispDist_EgoSize_tab = data.frame(imp = Perid_6KmDispDist_EgoSize_vec)
Perid_6KmDispDist_strength_tab = data.frame(imp = Perid_6KmDispDist_strength_vec)
Perid_6KmDispDist_deg_tab = data.frame(imp = Perid_6KmDispDist_deg_vec)
Perid_6KmDispDist_habAv_tab = data.frame(imp = Perid_6KmDispDist_habAv_vec)
Perid_6KmDispDist_unw_b_c_tab = data.frame(imp = Perid_6KmDispDist_unw_b_c_vec)
Perid_6KmDispDist_Patch_Area_tab = data.frame(imp = Perid_6KmDispDist_Patch_Area_vec)

Perid_6KmDispDist_HSI_tab$variable = "HSI"
Perid_6KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Perid_6KmDispDist_strength_tab$variable = "Strength"
Perid_6KmDispDist_deg_tab$variable = "Degree"
Perid_6KmDispDist_habAv_tab$variable = "Hab. Av."
Perid_6KmDispDist_unw_b_c_tab$variable = "B.C."
Perid_6KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Perid_6KmDispDist_AUC_cv_tab = data.frame(value = Perid_6KmDispDist_AUC_cv_vec)
Perid_6KmDispDist_AUC_train_tab = data.frame(value = Perid_6KmDispDist_AUC_vec)
#Make label of model for plot
Perid_6KmDispDist_AUC_cv_tab$model = "Perid_6KmDispDist"
Perid_6KmDispDist_AUC_train_tab$model = "Perid_6KmDispDist"

#combine pred. vars. into new data frame 
Perid_6KmDispDist_var_imp_tab = rbind(Perid_6KmDispDist_HSI_tab,Perid_6KmDispDist_EgoSize_tab,Perid_6KmDispDist_strength_tab,Perid_6KmDispDist_habAv_tab,Perid_6KmDispDist_deg_tab,Perid_6KmDispDist_unw_b_c_tab,Perid_6KmDispDist_Patch_Area_tab)

ggplot(Perid_6KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Perid_6KmDispDist_var_imp_tab$imp~Perid_6KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "6 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Perid_6KmDispDist_HSI_vec)
mean(Perid_6KmDispDist_EgoSize_vec)
mean(Perid_6KmDispDist_strength_vec)
mean(Perid_6KmDispDist_deg_vec)
mean(Perid_6KmDispDist_habAv_vec)
mean(Perid_6KmDispDist_unw_b_c_vec)
mean(Perid_6KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Perid_6KmDispDist_output_tab$AUC_cv)
summary(Perid_6KmDispDist_output_tab$AUC_train)


####################################### 4Km ####################################

# Create the output table that contains all the values of each of the runs
Perid_4KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid

# Create vectors to record performance measures
Perid_4KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Perid_4KmDispDist_AUC_vec = vector() #

Perid_4KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Perid_4KmDispDist_HSI_vec = vector()
Perid_4KmDispDist_EgoSize_vec = vector()
Perid_4KmDispDist_strength_vec  = vector()
Perid_4KmDispDist_deg_vec  = vector()
Perid_4KmDispDist_habAv_vec  = vector()
Perid_4KmDispDist_unw_b_c_vec  = vector()
Perid_4KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Perid_4KmDispDist_predict_mat = matrix(nrow = length(Perid_4KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Perid_4KmDispDist = gbm.step(data=Perid_4KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Perid_4KmDispDist)$var, imp = summary(gbm_mod_Perid_4KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Perid_4KmDispDist_HSI_vec = append(Perid_4KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Perid_4KmDispDist_EgoSize_vec = append(Perid_4KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Perid_4KmDispDist_strength_vec  = append(Perid_4KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Perid_4KmDispDist_deg_vec  = append(Perid_4KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Perid_4KmDispDist_habAv_vec  = append(Perid_4KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Perid_4KmDispDist_unw_b_c_vec  = append(Perid_4KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Perid_4KmDispDist_Patch_Area_vec  = append(Perid_4KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Perid_4KmDispDist_AUC_cv_vec = append(Perid_4KmDispDist_AUC_cv_vec, gbm_mod_Perid_4KmDispDist$cv.statistics$discrimination.mean)
  Perid_4KmDispDist_AUC_vec = append(Perid_4KmDispDist_AUC_vec, gbm_mod_Perid_4KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Perid_4KmDispDist$n.trees
  Perid_4KmDispDist_nt_vec = append(Perid_4KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Perid_4KmDispDist_predict_mat[,i] = predict(gbm_mod_Perid_4KmDispDist, Perid_4KmDispDist_stattest, gbm_mod_Perid_4KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Perid_4KmDispDist_output_tab[,"AUC_train"] = Perid_4KmDispDist_AUC_vec
Perid_4KmDispDist_output_tab[,"AUC_cv"] = Perid_4KmDispDist_AUC_cv_vec
Perid_4KmDispDist_output_tab[,"ntrees"] = Perid_4KmDispDist_nt_vec

#Var. importance columns
Perid_4KmDispDist_output_tab[,"HSI_imp"] = Perid_4KmDispDist_HSI_vec
Perid_4KmDispDist_output_tab[,"EgoSize_imp"] = Perid_4KmDispDist_EgoSize_vec
Perid_4KmDispDist_output_tab[,"strength_imp"] = Perid_4KmDispDist_strength_vec
Perid_4KmDispDist_output_tab[,"deg_imp"] = Perid_4KmDispDist_deg_vec
Perid_4KmDispDist_output_tab[,"habAv_imp"] = Perid_4KmDispDist_habAv_vec
Perid_4KmDispDist_output_tab[,"unw_b_c_imp"] = Perid_4KmDispDist_unw_b_c_vec
Perid_4KmDispDist_output_tab[,"Patch_Area_imp"] = Perid_4KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Perid_4KmDispDist_HSI_tab = data.frame(imp = Perid_4KmDispDist_HSI_vec)
Perid_4KmDispDist_EgoSize_tab = data.frame(imp = Perid_4KmDispDist_EgoSize_vec)
Perid_4KmDispDist_strength_tab = data.frame(imp = Perid_4KmDispDist_strength_vec)
Perid_4KmDispDist_deg_tab = data.frame(imp = Perid_4KmDispDist_deg_vec)
Perid_4KmDispDist_habAv_tab = data.frame(imp = Perid_4KmDispDist_habAv_vec)
Perid_4KmDispDist_unw_b_c_tab = data.frame(imp = Perid_4KmDispDist_unw_b_c_vec)
Perid_4KmDispDist_Patch_Area_tab = data.frame(imp = Perid_4KmDispDist_Patch_Area_vec)

Perid_4KmDispDist_HSI_tab$variable = "HSI"
Perid_4KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Perid_4KmDispDist_strength_tab$variable = "Strength"
Perid_4KmDispDist_deg_tab$variable = "Degree"
Perid_4KmDispDist_habAv_tab$variable = "Hab. Av."
Perid_4KmDispDist_unw_b_c_tab$variable = "B.C."
Perid_4KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Perid_4KmDispDist_AUC_cv_tab = data.frame(value = Perid_4KmDispDist_AUC_cv_vec)
Perid_4KmDispDist_AUC_train_tab = data.frame(value = Perid_4KmDispDist_AUC_vec)
#Make label of model for plot
Perid_4KmDispDist_AUC_cv_tab$model = "Perid_4KmDispDist"
Perid_4KmDispDist_AUC_train_tab$model = "Perid_4KmDispDist"

#combine pred. vars. into new data frame
Perid_4KmDispDist_var_imp_tab = rbind(Perid_4KmDispDist_HSI_tab,Perid_4KmDispDist_EgoSize_tab,Perid_4KmDispDist_strength_tab,Perid_4KmDispDist_habAv_tab,Perid_4KmDispDist_deg_tab,Perid_4KmDispDist_unw_b_c_tab,Perid_4KmDispDist_Patch_Area_tab)

ggplot(Perid_4KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Perid_4KmDispDist_var_imp_tab$imp~Perid_4KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "4 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars.
mean(Perid_4KmDispDist_HSI_vec)
mean(Perid_4KmDispDist_EgoSize_vec)
mean(Perid_4KmDispDist_strength_vec)
mean(Perid_4KmDispDist_deg_vec)
mean(Perid_4KmDispDist_habAv_vec)
mean(Perid_4KmDispDist_unw_b_c_vec)
mean(Perid_4KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy
summary(Perid_4KmDispDist_output_tab$AUC_cv)
summary(Perid_4KmDispDist_output_tab$AUC_train)


####################################### 8Km ####################################

# Create the output table that contains all the values of each of the runs
Perid_8KmDispDist_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid

# Create vectors to record performance measures
Perid_8KmDispDist_AUC_cv_vec = vector() #Cross-validated AUC
Perid_8KmDispDist_AUC_vec = vector() #

Perid_8KmDispDist_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Perid_8KmDispDist_HSI_vec = vector() 
Perid_8KmDispDist_EgoSize_vec = vector() 
Perid_8KmDispDist_strength_vec  = vector()
Perid_8KmDispDist_deg_vec  = vector()
Perid_8KmDispDist_habAv_vec  = vector()
Perid_8KmDispDist_unw_b_c_vec  = vector()
Perid_8KmDispDist_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Perid_8KmDispDist_predict_mat = matrix(nrow = length(Perid_8KmDispDist_stattest$PatchID), ncol = n_repeats, byrow = FALSE)

### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Perid_8KmDispDist = gbm.step(data=Perid_8KmDispDist_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Perid_8KmDispDist)$var, imp = summary(gbm_mod_Perid_8KmDispDist)$rel.inf)
  
  #Write var_imp results to the vectors
  Perid_8KmDispDist_HSI_vec = append(Perid_8KmDispDist_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Perid_8KmDispDist_EgoSize_vec = append(Perid_8KmDispDist_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Perid_8KmDispDist_strength_vec  = append(Perid_8KmDispDist_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Perid_8KmDispDist_deg_vec  = append(Perid_8KmDispDist_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Perid_8KmDispDist_habAv_vec  = append(Perid_8KmDispDist_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Perid_8KmDispDist_unw_b_c_vec  = append(Perid_8KmDispDist_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Perid_8KmDispDist_Patch_Area_vec  = append(Perid_8KmDispDist_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Perid_8KmDispDist_AUC_cv_vec = append(Perid_8KmDispDist_AUC_cv_vec, gbm_mod_Perid_8KmDispDist$cv.statistics$discrimination.mean)
  Perid_8KmDispDist_AUC_vec = append(Perid_8KmDispDist_AUC_vec, gbm_mod_Perid_8KmDispDist$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Perid_8KmDispDist$n.trees
  Perid_8KmDispDist_nt_vec = append(Perid_8KmDispDist_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Perid_8KmDispDist_predict_mat[,i] = predict(gbm_mod_Perid_8KmDispDist, Perid_8KmDispDist_stattest, gbm_mod_Perid_8KmDispDist$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}

Perid_8KmDispDist_output_tab[,"AUC_train"] = Perid_8KmDispDist_AUC_vec
Perid_8KmDispDist_output_tab[,"AUC_cv"] = Perid_8KmDispDist_AUC_cv_vec
Perid_8KmDispDist_output_tab[,"ntrees"] = Perid_8KmDispDist_nt_vec

#Var. importance columns
Perid_8KmDispDist_output_tab[,"HSI_imp"] = Perid_8KmDispDist_HSI_vec
Perid_8KmDispDist_output_tab[,"EgoSize_imp"] = Perid_8KmDispDist_EgoSize_vec
Perid_8KmDispDist_output_tab[,"strength_imp"] = Perid_8KmDispDist_strength_vec
Perid_8KmDispDist_output_tab[,"deg_imp"] = Perid_8KmDispDist_deg_vec
Perid_8KmDispDist_output_tab[,"habAv_imp"] = Perid_8KmDispDist_habAv_vec
Perid_8KmDispDist_output_tab[,"unw_b_c_imp"] = Perid_8KmDispDist_unw_b_c_vec
Perid_8KmDispDist_output_tab[,"Patch_Area_imp"] = Perid_8KmDispDist_Patch_Area_vec
#is.data.frame(output_tab)

# Make a dataframe for plotting overlaying histrograms in R
Perid_8KmDispDist_HSI_tab = data.frame(imp = Perid_8KmDispDist_HSI_vec)
Perid_8KmDispDist_EgoSize_tab = data.frame(imp = Perid_8KmDispDist_EgoSize_vec)
Perid_8KmDispDist_strength_tab = data.frame(imp = Perid_8KmDispDist_strength_vec)
Perid_8KmDispDist_deg_tab = data.frame(imp = Perid_8KmDispDist_deg_vec)
Perid_8KmDispDist_habAv_tab = data.frame(imp = Perid_8KmDispDist_habAv_vec)
Perid_8KmDispDist_unw_b_c_tab = data.frame(imp = Perid_8KmDispDist_unw_b_c_vec)
Perid_8KmDispDist_Patch_Area_tab = data.frame(imp = Perid_8KmDispDist_Patch_Area_vec)

Perid_8KmDispDist_HSI_tab$variable = "HSI"
Perid_8KmDispDist_EgoSize_tab$variable = "3rd. ord. neigh."
Perid_8KmDispDist_strength_tab$variable = "Strength"
Perid_8KmDispDist_deg_tab$variable = "Degree"
Perid_8KmDispDist_habAv_tab$variable = "Hab. Av."
Perid_8KmDispDist_unw_b_c_tab$variable = "B.C."
Perid_8KmDispDist_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Perid_8KmDispDist_AUC_cv_tab = data.frame(value = Perid_8KmDispDist_AUC_cv_vec)
Perid_8KmDispDist_AUC_train_tab = data.frame(value = Perid_8KmDispDist_AUC_vec)
#Make label of model for plot
Perid_8KmDispDist_AUC_cv_tab$model = "Perid_8KmDispDist"
Perid_8KmDispDist_AUC_train_tab$model = "Perid_8KmDispDist"

#combine pred. vars. into new data frame 
Perid_8KmDispDist_var_imp_tab = rbind(Perid_8KmDispDist_HSI_tab,Perid_8KmDispDist_EgoSize_tab,Perid_8KmDispDist_strength_tab,Perid_8KmDispDist_habAv_tab,Perid_8KmDispDist_deg_tab,Perid_8KmDispDist_unw_b_c_tab,Perid_8KmDispDist_Patch_Area_tab)

ggplot(Perid_8KmDispDist_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Perid_8KmDispDist_var_imp_tab$imp~Perid_8KmDispDist_var_imp_tab$variable,
        xlab = NULL, ylab= "Variable importance", 
        main = "8 km maximum dispersal distance",
        cex.axis = 1.25, cex.lab = 1.2)

#Get mean var. importance of all of the vars. 
mean(Perid_8KmDispDist_HSI_vec)
mean(Perid_8KmDispDist_EgoSize_vec)
mean(Perid_8KmDispDist_strength_vec)
mean(Perid_8KmDispDist_deg_vec)
mean(Perid_8KmDispDist_habAv_vec)
mean(Perid_8KmDispDist_unw_b_c_vec)
mean(Perid_8KmDispDist_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Perid_8KmDispDist_output_tab$AUC_cv)
summary(Perid_8KmDispDist_output_tab$AUC_train)


####################################################################################################
#### Compare scores between networks w/d0 variations of the same species ###########################
### Perid ##################

#Import evaluation df's of original run with species-specific dispersal distance
Perid_DefaultDispDist_AUC_cv_tab <- read.csv("Perid_DefaultDispDist_AUC_cv_tab.csv")
Perid_DefaultDispDist_AUC_train_tab <- read.csv("Perid_DefaultDispDist_AUC_train_tab.csv")
Perid_DefaultDispDist_var_imp_tab <- read.csv("Perid_DefaultDispDist_var_imp_tab.csv")
Perid_DefaultDispDist_output_tab <- read.csv("Perid_DefaultDispDist_BRToutput_tab.csv")
Perid_DefaultDispDist_predict_df <- read.csv("Perid_DefaultDispDist_predict_df.csv")

head(Perid_DefaultDispDist_AUC_cv_tab)
Perid_DefaultDispDist_AUC_cv_tab$X <- NULL
Perid_DefaultDispDist_AUC_train_tab$X <- NULL

####  Make a dataframe for plotting overlaying histograms in R
### cv AUC
Perid_AUC_cv_tab = rbind(Perid_DefaultDispDist_AUC_cv_tab, Perid_300mDispDist_AUC_cv_tab, Perid_1KmDispDist_AUC_cv_tab, Perid_2KmDispDist_AUC_cv_tab, 
                         Perid_4KmDispDist_AUC_cv_tab,
                         Perid_6KmDispDist_AUC_cv_tab, Perid_8KmDispDist_AUC_cv_tab, Perid_10KmDispDist_AUC_cv_tab)

#Change order of factors to display noTopo at the edge
Perid_AUC_cv_tab$model<- as.factor(Perid_AUC_cv_tab$model)
levels(Perid_AUC_cv_tab$model)
Perid_AUC_cv_tab$model<-factor(Perid_AUC_cv_tab$model, levels=c("Perid_300mDispDist", "Perid_1KmDispDist", "Perid", "Perid_2KmDispDist",
                                                                "Perid_4KmDispDist", 
                                                                "Perid_6KmDispDist", "Perid_8KmDispDist", "Perid_10KmDispDist"))

#Plot
ggplot(Perid_AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Perid_AUC_cv_tab$value~Perid_AUC_cv_tab$model, ylab= "Cross-validated AUC")

### train AUC
Perid_AUC_train_tab = rbind(Perid_DefaultDispDist_AUC_train_tab, Perid_300mDispDist_AUC_train_tab, Perid_1KmDispDist_AUC_train_tab, Perid_2KmDispDist_AUC_train_tab, 
                            Perid_4KmDispDist_AUC_train_tab,
                            Perid_6KmDispDist_AUC_train_tab, Perid_8KmDispDist_AUC_train_tab, Perid_10KmDispDist_AUC_train_tab)

levels(Perid_AUC_train_tab$model)
Perid_AUC_train_tab$model<-factor(Perid_AUC_train_tab$model, levels=c("Perid_300mDispDist", "Perid_1KmDispDist", "Perid", "Perid_2KmDispDist",
                                                                      "Perid_4KmDispDist", 
                                                                      "Perid_6KmDispDist", "Perid_8KmDispDist", "Perid_10KmDispDist"))

ggplot(Perid_AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(Perid_AUC_train_tab$value~Perid_AUC_train_tab$model, ylab= "Training AUC")

summary(Perid_DefaultDispDist_output_tab$AUC_cv)
summary(Perid_300mDispDist_output_tab$AUC_cv)
summary(Perid_1KmDispDist_output_tab$AUC_cv)
summary(Perid_2KmDispDist_output_tab$AUC_cv)
summary(Perid_4KmDispDist_output_tab$AUC_cv)
summary(Perid_6KmDispDist_output_tab$AUC_cv)
summary(Perid_8KmDispDist_output_tab$AUC_cv)
summary(Perid_10KmDispDist_output_tab$AUC_cv)

summary(Perid_DefaultDispDist_output_tab$AUC_train)
summary(Perid_300mDispDist_output_tab$AUC_train)
summary(Perid_1KmDispDist_output_tab$AUC_train)
summary(Perid_2KmDispDist_output_tab$AUC_train)
summary(Perid_4KmDispDist_output_tab$AUC_train)
summary(Perid_6KmDispDist_output_tab$AUC_train)
summary(Perid_8KmDispDist_output_tab$AUC_train)
summary(Perid_10KmDispDist_output_tab$AUC_train)

summary(Perid_DefaultDispDist_output_tab$ntrees)
summary(Perid_300mDispDist_output_tab$ntrees)
summary(Perid_1KmDispDist_output_tab$ntrees)
summary(Perid_2KmDispDist_output_tab$ntrees)
summary(Perid_4KmDispDist_output_tab$ntrees)
summary(Perid_6KmDispDist_output_tab$ntrees)
summary(Perid_8KmDispDist_output_tab$ntrees)
summary(Perid_10KmDispDist_output_tab$ntrees)




################################################################################
######## Plot response curves ##################################################
################################################################################

######### Get default disp. dist. var. imp. boxplots ###########################

boxplot(Bovar_DefaultDispDist_var_imp_tab$imp~Bovar_DefaultDispDist_var_imp_tab$variable, ylab= "Var. Importance", main = "Bovar_DefaultDispDist")
boxplot(Hyarb_DefaultDispDist_var_imp_tab$imp~Hyarb_DefaultDispDist_var_imp_tab$variable, ylab= "Var. Importance", main = "Hyarb_DefaultDispDist")
boxplot(Alobs_DefaultDispDist_var_imp_tab$imp~Alobs_DefaultDispDist_var_imp_tab$variable, ylab= "Var. Importance", main = "Alobs_DefaultDispDist")
boxplot(Epcal_DefaultDispDist_var_imp_tab$imp~Epcal_DefaultDispDist_var_imp_tab$variable, ylab= "Var. Importance", main = "Epcal_DefaultDispDist")
boxplot(Peagg_DefaultDispDist_var_imp_tab$imp~Peagg_DefaultDispDist_var_imp_tab$variable, ylab= "Var. Importance", main = "Peagg_DefaultDispDist")
boxplot(Perid_DefaultDispDist_var_imp_tab$imp~Perid_DefaultDispDist_var_imp_tab$variable, ylab= "Var. Importance", main = "Perid_DefaultDispDist")

###################################################################################################
#### Plot default (species-specificI) disp. dist. sample models' response curves 
#### (tagged with the kind of network they are) ####

### Import the response curves of the default (Sample) models from their original workspaces ######
setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")

gbm_mod_Hyarb_DefaultDispDist <- readRDS(file = "gbm_mod_Hyarb_DefaultDispDist.rds", refhook = NULL)
infoRDS("gbm_mod_Hyarb_DefaultDispDist.rds")
summary(gbm_mod_Hyarb_DefaultDispDist)
names(gbm_mod_Hyarb_DefaultDispDist)

gbm.plot.fits(gbm_mod_Hyarb_DefaultDispDist)
#plot different variables
gbm.plot(gbm_mod_Hyarb_DefaultDispDist, variable.no = 2, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,1), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)

gbm_mod_Alobs_DefaultDispDist <- readRDS(file = "gbm_mod_Alobs_DefaultDispDist.rds", refhook = NULL)
infoRDS("gbm_mod_Alobs_DefaultDispDist.rds")
summary(gbm_mod_Alobs_DefaultDispDist)

gbm.plot.fits(gbm_mod_Alobs_DefaultDispDist)
gbm.plot(gbm_mod_Alobs_DefaultDispDist, variable.no = 2, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,1), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)

setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
gbm_mod_Bovar_DefaultDispDist <- readRDS(file = "gbm_mod_Bovar_DefaultDispDist.rds", refhook = NULL)
infoRDS("gbm_mod_Bovar_DefaultDispDist.rds")
summary(gbm_mod_Bovar_DefaultDispDist)

setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
gbm_mod_Epcal_DefaultDispDist <- readRDS(file = "gbm_mod_Epcal_DefaultDispDist.rds", refhook = NULL)
infoRDS("gbm_mod_Epcal_DefaultDispDist.rds")
summary(gbm_mod_Epcal_DefaultDispDist)

setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
gbm_mod_Peagg_DefaultDispDist <- readRDS(file = "gbm_mod_Peagg_DefaultDispDist.rds", refhook = NULL)
infoRDS("gbm_mod_Peagg_DefaultDispDist.rds")
summary(gbm_mod_Peagg_DefaultDispDist)

setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
gbm_mod_Perid_DefaultDispDist <- readRDS(file = "gbm_mod_Perid_DefaultDispDist.rds", refhook = NULL)
infoRDS("gbm_mod_Perid_DefaultDispDist.rds")
summary(gbm_mod_Perid_DefaultDispDist)



################################################################################
#### Get Response curves of all the other models (the rest of the disp. dists.)#
################################################################################

### Turn the gbm.plot objects into dfs for alternetive/standalone presentation
# The code below is inspired by: https://stats.stackexchange.com/questions/122721/r-partial-dependency-plots-from-gbm-package-values-and-y-axis
# Name of predictors by number (changes in at least the default for Bovar)
## 1) Patch area
# 2) deg
# 3) unw_b_c
# 4) strength
# 5) habAv
# 6) EgoSize (3rd Order neighborhood)
# 7) HSI

#### Response curves for HSI, HabAv and EgoSize for Alytes obstetricans & Bombina variegata


################################################################################
#### Bovar #####################################################################


### 3rd Order neighborhood (EgoSize) ###########################################

#300m
gbm_mod_Bovar_300mDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Bovar_300mDispDist = plot(gbm_mod_Bovar_300mDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Bovar_300mDispDist$y = scale(resp_curve_EgoSize_Bovar_300mDispDist$y, scale = FALSE)
## standalone plot
# plot(resp_curve_Bovar_300mDispDist$EgoSize, resp_curve_Bovar_300mDispDist$y, main  = "Bovar_300mDispDist")
# lines(resp_curve_Bovar_300mDispDist$EgoSize, resp_curve_Bovar_300mDispDist$y)
## to Plot jointly
#Set identifier tags
resp_curve_EgoSize_Bovar_300mDispDist$model = "Bovar_300mDispDist"
resp_curve_EgoSize_Bovar_300mDispDist$Species = "Bovar"
resp_curve_EgoSize_Bovar_300mDispDist$Max_disp_dist_m = "300"
head(resp_curve_EgoSize_Bovar_300mDispDist)

# 1Km
gbm_mod_Bovar_1KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Bovar_1KmDispDist = plot(gbm_mod_Bovar_1KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Bovar_1KmDispDist$y = scale(resp_curve_EgoSize_Bovar_1KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Bovar_1KmDispDist$model = "Bovar_1KmDispDist"
resp_curve_EgoSize_Bovar_1KmDispDist$Species = "Bovar"
resp_curve_EgoSize_Bovar_1KmDispDist$Max_disp_dist_m = "1000"
head(resp_curve_EgoSize_Bovar_1KmDispDist)

# 2Km
gbm_mod_Bovar_2KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Bovar_2KmDispDist = plot(gbm_mod_Bovar_2KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Bovar_2KmDispDist$y = scale(resp_curve_EgoSize_Bovar_2KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Bovar_2KmDispDist$model = "Bovar_2KmDispDist"
resp_curve_EgoSize_Bovar_2KmDispDist$Species = "Bovar"
resp_curve_EgoSize_Bovar_2KmDispDist$Max_disp_dist_m = "2000"
head(resp_curve_EgoSize_Bovar_2KmDispDist)

# 4Km: species-specific maximum dispersal distance (ssmdd) for Bovar, so I don't run this
# gbm_mod_Bovar_4KmDispDist$var.names #query the model to get position of relevant variable
# resp_curve_EgoSize_Bovar_4KmDispDist = plot(gbm_mod_Bovar_4KmDispDist, i.var = 5, return.grid=TRUE)
# resp_curve_EgoSize_Bovar_4KmDispDist$y = scale(resp_curve_EgoSize_Bovar_4KmDispDist$y, scale = FALSE)
# #Set identifier tags
# resp_curve_EgoSize_Bovar_4KmDispDist$model = "Bovar_4KmDispDist"
# resp_curve_EgoSize_Bovar_4KmDispDist$Species = "Bovar"
# resp_curve_EgoSize_Bovar_4KmDispDist$Max_disp_dist_m = "1000"
# head(resp_curve_EgoSize_Bovar_4KmDispDist)

#Default (ssmdd)
gbm_mod_Bovar_DefaultDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Bovar_defaultDispDist = plot(gbm_mod_Bovar_DefaultDispDist, i.var = 6, return.grid=TRUE)
resp_curve_EgoSize_Bovar_defaultDispDist$y = scale(resp_curve_EgoSize_Bovar_defaultDispDist$y, scale = FALSE)
## standalone plot
# plot(resp_curve_Bovar_DefaultDispDist$Patch_Area, resp_curve_Bovar_DefaultDispDist$y, main  = "Bovar_DefaultDispDist")
# lines(resp_curve_Bovar_DefaultDispDist$Patch_Area, resp_curve_Bovar_DefaultDispDist$y)
## to Plot jointly
#Set identifier tags
resp_curve_EgoSize_Bovar_defaultDispDist$model = "Bovar_defaultDispDist"
resp_curve_EgoSize_Bovar_defaultDispDist$Species = "Bovar"
resp_curve_EgoSize_Bovar_defaultDispDist$Max_disp_dist_m = "4000"
head(resp_curve_EgoSize_Bovar_defaultDispDist)

# 6Km
gbm_mod_Bovar_6KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Bovar_6KmDispDist = plot(gbm_mod_Bovar_6KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Bovar_6KmDispDist$y = scale(resp_curve_EgoSize_Bovar_6KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Bovar_6KmDispDist$model = "Bovar_6KmDispDist"
resp_curve_EgoSize_Bovar_6KmDispDist$Species = "Bovar"
resp_curve_EgoSize_Bovar_6KmDispDist$Max_disp_dist_m = "6000"
head(resp_curve_EgoSize_Bovar_6KmDispDist)

# 8Km
gbm_mod_Bovar_8KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Bovar_8KmDispDist = plot(gbm_mod_Bovar_8KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Bovar_8KmDispDist$y = scale(resp_curve_EgoSize_Bovar_8KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Bovar_8KmDispDist$model = "Bovar_8KmDispDist"
resp_curve_EgoSize_Bovar_8KmDispDist$Species = "Bovar"
resp_curve_EgoSize_Bovar_8KmDispDist$Max_disp_dist_m = "8000"
head(resp_curve_EgoSize_Bovar_8KmDispDist)

# 10Km
gbm_mod_Bovar_10KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Bovar_10KmDispDist = plot(gbm_mod_Bovar_10KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Bovar_10KmDispDist$y = scale(resp_curve_EgoSize_Bovar_10KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Bovar_10KmDispDist$model = "Bovar_10KmDispDist"
resp_curve_EgoSize_Bovar_10KmDispDist$Species = "Bovar"
resp_curve_EgoSize_Bovar_10KmDispDist$Max_disp_dist_m = "10000"
head(resp_curve_EgoSize_Bovar_10KmDispDist)

##join to get joint response curves
resp_curve_EgoSize_Bovar_tab = rbind(resp_curve_EgoSize_Bovar_300mDispDist, resp_curve_EgoSize_Bovar_1KmDispDist,
                                     resp_curve_EgoSize_Bovar_2KmDispDist, resp_curve_EgoSize_Bovar_defaultDispDist,
                                     resp_curve_EgoSize_Bovar_6KmDispDist, resp_curve_EgoSize_Bovar_8KmDispDist,
                                     resp_curve_EgoSize_Bovar_10KmDispDist)

#Turn distances into numeric
resp_curve_EgoSize_Bovar_tab$Max_disp_dist_m <- as.numeric(resp_curve_EgoSize_Bovar_tab$Max_disp_dist_m)

#Reorder joint response curves tab
resp_curve_EgoSize_Bovar_tab$model <- as.factor(resp_curve_EgoSize_Bovar_tab$model)
resp_curve_EgoSize_Bovar_tab$model <-factor(resp_curve_EgoSize_Bovar_tab$model, 
                                            levels=c("Bovar_300mDispDist", "Bovar_1KmDispDist", "Bovar_2KmDispDist",
                                                     "Bovar_defaultDispDist", "Bovar_6KmDispDist", 
                                                     "Bovar_8KmDispDist","Bovar_10KmDispDist"))

#Plot
ggplot(resp_curve_EgoSize_Bovar_tab, aes(x=EgoSize, y = y, color = model, group = model)) +
  geom_point(size = 2) + geom_line(size = 0.5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "3rd order neighborhood", y= "Occurrence-state")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))

ggplot(resp_curve_EgoSize_Bovar_tab, aes(x=EgoSize, y = y, color = Max_disp_dist_m, group = Max_disp_dist_m)) +
  geom_point() + geom_line()


### HSI ########################################################################

#300m
gbm_mod_Bovar_300mDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Bovar_300mDispDist = plot(gbm_mod_Bovar_300mDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Bovar_300mDispDist$y = scale(resp_curve_HSI_Bovar_300mDispDist$y, scale = FALSE)
## standalone plot
# plot(resp_curve_Bovar_300mDispDist$HSI, resp_curve_Bovar_300mDispDist$y, main  = "Bovar_300mDispDist")
# lines(resp_curve_Bovar_300mDispDist$HSI, resp_curve_Bovar_300mDispDist$y)
## to Plot jointly
#Set identifier tags
resp_curve_HSI_Bovar_300mDispDist$model = "Bovar_300mDispDist"
resp_curve_HSI_Bovar_300mDispDist$Species = "Bovar"
resp_curve_HSI_Bovar_300mDispDist$Max_disp_dist_m = "300"
head(resp_curve_HSI_Bovar_300mDispDist)

# 1Km
gbm_mod_Bovar_1KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Bovar_1KmDispDist = plot(gbm_mod_Bovar_1KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Bovar_1KmDispDist$y = scale(resp_curve_HSI_Bovar_1KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Bovar_1KmDispDist$model = "Bovar_1KmDispDist"
resp_curve_HSI_Bovar_1KmDispDist$Species = "Bovar"
resp_curve_HSI_Bovar_1KmDispDist$Max_disp_dist_m = "1000"
head(resp_curve_HSI_Bovar_1KmDispDist)

# 2Km
gbm_mod_Bovar_2KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Bovar_2KmDispDist = plot(gbm_mod_Bovar_2KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Bovar_2KmDispDist$y = scale(resp_curve_HSI_Bovar_2KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Bovar_2KmDispDist$model = "Bovar_2KmDispDist"
resp_curve_HSI_Bovar_2KmDispDist$Species = "Bovar"
resp_curve_HSI_Bovar_2KmDispDist$Max_disp_dist_m = "2000"
head(resp_curve_HSI_Bovar_2KmDispDist)

# 4Km: ssmdd for Bovar
# gbm_mod_Bovar_4KmDispDist$var.names #query the model to get position of relevant variable
# resp_curve_HSI_Bovar_4KmDispDist = plot(gbm_mod_Bovar_4KmDispDist, i.var = 5, return.grid=TRUE)
# resp_curve_HSI_Bovar_4KmDispDist$y = scale(resp_curve_HSI_Bovar_4KmDispDist$y, scale = FALSE)
# #Set identifier tags
# resp_curve_HSI_Bovar_4KmDispDist$model = "Bovar_4KmDispDist"
# resp_curve_HSI_Bovar_4KmDispDist$Species = "Bovar"
# resp_curve_HSI_Bovar_4KmDispDist$Max_disp_dist_m = "1000"
# head(resp_curve_HSI_Bovar_4KmDispDist)

#Default
gbm_mod_Bovar_DefaultDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Bovar_defaultDispDist = plot(gbm_mod_Bovar_DefaultDispDist, i.var = 7, return.grid=TRUE)
resp_curve_HSI_Bovar_defaultDispDist$y = scale(resp_curve_HSI_Bovar_defaultDispDist$y, scale = FALSE)
## standalone plot
# plot(resp_curve_Bovar_DefaultDispDist$Patch_Area, resp_curve_Bovar_DefaultDispDist$y, main  = "Bovar_DefaultDispDist")
# lines(resp_curve_Bovar_DefaultDispDist$Patch_Area, resp_curve_Bovar_DefaultDispDist$y)
## to Plot jointly
#Set identifier tags
resp_curve_HSI_Bovar_defaultDispDist$model = "Bovar_defaultDispDist"
resp_curve_HSI_Bovar_defaultDispDist$Species = "Bovar"
resp_curve_HSI_Bovar_defaultDispDist$Max_disp_dist_m = "4000"
head(resp_curve_HSI_Bovar_defaultDispDist)

# 6Km
gbm_mod_Bovar_6KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Bovar_6KmDispDist = plot(gbm_mod_Bovar_6KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Bovar_6KmDispDist$y = scale(resp_curve_HSI_Bovar_6KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Bovar_6KmDispDist$model = "Bovar_6KmDispDist"
resp_curve_HSI_Bovar_6KmDispDist$Species = "Bovar"
resp_curve_HSI_Bovar_6KmDispDist$Max_disp_dist_m = "6000"
head(resp_curve_HSI_Bovar_6KmDispDist)

# 8Km
gbm_mod_Bovar_8KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Bovar_8KmDispDist = plot(gbm_mod_Bovar_8KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Bovar_8KmDispDist$y = scale(resp_curve_HSI_Bovar_8KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Bovar_8KmDispDist$model = "Bovar_8KmDispDist"
resp_curve_HSI_Bovar_8KmDispDist$Species = "Bovar"
resp_curve_HSI_Bovar_8KmDispDist$Max_disp_dist_m = "8000"
head(resp_curve_HSI_Bovar_8KmDispDist)

# 10Km
gbm_mod_Bovar_10KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Bovar_10KmDispDist = plot(gbm_mod_Bovar_10KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Bovar_10KmDispDist$y = scale(resp_curve_HSI_Bovar_10KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Bovar_10KmDispDist$model = "Bovar_10KmDispDist"
resp_curve_HSI_Bovar_10KmDispDist$Species = "Bovar"
resp_curve_HSI_Bovar_10KmDispDist$Max_disp_dist_m = "10000"
head(resp_curve_HSI_Bovar_10KmDispDist)

##join to get joint response curves
resp_curve_HSI_Bovar_tab = rbind(resp_curve_HSI_Bovar_300mDispDist, resp_curve_HSI_Bovar_1KmDispDist,
                                 resp_curve_HSI_Bovar_2KmDispDist, resp_curve_HSI_Bovar_defaultDispDist,
                                 resp_curve_HSI_Bovar_6KmDispDist, resp_curve_HSI_Bovar_8KmDispDist,
                                 resp_curve_HSI_Bovar_10KmDispDist)

#Turn distances into numeric
resp_curve_HSI_Bovar_tab$Max_disp_dist_m <- as.numeric(resp_curve_HSI_Bovar_tab$Max_disp_dist_m)

#Reorder joint response curves tab
resp_curve_HSI_Bovar_tab$model <- as.factor(resp_curve_HSI_Bovar_tab$model)
resp_curve_HSI_Bovar_tab$model <-factor(resp_curve_HSI_Bovar_tab$model, 
                                        levels=c("Bovar_300mDispDist", "Bovar_1KmDispDist", "Bovar_2KmDispDist",
                                                 "Bovar_defaultDispDist", "Bovar_6KmDispDist", 
                                                 "Bovar_8KmDispDist","Bovar_10KmDispDist"))

#Plot
ggplot(resp_curve_HSI_Bovar_tab, aes(x=HSI, y = y, color = model, group = model)) +
  geom_point(size = 2) + geom_line(size = 0.5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Habitat suitability index", y= "Occurrence-state")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))

ggplot(resp_curve_HSI_Bovar_tab, aes(x=HSI, y = y, color = Max_disp_dist_m, group = Max_disp_dist_m)) +
  geom_point() + geom_line()


### habAv ######################################################################

#300m
gbm_mod_Bovar_300mDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Bovar_300mDispDist = plot(gbm_mod_Bovar_300mDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Bovar_300mDispDist$y = scale(resp_curve_habAv_Bovar_300mDispDist$y, scale = FALSE)
## standalone plot
# plot(resp_curve_Bovar_300mDispDist$habAv, resp_curve_Bovar_300mDispDist$y, main  = "Bovar_300mDispDist")
# lines(resp_curve_Bovar_300mDispDist$habAv, resp_curve_Bovar_300mDispDist$y)
## to Plot jointly
#Set identifier tags
resp_curve_habAv_Bovar_300mDispDist$model = "Bovar_300mDispDist"
resp_curve_habAv_Bovar_300mDispDist$Species = "Bovar"
resp_curve_habAv_Bovar_300mDispDist$Max_disp_dist_m = "300"
head(resp_curve_habAv_Bovar_300mDispDist)

# 1Km
gbm_mod_Bovar_1KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Bovar_1KmDispDist = plot(gbm_mod_Bovar_1KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Bovar_1KmDispDist$y = scale(resp_curve_habAv_Bovar_1KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Bovar_1KmDispDist$model = "Bovar_1KmDispDist"
resp_curve_habAv_Bovar_1KmDispDist$Species = "Bovar"
resp_curve_habAv_Bovar_1KmDispDist$Max_disp_dist_m = "1000"
head(resp_curve_habAv_Bovar_1KmDispDist)

# 2Km
gbm_mod_Bovar_2KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Bovar_2KmDispDist = plot(gbm_mod_Bovar_2KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Bovar_2KmDispDist$y = scale(resp_curve_habAv_Bovar_2KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Bovar_2KmDispDist$model = "Bovar_2KmDispDist"
resp_curve_habAv_Bovar_2KmDispDist$Species = "Bovar"
resp_curve_habAv_Bovar_2KmDispDist$Max_disp_dist_m = "2000"
head(resp_curve_habAv_Bovar_2KmDispDist)

# 4Km: ssmdd for Bovar
# gbm_mod_Bovar_4KmDispDist$var.names #query the model to get position of relevant variable
# resp_curve_habAv_Bovar_4KmDispDist = plot(gbm_mod_Bovar_4KmDispDist, i.var = 5, return.grid=TRUE)
# resp_curve_habAv_Bovar_4KmDispDist$y = scale(resp_curve_habAv_Bovar_4KmDispDist$y, scale = FALSE)
# #Set identifier tags
# resp_curve_habAv_Bovar_4KmDispDist$model = "Bovar_4KmDispDist"
# resp_curve_habAv_Bovar_4KmDispDist$Species = "Bovar"
# resp_curve_habAv_Bovar_4KmDispDist$Max_disp_dist_m = "1000"
# head(resp_curve_habAv_Bovar_4KmDispDist)

#Default
gbm_mod_Bovar_DefaultDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Bovar_defaultDispDist = plot(gbm_mod_Bovar_DefaultDispDist, i.var = 5, return.grid=TRUE)
resp_curve_habAv_Bovar_defaultDispDist$y = scale(resp_curve_habAv_Bovar_defaultDispDist$y, scale = FALSE)
## standalone plot
# plot(resp_curve_Bovar_DefaultDispDist$Patch_Area, resp_curve_Bovar_DefaultDispDist$y, main  = "Bovar_DefaultDispDist")
# lines(resp_curve_Bovar_DefaultDispDist$Patch_Area, resp_curve_Bovar_DefaultDispDist$y)
## to Plot jointly
#Set identifier tags
resp_curve_habAv_Bovar_defaultDispDist$model = "Bovar_defaultDispDist"
resp_curve_habAv_Bovar_defaultDispDist$Species = "Bovar"
resp_curve_habAv_Bovar_defaultDispDist$Max_disp_dist_m = "4000"
head(resp_curve_habAv_Bovar_defaultDispDist)

# 6Km
gbm_mod_Bovar_6KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Bovar_6KmDispDist = plot(gbm_mod_Bovar_6KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Bovar_6KmDispDist$y = scale(resp_curve_habAv_Bovar_6KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Bovar_6KmDispDist$model = "Bovar_6KmDispDist"
resp_curve_habAv_Bovar_6KmDispDist$Species = "Bovar"
resp_curve_habAv_Bovar_6KmDispDist$Max_disp_dist_m = "6000"
head(resp_curve_habAv_Bovar_6KmDispDist)

# 8Km
gbm_mod_Bovar_8KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Bovar_8KmDispDist = plot(gbm_mod_Bovar_8KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Bovar_8KmDispDist$y = scale(resp_curve_habAv_Bovar_8KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Bovar_8KmDispDist$model = "Bovar_8KmDispDist"
resp_curve_habAv_Bovar_8KmDispDist$Species = "Bovar"
resp_curve_habAv_Bovar_8KmDispDist$Max_disp_dist_m = "8000"
head(resp_curve_habAv_Bovar_8KmDispDist)

# 10Km
gbm_mod_Bovar_10KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Bovar_10KmDispDist = plot(gbm_mod_Bovar_10KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Bovar_10KmDispDist$y = scale(resp_curve_habAv_Bovar_10KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Bovar_10KmDispDist$model = "Bovar_10KmDispDist"
resp_curve_habAv_Bovar_10KmDispDist$Species = "Bovar"
resp_curve_habAv_Bovar_10KmDispDist$Max_disp_dist_m = "10000"
head(resp_curve_habAv_Bovar_10KmDispDist)

##join to get joint response curves
resp_curve_habAv_Bovar_tab = rbind(resp_curve_habAv_Bovar_300mDispDist, resp_curve_habAv_Bovar_1KmDispDist,
                                   resp_curve_habAv_Bovar_2KmDispDist, resp_curve_habAv_Bovar_defaultDispDist,
                                   resp_curve_habAv_Bovar_6KmDispDist, resp_curve_habAv_Bovar_8KmDispDist,
                                   resp_curve_habAv_Bovar_10KmDispDist)

#Turn distances into numeric
resp_curve_habAv_Bovar_tab$Max_disp_dist_m <- as.numeric(resp_curve_habAv_Bovar_tab$Max_disp_dist_m)

#Reorder joint response curves tab
resp_curve_habAv_Bovar_tab$model <- as.factor(resp_curve_habAv_Bovar_tab$model)
resp_curve_habAv_Bovar_tab$model <-factor(resp_curve_habAv_Bovar_tab$model, 
                                          levels=c("Bovar_300mDispDist", "Bovar_1KmDispDist", "Bovar_2KmDispDist",
                                                   "Bovar_defaultDispDist", "Bovar_6KmDispDist", 
                                                   "Bovar_8KmDispDist","Bovar_10KmDispDist"))

#Plot
ggplot(resp_curve_habAv_Bovar_tab, aes(x=habAv, y = y, color = model, group = model)) +
  geom_point(size = 2) + geom_line(size = 0.5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Habitat availability", y= "Occurrence-state")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))

ggplot(resp_curve_habAv_Bovar_tab, aes(x=habAv, y = y, color = Max_disp_dist_m, group = Max_disp_dist_m)) +
  geom_point() + geom_line()



################################################################################
### Alobs ######################################################################

## 3rd Order neighborhood (EgoSize) ############################################

#300m
gbm_mod_Alobs_300mDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Alobs_300mDispDist = plot(gbm_mod_Alobs_300mDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Alobs_300mDispDist$y = scale(resp_curve_EgoSize_Alobs_300mDispDist$y, scale = FALSE)
## to Plot jointly
#Set identifier tags
resp_curve_EgoSize_Alobs_300mDispDist$model = "Alobs_300mDispDist"
resp_curve_EgoSize_Alobs_300mDispDist$Species = "Alobs"
resp_curve_EgoSize_Alobs_300mDispDist$Max_disp_dist_m = "300"
head(resp_curve_EgoSize_Alobs_300mDispDist)

# 1Km
gbm_mod_Alobs_1KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Alobs_1KmDispDist = plot(gbm_mod_Alobs_1KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Alobs_1KmDispDist$y = scale(resp_curve_EgoSize_Alobs_1KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Alobs_1KmDispDist$model = "Alobs_1KmDispDist"
resp_curve_EgoSize_Alobs_1KmDispDist$Species = "Alobs"
resp_curve_EgoSize_Alobs_1KmDispDist$Max_disp_dist_m = "1000"
head(resp_curve_EgoSize_Alobs_1KmDispDist)

#Default
gbm_mod_Alobs_DefaultDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Alobs_defaultDispDist = plot(gbm_mod_Alobs_DefaultDispDist, i.var = 6, return.grid=TRUE)
resp_curve_EgoSize_Alobs_defaultDispDist$y = scale(resp_curve_EgoSize_Alobs_defaultDispDist$y, scale = FALSE)
## to Plot jointly
#Set identifier tags
resp_curve_EgoSize_Alobs_defaultDispDist$model = "Alobs_defaultDispDist"
resp_curve_EgoSize_Alobs_defaultDispDist$Species = "Alobs"
resp_curve_EgoSize_Alobs_defaultDispDist$Max_disp_dist_m = "1500"
head(resp_curve_EgoSize_Alobs_defaultDispDist)

# 2Km
gbm_mod_Alobs_2KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Alobs_2KmDispDist = plot(gbm_mod_Alobs_2KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Alobs_2KmDispDist$y = scale(resp_curve_EgoSize_Alobs_2KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Alobs_2KmDispDist$model = "Alobs_2KmDispDist"
resp_curve_EgoSize_Alobs_2KmDispDist$Species = "Alobs"
resp_curve_EgoSize_Alobs_2KmDispDist$Max_disp_dist_m = "2000"
head(resp_curve_EgoSize_Alobs_2KmDispDist)

# 4km
gbm_mod_Alobs_4KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Alobs_4KmDispDist = plot(gbm_mod_Alobs_4KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Alobs_4KmDispDist$y = scale(resp_curve_EgoSize_Alobs_4KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Alobs_4KmDispDist$model = "Alobs_4KmDispDist"
resp_curve_EgoSize_Alobs_4KmDispDist$Species = "Alobs"
resp_curve_EgoSize_Alobs_4KmDispDist$Max_disp_dist_m = "4000"
head(resp_curve_EgoSize_Alobs_4KmDispDist)

# 6Km
gbm_mod_Alobs_6KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Alobs_6KmDispDist = plot(gbm_mod_Alobs_6KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Alobs_6KmDispDist$y = scale(resp_curve_EgoSize_Alobs_6KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Alobs_6KmDispDist$model = "Alobs_6KmDispDist"
resp_curve_EgoSize_Alobs_6KmDispDist$Species = "Alobs"
resp_curve_EgoSize_Alobs_6KmDispDist$Max_disp_dist_m = "6000"
head(resp_curve_EgoSize_Alobs_6KmDispDist)

# 8Km
gbm_mod_Alobs_8KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Alobs_8KmDispDist = plot(gbm_mod_Alobs_8KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Alobs_8KmDispDist$y = scale(resp_curve_EgoSize_Alobs_8KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Alobs_8KmDispDist$model = "Alobs_8KmDispDist"
resp_curve_EgoSize_Alobs_8KmDispDist$Species = "Alobs"
resp_curve_EgoSize_Alobs_8KmDispDist$Max_disp_dist_m = "8000"
head(resp_curve_EgoSize_Alobs_8KmDispDist)

# 10Km
gbm_mod_Alobs_10KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_EgoSize_Alobs_10KmDispDist = plot(gbm_mod_Alobs_10KmDispDist, i.var = 5, return.grid=TRUE)
resp_curve_EgoSize_Alobs_10KmDispDist$y = scale(resp_curve_EgoSize_Alobs_10KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_EgoSize_Alobs_10KmDispDist$model = "Alobs_10KmDispDist"
resp_curve_EgoSize_Alobs_10KmDispDist$Species = "Alobs"
resp_curve_EgoSize_Alobs_10KmDispDist$Max_disp_dist_m = "10000"
head(resp_curve_EgoSize_Alobs_10KmDispDist)

##join to get joint response curves
resp_curve_EgoSize_Alobs_tab = rbind(resp_curve_EgoSize_Alobs_300mDispDist, resp_curve_EgoSize_Alobs_1KmDispDist,
                                     resp_curve_EgoSize_Alobs_defaultDispDist, resp_curve_EgoSize_Alobs_2KmDispDist, resp_curve_EgoSize_Alobs_4KmDispDist, 
                                     resp_curve_EgoSize_Alobs_4KmDispDist, resp_curve_EgoSize_Alobs_6KmDispDist, 
                                     resp_curve_EgoSize_Alobs_8KmDispDist, resp_curve_EgoSize_Alobs_10KmDispDist)

#Turn distances into numeric
resp_curve_EgoSize_Alobs_tab$Max_disp_dist_m <- as.numeric(resp_curve_EgoSize_Alobs_tab$Max_disp_dist_m)

#Reorder joint response curves tab
resp_curve_EgoSize_Alobs_tab$model <- as.factor(resp_curve_EgoSize_Alobs_tab$model)
resp_curve_EgoSize_Alobs_tab$model <-factor(resp_curve_EgoSize_Alobs_tab$model, 
                                            levels=c("Alobs_300mDispDist", "Alobs_1KmDispDist", "Alobs_defaultDispDist",
                                                     "Alobs_2KmDispDist", "Alobs_4KmDispDist", "Alobs_6KmDispDist", 
                                                     "Alobs_8KmDispDist","Alobs_10KmDispDist")) 

#Plot
ggplot(resp_curve_EgoSize_Alobs_tab, aes(x=EgoSize, y = y, color = model, group = model)) +
  geom_point(size = 2) + geom_line(size = 0.5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "3rd order neighborhood", y= "Occurrence-state")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))

ggplot(resp_curve_EgoSize_Alobs_tab, aes(x=EgoSize, y = y, color = Max_disp_dist_m, group = Max_disp_dist_m)) +
  geom_point() + geom_line()


## HSI #########################################################################

#300m
gbm_mod_Alobs_300mDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Alobs_300mDispDist = plot(gbm_mod_Alobs_300mDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Alobs_300mDispDist$y = scale(resp_curve_HSI_Alobs_300mDispDist$y, scale = FALSE)
## to Plot jointly
#Set identifier tags
resp_curve_HSI_Alobs_300mDispDist$model = "Alobs_300mDispDist"
resp_curve_HSI_Alobs_300mDispDist$Species = "Alobs"
resp_curve_HSI_Alobs_300mDispDist$Max_disp_dist_m = "300"
head(resp_curve_HSI_Alobs_300mDispDist)

# 1Km
gbm_mod_Alobs_1KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Alobs_1KmDispDist = plot(gbm_mod_Alobs_1KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Alobs_1KmDispDist$y = scale(resp_curve_HSI_Alobs_1KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Alobs_1KmDispDist$model = "Alobs_1KmDispDist"
resp_curve_HSI_Alobs_1KmDispDist$Species = "Alobs"
resp_curve_HSI_Alobs_1KmDispDist$Max_disp_dist_m = "1000"
head(resp_curve_HSI_Alobs_1KmDispDist)

#Default
gbm_mod_Alobs_DefaultDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Alobs_defaultDispDist = plot(gbm_mod_Alobs_DefaultDispDist, i.var = 7, return.grid=TRUE)
resp_curve_HSI_Alobs_defaultDispDist$y = scale(resp_curve_HSI_Alobs_defaultDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Alobs_defaultDispDist$model = "Alobs_defaultDispDist"
resp_curve_HSI_Alobs_defaultDispDist$Species = "Alobs"
resp_curve_HSI_Alobs_defaultDispDist$Max_disp_dist_m = "1500"
head(resp_curve_HSI_Alobs_defaultDispDist)

# 2Km
gbm_mod_Alobs_2KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Alobs_2KmDispDist = plot(gbm_mod_Alobs_2KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Alobs_2KmDispDist$y = scale(resp_curve_HSI_Alobs_2KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Alobs_2KmDispDist$model = "Alobs_2KmDispDist"
resp_curve_HSI_Alobs_2KmDispDist$Species = "Alobs"
resp_curve_HSI_Alobs_2KmDispDist$Max_disp_dist_m = "2000"
head(resp_curve_HSI_Alobs_2KmDispDist)

# 4km
gbm_mod_Alobs_4KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Alobs_4KmDispDist = plot(gbm_mod_Alobs_4KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Alobs_4KmDispDist$y = scale(resp_curve_HSI_Alobs_4KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Alobs_4KmDispDist$model = "Alobs_4KmDispDist"
resp_curve_HSI_Alobs_4KmDispDist$Species = "Alobs"
resp_curve_HSI_Alobs_4KmDispDist$Max_disp_dist_m = "1000"
head(resp_curve_HSI_Alobs_4KmDispDist)

# 6Km
gbm_mod_Alobs_6KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Alobs_6KmDispDist = plot(gbm_mod_Alobs_6KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Alobs_6KmDispDist$y = scale(resp_curve_HSI_Alobs_6KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Alobs_6KmDispDist$model = "Alobs_6KmDispDist"
resp_curve_HSI_Alobs_6KmDispDist$Species = "Alobs"
resp_curve_HSI_Alobs_6KmDispDist$Max_disp_dist_m = "6000"
head(resp_curve_HSI_Alobs_6KmDispDist)

# 8Km
gbm_mod_Alobs_8KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Alobs_8KmDispDist = plot(gbm_mod_Alobs_8KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Alobs_8KmDispDist$y = scale(resp_curve_HSI_Alobs_8KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Alobs_8KmDispDist$model = "Alobs_8KmDispDist"
resp_curve_HSI_Alobs_8KmDispDist$Species = "Alobs"
resp_curve_HSI_Alobs_8KmDispDist$Max_disp_dist_m = "8000"
head(resp_curve_HSI_Alobs_8KmDispDist)

# 10Km
gbm_mod_Alobs_10KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_HSI_Alobs_10KmDispDist = plot(gbm_mod_Alobs_10KmDispDist, i.var = 6, return.grid=TRUE)
resp_curve_HSI_Alobs_10KmDispDist$y = scale(resp_curve_HSI_Alobs_10KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_HSI_Alobs_10KmDispDist$model = "Alobs_10KmDispDist"
resp_curve_HSI_Alobs_10KmDispDist$Species = "Alobs"
resp_curve_HSI_Alobs_10KmDispDist$Max_disp_dist_m = "10000"
head(resp_curve_HSI_Alobs_10KmDispDist)

##join to get joint response curves
resp_curve_HSI_Alobs_tab = rbind(resp_curve_HSI_Alobs_300mDispDist, resp_curve_HSI_Alobs_1KmDispDist,
                                 resp_curve_HSI_Alobs_defaultDispDist, resp_curve_HSI_Alobs_2KmDispDist, 
                                 resp_curve_HSI_Alobs_4KmDispDist, resp_curve_HSI_Alobs_6KmDispDist,
                                 resp_curve_HSI_Alobs_8KmDispDist, resp_curve_HSI_Alobs_10KmDispDist)

#Turn distances into numeric
resp_curve_HSI_Alobs_tab$Max_disp_dist_m <- as.numeric(resp_curve_HSI_Alobs_tab$Max_disp_dist_m)

#Reorder joint response curves tab
resp_curve_HSI_Alobs_tab$model <- as.factor(resp_curve_HSI_Alobs_tab$model)
resp_curve_HSI_Alobs_tab$model <-factor(resp_curve_HSI_Alobs_tab$model, 
                                        levels=c("Alobs_300mDispDist", "Alobs_1KmDispDist", "Alobs_defaultDispDist",
                                                 "Alobs_2KmDispDist", "Alobs_4KmDispDist", "Alobs_6KmDispDist", 
                                                 "Alobs_8KmDispDist","Alobs_10KmDispDist")) 

#Plot
ggplot(resp_curve_HSI_Alobs_tab, aes(x=HSI, y = y, color = model, group = model)) +
  geom_point(size = 2) + geom_line(size = 0.5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Habitat suitability index", y= "Occurrence-state")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))

ggplot(resp_curve_HSI_Alobs_tab, aes(x=HSI, y = y, color = Max_disp_dist_m, group = Max_disp_dist_m)) +
  geom_point() + geom_line()


### habAv ######################################################################

#300m
gbm_mod_Alobs_300mDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Alobs_300mDispDist = plot(gbm_mod_Alobs_300mDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Alobs_300mDispDist$y = scale(resp_curve_habAv_Alobs_300mDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Alobs_300mDispDist$model = "Alobs_300mDispDist"
resp_curve_habAv_Alobs_300mDispDist$Species = "Alobs"
resp_curve_habAv_Alobs_300mDispDist$Max_disp_dist_m = "300"
head(resp_curve_habAv_Alobs_300mDispDist)

# 1Km
gbm_mod_Alobs_1KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Alobs_1KmDispDist = plot(gbm_mod_Alobs_1KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Alobs_1KmDispDist$y = scale(resp_curve_habAv_Alobs_1KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Alobs_1KmDispDist$model = "Alobs_1KmDispDist"
resp_curve_habAv_Alobs_1KmDispDist$Species = "Alobs"
resp_curve_habAv_Alobs_1KmDispDist$Max_disp_dist_m = "1000"
head(resp_curve_habAv_Alobs_1KmDispDist)

#Default
gbm_mod_Alobs_DefaultDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Alobs_defaultDispDist = plot(gbm_mod_Alobs_DefaultDispDist, i.var = 5, return.grid=TRUE)
resp_curve_habAv_Alobs_defaultDispDist$y = scale(resp_curve_habAv_Alobs_defaultDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Alobs_defaultDispDist$model = "Alobs_defaultDispDist"
resp_curve_habAv_Alobs_defaultDispDist$Species = "Alobs"
resp_curve_habAv_Alobs_defaultDispDist$Max_disp_dist_m = "1500"
head(resp_curve_habAv_Alobs_defaultDispDist)

# 2Km
gbm_mod_Alobs_2KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Alobs_2KmDispDist = plot(gbm_mod_Alobs_2KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Alobs_2KmDispDist$y = scale(resp_curve_habAv_Alobs_2KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Alobs_2KmDispDist$model = "Alobs_2KmDispDist"
resp_curve_habAv_Alobs_2KmDispDist$Species = "Alobs"
resp_curve_habAv_Alobs_2KmDispDist$Max_disp_dist_m = "2000"
head(resp_curve_habAv_Alobs_2KmDispDist)

# 4km
gbm_mod_Alobs_4KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Alobs_4KmDispDist = plot(gbm_mod_Alobs_4KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Alobs_4KmDispDist$y = scale(resp_curve_habAv_Alobs_4KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Alobs_4KmDispDist$model = "Alobs_4KmDispDist"
resp_curve_habAv_Alobs_4KmDispDist$Species = "Alobs"
resp_curve_habAv_Alobs_4KmDispDist$Max_disp_dist_m = "1000"
head(resp_curve_habAv_Alobs_4KmDispDist)

# 6Km
gbm_mod_Alobs_6KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Alobs_6KmDispDist = plot(gbm_mod_Alobs_6KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Alobs_6KmDispDist$y = scale(resp_curve_habAv_Alobs_6KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Alobs_6KmDispDist$model = "Alobs_6KmDispDist"
resp_curve_habAv_Alobs_6KmDispDist$Species = "Alobs"
resp_curve_habAv_Alobs_6KmDispDist$Max_disp_dist_m = "6000"
head(resp_curve_habAv_Alobs_6KmDispDist)

# 8Km
gbm_mod_Alobs_8KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Alobs_8KmDispDist = plot(gbm_mod_Alobs_8KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Alobs_8KmDispDist$y = scale(resp_curve_habAv_Alobs_8KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Alobs_8KmDispDist$model = "Alobs_8KmDispDist"
resp_curve_habAv_Alobs_8KmDispDist$Species = "Alobs"
resp_curve_habAv_Alobs_8KmDispDist$Max_disp_dist_m = "8000"
head(resp_curve_habAv_Alobs_8KmDispDist)

# 10Km
gbm_mod_Alobs_10KmDispDist$var.names #query the model to get position of relevant variable
resp_curve_habAv_Alobs_10KmDispDist = plot(gbm_mod_Alobs_10KmDispDist, i.var = 7, return.grid=TRUE)
resp_curve_habAv_Alobs_10KmDispDist$y = scale(resp_curve_habAv_Alobs_10KmDispDist$y, scale = FALSE)
#Set identifier tags
resp_curve_habAv_Alobs_10KmDispDist$model = "Alobs_10KmDispDist"
resp_curve_habAv_Alobs_10KmDispDist$Species = "Alobs"
resp_curve_habAv_Alobs_10KmDispDist$Max_disp_dist_m = "10000"
head(resp_curve_habAv_Alobs_10KmDispDist)

##join to get joint response curves
resp_curve_habAv_Alobs_tab = rbind(resp_curve_habAv_Alobs_300mDispDist, resp_curve_habAv_Alobs_1KmDispDist,
                                   resp_curve_habAv_Alobs_defaultDispDist, resp_curve_habAv_Alobs_2KmDispDist, 
                                   resp_curve_habAv_Alobs_4KmDispDist, resp_curve_habAv_Alobs_6KmDispDist,
                                   resp_curve_habAv_Alobs_8KmDispDist,resp_curve_habAv_Alobs_10KmDispDist)

#Turn distances into numeric
resp_curve_habAv_Alobs_tab$Max_disp_dist_m <- as.numeric(resp_curve_habAv_Alobs_tab$Max_disp_dist_m)

#Reorder joint response curves tab
resp_curve_habAv_Alobs_tab$model <- as.factor(resp_curve_habAv_Alobs_tab$model)
resp_curve_habAv_Alobs_tab$model <-factor(resp_curve_habAv_Alobs_tab$model, 
                                          levels=c("Alobs_300mDispDist", "Alobs_1KmDispDist", "Alobs_defaultDispDist",
                                                   "Alobs_2KmDispDist", "Alobs_4KmDispDist", "Alobs_6KmDispDist", 
                                                   "Alobs_8KmDispDist","Alobs_10KmDispDist"))

#Plot
ggplot(resp_curve_habAv_Alobs_tab, aes(x=habAv, y = y, color = model, group = model)) +
  geom_point(size = 2) + geom_line(size = 0.5)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Habitat availability", y= "Occurrence-state")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))

ggplot(resp_curve_habAv_Alobs_tab, aes(x=habAv, y = y, color = Max_disp_dist_m, group = Max_disp_dist_m)) +
  geom_point() + geom_line()




#########################################################################################################
#### Generate plots relating dispersal distnace. predictive performance, components, number of trees, variable importance
# For all the different distance settings, both the species-specific ones and the uniformly-set ones ####
#########################################################################################################

#### Plot the different distance settings
setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/Network_setup")
Dispersal_Distances_table <- read.csv("Dispersal_Distances_table.csv")
Dispersal_Distances_table <-rename(Dispersal_Distances_table, "Max_DispDist" = "Max_DispDist..m." )
Dispersal_Distances_table$Distance_Setting <- as.factor(Dispersal_Distances_table$Distance_Setting)
Dispersal_Distances_table$Distance_Setting <-factor(Dispersal_Distances_table$Distance_Setting, levels=c("300m", "1Km", "Alobs", "Peagg", "Perid", "2Km", "Hyarb", "Bovar", "4Km", "Epcal", "6Km", "8Km","10Km"))

plot(Dispersal_Distances_table$Max_DispDist~Dispersal_Distances_table$Distance_Setting, xlab= " Distance Setting", ylab= "Max. Dispersal Distance", cex=5)



################################################################################
####### Get plots with all the mean AUC_cv's ###################################
#### To get plot AUC-cv ~ Max. disp. dist. & related ####

setwd('C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/Mega_df_allmodels')

#Export output dataframes of BRT's for every setting/species
head(Hyarb_DefaultDispDist_AUC_cv_tab)
write.csv(Hyarb_DefaultDispDist_AUC_cv_tab, "Hyarb_DefaultDispDist_AUC_cv_tab.csv")


################################################################################
#Do invidual df's with only the mean values of each setting and species ########

### Hyarb ######################################################################

Hyarb_defaultDispDist_meanAUCcv_tab = data.frame(mean(Hyarb_DefaultDispDist_AUC_cv_tab$value))
Hyarb_defaultDispDist_meanAUCcv_tab <- setnames(Hyarb_defaultDispDist_meanAUCcv_tab, "mean.Hyarb_DefaultDispDist_AUC_cv_tab.value.",
                                                "Mean_AUC_cv")
#Set identifier labels
Hyarb_defaultDispDist_meanAUCcv_tab$model = "Hyarb_defaultDispDist"
Hyarb_defaultDispDist_meanAUCcv_tab$Species = "Hyarb"
Hyarb_defaultDispDist_meanAUCcv_tab$Max_disp_dist_m = "2657"

Hyarb_300mDispDist_meanAUCcv_tab = data.frame(mean(Hyarb_300mDispDist_AUC_cv_tab$value))
Hyarb_300mDispDist_meanAUCcv_tab <- setnames(Hyarb_300mDispDist_meanAUCcv_tab, "mean.Hyarb_300mDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Hyarb_300mDispDist_meanAUCcv_tab$model = "Hyarb_300mDispDist"
Hyarb_300mDispDist_meanAUCcv_tab$Species = "Hyarb"
Hyarb_300mDispDist_meanAUCcv_tab$Max_disp_dist_m = "300"

Hyarb_1KmDispDist_meanAUCcv_tab = data.frame(mean(Hyarb_1KmDispDist_AUC_cv_tab$value))
Hyarb_1KmDispDist_meanAUCcv_tab <- setnames(Hyarb_1KmDispDist_meanAUCcv_tab, "mean.Hyarb_1KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Hyarb_1KmDispDist_meanAUCcv_tab$model = "Hyarb_1KmDispDist"
Hyarb_1KmDispDist_meanAUCcv_tab$Species = "Hyarb"
Hyarb_1KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "1000"

Hyarb_2KmDispDist_meanAUCcv_tab = data.frame(mean(Hyarb_2KmDispDist_AUC_cv_tab$value))
Hyarb_2KmDispDist_meanAUCcv_tab <- setnames(Hyarb_2KmDispDist_meanAUCcv_tab, "mean.Hyarb_2KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Hyarb_2KmDispDist_meanAUCcv_tab$model = "Hyarb_2KmDispDist"
Hyarb_2KmDispDist_meanAUCcv_tab$Species = "Hyarb"
Hyarb_2KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "2000"

Hyarb_4KmDispDist_meanAUCcv_tab = data.frame(mean(Hyarb_4KmDispDist_AUC_cv_tab$value))
Hyarb_4KmDispDist_meanAUCcv_tab <- setnames(Hyarb_4KmDispDist_meanAUCcv_tab, "mean.Hyarb_4KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Hyarb_4KmDispDist_meanAUCcv_tab$model = "Hyarb_4KmDispDist"
Hyarb_4KmDispDist_meanAUCcv_tab$Species = "Hyarb"
Hyarb_4KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "4000"

Hyarb_6KmDispDist_meanAUCcv_tab = data.frame(mean(Hyarb_6KmDispDist_AUC_cv_tab$value))
Hyarb_6KmDispDist_meanAUCcv_tab <- setnames(Hyarb_6KmDispDist_meanAUCcv_tab, "mean.Hyarb_6KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Hyarb_6KmDispDist_meanAUCcv_tab$model = "Hyarb_6KmDispDist"
Hyarb_6KmDispDist_meanAUCcv_tab$Species = "Hyarb"
Hyarb_6KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "6000"

Hyarb_8KmDispDist_meanAUCcv_tab = data.frame(mean(Hyarb_8KmDispDist_AUC_cv_tab$value))
Hyarb_8KmDispDist_meanAUCcv_tab <- setnames(Hyarb_8KmDispDist_meanAUCcv_tab, "mean.Hyarb_8KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Hyarb_8KmDispDist_meanAUCcv_tab$model = "Hyarb_8KmDispDist"
Hyarb_8KmDispDist_meanAUCcv_tab$Species = "Hyarb"
Hyarb_8KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "8000"

Hyarb_10KmDispDist_meanAUCcv_tab = data.frame(mean(Hyarb_10KmDispDist_AUC_cv_tab$value))
Hyarb_10KmDispDist_meanAUCcv_tab <- setnames(Hyarb_10KmDispDist_meanAUCcv_tab, "mean.Hyarb_10KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Hyarb_10KmDispDist_meanAUCcv_tab$model = "Hyarb_10KmDispDist"
Hyarb_10KmDispDist_meanAUCcv_tab$Species = "Hyarb"
Hyarb_10KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "10000"

##join them
Hyarb_meanAUCcv_tab = rbind(Hyarb_300mDispDist_meanAUCcv_tab, Hyarb_1KmDispDist_meanAUCcv_tab, Hyarb_2KmDispDist_meanAUCcv_tab,
                            Hyarb_defaultDispDist_meanAUCcv_tab, Hyarb_4KmDispDist_meanAUCcv_tab, Hyarb_6KmDispDist_meanAUCcv_tab,
                            Hyarb_8KmDispDist_meanAUCcv_tab, Hyarb_10KmDispDist_meanAUCcv_tab)


### Alobs ######################################################################

Alobs_defaultDispDist_meanAUCcv_tab = data.frame(mean(Alobs_DefaultDispDist_AUC_cv_tab$value))
Alobs_defaultDispDist_meanAUCcv_tab <- setnames(Alobs_defaultDispDist_meanAUCcv_tab, "mean.Alobs_DefaultDispDist_AUC_cv_tab.value.",
                                                "Mean_AUC_cv")
#Set identifier labels
Alobs_defaultDispDist_meanAUCcv_tab$model = "Alobs_defaultDispDist"
Alobs_defaultDispDist_meanAUCcv_tab$Species = "Alobs"
Alobs_defaultDispDist_meanAUCcv_tab$Max_disp_dist_m = "1500"

Alobs_300mDispDist_meanAUCcv_tab = data.frame(mean(Alobs_300mDispDist_AUC_cv_tab$value))
Alobs_300mDispDist_meanAUCcv_tab <- setnames(Alobs_300mDispDist_meanAUCcv_tab, "mean.Alobs_300mDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Alobs_300mDispDist_meanAUCcv_tab$model = "Alobs_300mDispDist"
Alobs_300mDispDist_meanAUCcv_tab$Species = "Alobs"
Alobs_300mDispDist_meanAUCcv_tab$Max_disp_dist_m = "300"

Alobs_1KmDispDist_meanAUCcv_tab = data.frame(mean(Alobs_1KmDispDist_AUC_cv_tab$value))
Alobs_1KmDispDist_meanAUCcv_tab <- setnames(Alobs_1KmDispDist_meanAUCcv_tab, "mean.Alobs_1KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Alobs_1KmDispDist_meanAUCcv_tab$model = "Alobs_1KmDispDist"
Alobs_1KmDispDist_meanAUCcv_tab$Species = "Alobs"
Alobs_1KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "1000"

Alobs_2KmDispDist_meanAUCcv_tab = data.frame(mean(Alobs_2KmDispDist_AUC_cv_tab$value))
Alobs_2KmDispDist_meanAUCcv_tab <- setnames(Alobs_2KmDispDist_meanAUCcv_tab, "mean.Alobs_2KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Alobs_2KmDispDist_meanAUCcv_tab$model = "Alobs_2KmDispDist"
Alobs_2KmDispDist_meanAUCcv_tab$Species = "Alobs"
Alobs_2KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "2000"

Alobs_4KmDispDist_meanAUCcv_tab = data.frame(mean(Alobs_4KmDispDist_AUC_cv_tab$value))
Alobs_4KmDispDist_meanAUCcv_tab <- setnames(Alobs_4KmDispDist_meanAUCcv_tab, "mean.Alobs_4KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Alobs_4KmDispDist_meanAUCcv_tab$model = "Alobs_4KmDispDist"
Alobs_4KmDispDist_meanAUCcv_tab$Species = "Alobs"
Alobs_4KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "4000"

Alobs_6KmDispDist_meanAUCcv_tab = data.frame(mean(Alobs_6KmDispDist_AUC_cv_tab$value))
Alobs_6KmDispDist_meanAUCcv_tab <- setnames(Alobs_6KmDispDist_meanAUCcv_tab, "mean.Alobs_6KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Alobs_6KmDispDist_meanAUCcv_tab$model = "Alobs_6KmDispDist"
Alobs_6KmDispDist_meanAUCcv_tab$Species = "Alobs"
Alobs_6KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "6000"

Alobs_8KmDispDist_meanAUCcv_tab = data.frame(mean(Alobs_8KmDispDist_AUC_cv_tab$value))
Alobs_8KmDispDist_meanAUCcv_tab <- setnames(Alobs_8KmDispDist_meanAUCcv_tab, "mean.Alobs_8KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Alobs_8KmDispDist_meanAUCcv_tab$model = "Alobs_8KmDispDist"
Alobs_8KmDispDist_meanAUCcv_tab$Species = "Alobs"
Alobs_8KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "8000"

Alobs_10KmDispDist_meanAUCcv_tab = data.frame(mean(Alobs_10KmDispDist_AUC_cv_tab$value))
Alobs_10KmDispDist_meanAUCcv_tab <- setnames(Alobs_10KmDispDist_meanAUCcv_tab, "mean.Alobs_10KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Alobs_10KmDispDist_meanAUCcv_tab$model = "Alobs_10KmDispDist"
Alobs_10KmDispDist_meanAUCcv_tab$Species = "Alobs"
Alobs_10KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "10000"

##join them
Alobs_meanAUCcv_tab = rbind(Alobs_300mDispDist_meanAUCcv_tab, Alobs_1KmDispDist_meanAUCcv_tab, Alobs_defaultDispDist_meanAUCcv_tab,
                            Alobs_2KmDispDist_meanAUCcv_tab, Alobs_4KmDispDist_meanAUCcv_tab, Alobs_6KmDispDist_meanAUCcv_tab,
                            Alobs_8KmDispDist_meanAUCcv_tab, Alobs_10KmDispDist_meanAUCcv_tab)


### Bovar ######################################################################

Bovar_defaultDispDist_meanAUCcv_tab = data.frame(mean(Bovar_DefaultDispDist_AUC_cv_tab$value))
Bovar_defaultDispDist_meanAUCcv_tab <- setnames(Bovar_defaultDispDist_meanAUCcv_tab, "mean.Bovar_DefaultDispDist_AUC_cv_tab.value.",
                                                "Mean_AUC_cv")
#Set identifier labels
Bovar_defaultDispDist_meanAUCcv_tab$model = "Bovar_defaultDispDist"
Bovar_defaultDispDist_meanAUCcv_tab$Species = "Bovar"
Bovar_defaultDispDist_meanAUCcv_tab$Max_disp_dist_m = "4000"

Bovar_300mDispDist_meanAUCcv_tab = data.frame(mean(Bovar_300mDispDist_AUC_cv_tab$value))
Bovar_300mDispDist_meanAUCcv_tab <- setnames(Bovar_300mDispDist_meanAUCcv_tab, "mean.Bovar_300mDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Bovar_300mDispDist_meanAUCcv_tab$model = "Bovar_300mDispDist"
Bovar_300mDispDist_meanAUCcv_tab$Species = "Bovar"
Bovar_300mDispDist_meanAUCcv_tab$Max_disp_dist_m = "300"

Bovar_1KmDispDist_meanAUCcv_tab = data.frame(mean(Bovar_1KmDispDist_AUC_cv_tab$value))
Bovar_1KmDispDist_meanAUCcv_tab <- setnames(Bovar_1KmDispDist_meanAUCcv_tab, "mean.Bovar_1KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Bovar_1KmDispDist_meanAUCcv_tab$model = "Bovar_1KmDispDist"
Bovar_1KmDispDist_meanAUCcv_tab$Species = "Bovar"
Bovar_1KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "1000"

Bovar_2KmDispDist_meanAUCcv_tab = data.frame(mean(Bovar_2KmDispDist_AUC_cv_tab$value))
Bovar_2KmDispDist_meanAUCcv_tab <- setnames(Bovar_2KmDispDist_meanAUCcv_tab, "mean.Bovar_2KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Bovar_2KmDispDist_meanAUCcv_tab$model = "Bovar_2KmDispDist"
Bovar_2KmDispDist_meanAUCcv_tab$Species = "Bovar"
Bovar_2KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "2000"

#4 Km = default for this species
# Bovar_4KmDispDist_meanAUCcv_tab = data.frame(mean(Bovar_4KmDispDist_AUC_cv_tab$value))
# Bovar_4KmDispDist_meanAUCcv_tab <- setnames(Bovar_4KmDispDist_meanAUCcv_tab, "mean.Bovar_4KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
# #Set identifier labels
# Bovar_4KmDispDist_meanAUCcv_tab$model = "Bovar_4KmDispDist"
# Bovar_4KmDispDist_meanAUCcv_tab$Species = "Bovar"
# Bovar_4KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "4000"

Bovar_6KmDispDist_meanAUCcv_tab = data.frame(mean(Bovar_6KmDispDist_AUC_cv_tab$value))
Bovar_6KmDispDist_meanAUCcv_tab <- setnames(Bovar_6KmDispDist_meanAUCcv_tab, "mean.Bovar_6KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Bovar_6KmDispDist_meanAUCcv_tab$model = "Bovar_6KmDispDist"
Bovar_6KmDispDist_meanAUCcv_tab$Species = "Bovar"
Bovar_6KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "6000"

Bovar_8KmDispDist_meanAUCcv_tab = data.frame(mean(Bovar_8KmDispDist_AUC_cv_tab$value))
Bovar_8KmDispDist_meanAUCcv_tab <- setnames(Bovar_8KmDispDist_meanAUCcv_tab, "mean.Bovar_8KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Bovar_8KmDispDist_meanAUCcv_tab$model = "Bovar_8KmDispDist"
Bovar_8KmDispDist_meanAUCcv_tab$Species = "Bovar"
Bovar_8KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "8000"

Bovar_10KmDispDist_meanAUCcv_tab = data.frame(mean(Bovar_10KmDispDist_AUC_cv_tab$value))
Bovar_10KmDispDist_meanAUCcv_tab <- setnames(Bovar_10KmDispDist_meanAUCcv_tab, "mean.Bovar_10KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Bovar_10KmDispDist_meanAUCcv_tab$model = "Bovar_10KmDispDist"
Bovar_10KmDispDist_meanAUCcv_tab$Species = "Bovar"
Bovar_10KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "10000"

##join them
Bovar_meanAUCcv_tab = rbind(Bovar_300mDispDist_meanAUCcv_tab, Bovar_1KmDispDist_meanAUCcv_tab, 
                            Bovar_2KmDispDist_meanAUCcv_tab, Bovar_defaultDispDist_meanAUCcv_tab, Bovar_6KmDispDist_meanAUCcv_tab,
                            Bovar_8KmDispDist_meanAUCcv_tab, Bovar_10KmDispDist_meanAUCcv_tab)


### Epcal ######################################################################

Epcal_defaultDispDist_meanAUCcv_tab = data.frame(mean(Epcal_DefaultDispDist_AUC_cv_tab$value))
Epcal_defaultDispDist_meanAUCcv_tab <- setnames(Epcal_defaultDispDist_meanAUCcv_tab, "mean.Epcal_DefaultDispDist_AUC_cv_tab.value.",
                                                "Mean_AUC_cv")
#Set identifier labels
Epcal_defaultDispDist_meanAUCcv_tab$model = "Epcal_defaultDispDist"
Epcal_defaultDispDist_meanAUCcv_tab$Species = "Epcal"
Epcal_defaultDispDist_meanAUCcv_tab$Max_disp_dist_m = "4411"

Epcal_300mDispDist_meanAUCcv_tab = data.frame(mean(Epcal_300mDispDist_AUC_cv_tab$value))
Epcal_300mDispDist_meanAUCcv_tab <- setnames(Epcal_300mDispDist_meanAUCcv_tab, "mean.Epcal_300mDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Epcal_300mDispDist_meanAUCcv_tab$model = "Epcal_300mDispDist"
Epcal_300mDispDist_meanAUCcv_tab$Species = "Epcal"
Epcal_300mDispDist_meanAUCcv_tab$Max_disp_dist_m = "300"

Epcal_1KmDispDist_meanAUCcv_tab = data.frame(mean(Epcal_1KmDispDist_AUC_cv_tab$value))
Epcal_1KmDispDist_meanAUCcv_tab <- setnames(Epcal_1KmDispDist_meanAUCcv_tab, "mean.Epcal_1KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Epcal_1KmDispDist_meanAUCcv_tab$model = "Epcal_1KmDispDist"
Epcal_1KmDispDist_meanAUCcv_tab$Species = "Epcal"
Epcal_1KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "1000"

Epcal_2KmDispDist_meanAUCcv_tab = data.frame(mean(Epcal_2KmDispDist_AUC_cv_tab$value))
Epcal_2KmDispDist_meanAUCcv_tab <- setnames(Epcal_2KmDispDist_meanAUCcv_tab, "mean.Epcal_2KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Epcal_2KmDispDist_meanAUCcv_tab$model = "Epcal_2KmDispDist"
Epcal_2KmDispDist_meanAUCcv_tab$Species = "Epcal"
Epcal_2KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "2000"

Epcal_4KmDispDist_meanAUCcv_tab = data.frame(mean(Epcal_4KmDispDist_AUC_cv_tab$value))
Epcal_4KmDispDist_meanAUCcv_tab <- setnames(Epcal_4KmDispDist_meanAUCcv_tab, "mean.Epcal_4KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Epcal_4KmDispDist_meanAUCcv_tab$model = "Epcal_4KmDispDist"
Epcal_4KmDispDist_meanAUCcv_tab$Species = "Epcal"
Epcal_4KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "4000"

Epcal_6KmDispDist_meanAUCcv_tab = data.frame(mean(Epcal_6KmDispDist_AUC_cv_tab$value))
Epcal_6KmDispDist_meanAUCcv_tab <- setnames(Epcal_6KmDispDist_meanAUCcv_tab, "mean.Epcal_6KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Epcal_6KmDispDist_meanAUCcv_tab$model = "Epcal_6KmDispDist"
Epcal_6KmDispDist_meanAUCcv_tab$Species = "Epcal"
Epcal_6KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "6000"

Epcal_8KmDispDist_meanAUCcv_tab = data.frame(mean(Epcal_8KmDispDist_AUC_cv_tab$value))
Epcal_8KmDispDist_meanAUCcv_tab <- setnames(Epcal_8KmDispDist_meanAUCcv_tab, "mean.Epcal_8KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Epcal_8KmDispDist_meanAUCcv_tab$model = "Epcal_8KmDispDist"
Epcal_8KmDispDist_meanAUCcv_tab$Species = "Epcal"
Epcal_8KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "8000"

Epcal_10KmDispDist_meanAUCcv_tab = data.frame(mean(Epcal_10KmDispDist_AUC_cv_tab$value))
Epcal_10KmDispDist_meanAUCcv_tab <- setnames(Epcal_10KmDispDist_meanAUCcv_tab, "mean.Epcal_10KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Epcal_10KmDispDist_meanAUCcv_tab$model = "Epcal_10KmDispDist"
Epcal_10KmDispDist_meanAUCcv_tab$Species = "Epcal"
Epcal_10KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "10000"

##join them
Epcal_meanAUCcv_tab = rbind(Epcal_300mDispDist_meanAUCcv_tab, Epcal_1KmDispDist_meanAUCcv_tab, Epcal_2KmDispDist_meanAUCcv_tab,
                            Epcal_4KmDispDist_meanAUCcv_tab, Epcal_defaultDispDist_meanAUCcv_tab,Epcal_6KmDispDist_meanAUCcv_tab,
                            Epcal_8KmDispDist_meanAUCcv_tab, Epcal_10KmDispDist_meanAUCcv_tab)


### Peagg ######################################################################

Peagg_defaultDispDist_meanAUCcv_tab = data.frame(mean(Peagg_DefaultDispDist_AUC_cv_tab$value))
Peagg_defaultDispDist_meanAUCcv_tab <- setnames(Peagg_defaultDispDist_meanAUCcv_tab, "mean.Peagg_DefaultDispDist_AUC_cv_tab.value.",
                                                "Mean_AUC_cv")
#Set identifier labels
Peagg_defaultDispDist_meanAUCcv_tab$model = "Peagg_defaultDispDist"
Peagg_defaultDispDist_meanAUCcv_tab$Species = "Peagg"
Peagg_defaultDispDist_meanAUCcv_tab$Max_disp_dist_m = "1760"

Peagg_300mDispDist_meanAUCcv_tab = data.frame(mean(Peagg_300mDispDist_AUC_cv_tab$value))
Peagg_300mDispDist_meanAUCcv_tab <- setnames(Peagg_300mDispDist_meanAUCcv_tab, "mean.Peagg_300mDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Peagg_300mDispDist_meanAUCcv_tab$model = "Peagg_300mDispDist"
Peagg_300mDispDist_meanAUCcv_tab$Species = "Peagg"
Peagg_300mDispDist_meanAUCcv_tab$Max_disp_dist_m = "300"

Peagg_1KmDispDist_meanAUCcv_tab = data.frame(mean(Peagg_1KmDispDist_AUC_cv_tab$value))
Peagg_1KmDispDist_meanAUCcv_tab <- setnames(Peagg_1KmDispDist_meanAUCcv_tab, "mean.Peagg_1KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Peagg_1KmDispDist_meanAUCcv_tab$model = "Peagg_1KmDispDist"
Peagg_1KmDispDist_meanAUCcv_tab$Species = "Peagg"
Peagg_1KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "1000"

Peagg_2KmDispDist_meanAUCcv_tab = data.frame(mean(Peagg_2KmDispDist_AUC_cv_tab$value))
Peagg_2KmDispDist_meanAUCcv_tab <- setnames(Peagg_2KmDispDist_meanAUCcv_tab, "mean.Peagg_2KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Peagg_2KmDispDist_meanAUCcv_tab$model = "Peagg_2KmDispDist"
Peagg_2KmDispDist_meanAUCcv_tab$Species = "Peagg"
Peagg_2KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "2000"

Peagg_4KmDispDist_meanAUCcv_tab = data.frame(mean(Peagg_4KmDispDist_AUC_cv_tab$value))
Peagg_4KmDispDist_meanAUCcv_tab <- setnames(Peagg_4KmDispDist_meanAUCcv_tab, "mean.Peagg_4KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Peagg_4KmDispDist_meanAUCcv_tab$model = "Peagg_4KmDispDist"
Peagg_4KmDispDist_meanAUCcv_tab$Species = "Peagg"
Peagg_4KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "4000"

Peagg_6KmDispDist_meanAUCcv_tab = data.frame(mean(Peagg_6KmDispDist_AUC_cv_tab$value))
Peagg_6KmDispDist_meanAUCcv_tab <- setnames(Peagg_6KmDispDist_meanAUCcv_tab, "mean.Peagg_6KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Peagg_6KmDispDist_meanAUCcv_tab$model = "Peagg_6KmDispDist"
Peagg_6KmDispDist_meanAUCcv_tab$Species = "Peagg"
Peagg_6KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "6000"

Peagg_8KmDispDist_meanAUCcv_tab = data.frame(mean(Peagg_8KmDispDist_AUC_cv_tab$value))
Peagg_8KmDispDist_meanAUCcv_tab <- setnames(Peagg_8KmDispDist_meanAUCcv_tab, "mean.Peagg_8KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Peagg_8KmDispDist_meanAUCcv_tab$model = "Peagg_8KmDispDist"
Peagg_8KmDispDist_meanAUCcv_tab$Species = "Peagg"
Peagg_8KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "8000"

Peagg_10KmDispDist_meanAUCcv_tab = data.frame(mean(Peagg_10KmDispDist_AUC_cv_tab$value))
Peagg_10KmDispDist_meanAUCcv_tab <- setnames(Peagg_10KmDispDist_meanAUCcv_tab, "mean.Peagg_10KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Peagg_10KmDispDist_meanAUCcv_tab$model = "Peagg_10KmDispDist"
Peagg_10KmDispDist_meanAUCcv_tab$Species = "Peagg"
Peagg_10KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "10000"

##join them
Peagg_meanAUCcv_tab = rbind(Peagg_300mDispDist_meanAUCcv_tab, Peagg_1KmDispDist_meanAUCcv_tab, Peagg_defaultDispDist_meanAUCcv_tab,
                            Peagg_2KmDispDist_meanAUCcv_tab, Peagg_4KmDispDist_meanAUCcv_tab, Peagg_6KmDispDist_meanAUCcv_tab,
                            Peagg_8KmDispDist_meanAUCcv_tab, Peagg_10KmDispDist_meanAUCcv_tab)


### Perid ######################################################################

Perid_defaultDispDist_meanAUCcv_tab = data.frame(mean(Perid_DefaultDispDist_AUC_cv_tab$value))
Perid_defaultDispDist_meanAUCcv_tab <- setnames(Perid_defaultDispDist_meanAUCcv_tab, "mean.Perid_DefaultDispDist_AUC_cv_tab.value.",
                                                "Mean_AUC_cv")
#Set identifier labels
Perid_defaultDispDist_meanAUCcv_tab$model = "Perid_defaultDispDist"
Perid_defaultDispDist_meanAUCcv_tab$Species = "Perid"
Perid_defaultDispDist_meanAUCcv_tab$Max_disp_dist_m = "1760"

Perid_300mDispDist_meanAUCcv_tab = data.frame(mean(Perid_300mDispDist_AUC_cv_tab$value))
Perid_300mDispDist_meanAUCcv_tab <- setnames(Perid_300mDispDist_meanAUCcv_tab, "mean.Perid_300mDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Perid_300mDispDist_meanAUCcv_tab$model = "Perid_300mDispDist"
Perid_300mDispDist_meanAUCcv_tab$Species = "Perid"
Perid_300mDispDist_meanAUCcv_tab$Max_disp_dist_m = "300"

Perid_1KmDispDist_meanAUCcv_tab = data.frame(mean(Perid_1KmDispDist_AUC_cv_tab$value))
Perid_1KmDispDist_meanAUCcv_tab <- setnames(Perid_1KmDispDist_meanAUCcv_tab, "mean.Perid_1KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Perid_1KmDispDist_meanAUCcv_tab$model = "Perid_1KmDispDist"
Perid_1KmDispDist_meanAUCcv_tab$Species = "Perid"
Perid_1KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "1000"

Perid_2KmDispDist_meanAUCcv_tab = data.frame(mean(Perid_2KmDispDist_AUC_cv_tab$value))
Perid_2KmDispDist_meanAUCcv_tab <- setnames(Perid_2KmDispDist_meanAUCcv_tab, "mean.Perid_2KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Perid_2KmDispDist_meanAUCcv_tab$model = "Perid_2KmDispDist"
Perid_2KmDispDist_meanAUCcv_tab$Species = "Perid"
Perid_2KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "2000"

Perid_4KmDispDist_meanAUCcv_tab = data.frame(mean(Perid_4KmDispDist_AUC_cv_tab$value))
Perid_4KmDispDist_meanAUCcv_tab <- setnames(Perid_4KmDispDist_meanAUCcv_tab, "mean.Perid_4KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Perid_4KmDispDist_meanAUCcv_tab$model = "Perid_4KmDispDist"
Perid_4KmDispDist_meanAUCcv_tab$Species = "Perid"
Perid_4KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "4000"

Perid_6KmDispDist_meanAUCcv_tab = data.frame(mean(Perid_6KmDispDist_AUC_cv_tab$value))
Perid_6KmDispDist_meanAUCcv_tab <- setnames(Perid_6KmDispDist_meanAUCcv_tab, "mean.Perid_6KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Perid_6KmDispDist_meanAUCcv_tab$model = "Perid_6KmDispDist"
Perid_6KmDispDist_meanAUCcv_tab$Species = "Perid"
Perid_6KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "6000"

Perid_8KmDispDist_meanAUCcv_tab = data.frame(mean(Perid_8KmDispDist_AUC_cv_tab$value))
Perid_8KmDispDist_meanAUCcv_tab <- setnames(Perid_8KmDispDist_meanAUCcv_tab, "mean.Perid_8KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Perid_8KmDispDist_meanAUCcv_tab$model = "Perid_8KmDispDist"
Perid_8KmDispDist_meanAUCcv_tab$Species = "Perid"
Perid_8KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "8000"

Perid_10KmDispDist_meanAUCcv_tab = data.frame(mean(Perid_10KmDispDist_AUC_cv_tab$value))
Perid_10KmDispDist_meanAUCcv_tab <- setnames(Perid_10KmDispDist_meanAUCcv_tab, "mean.Perid_10KmDispDist_AUC_cv_tab.value.","Mean_AUC_cv")
#Set identifier labels
Perid_10KmDispDist_meanAUCcv_tab$model = "Perid_10KmDispDist"
Perid_10KmDispDist_meanAUCcv_tab$Species = "Perid"
Perid_10KmDispDist_meanAUCcv_tab$Max_disp_dist_m = "10000"

##join them
Perid_meanAUCcv_tab = rbind(Perid_300mDispDist_meanAUCcv_tab, Perid_1KmDispDist_meanAUCcv_tab, Perid_defaultDispDist_meanAUCcv_tab,
                            Perid_2KmDispDist_meanAUCcv_tab, Perid_4KmDispDist_meanAUCcv_tab, Perid_6KmDispDist_meanAUCcv_tab,
                            Perid_8KmDispDist_meanAUCcv_tab, Perid_10KmDispDist_meanAUCcv_tab)


################################################################################
#### Join the AUC-cv dfs of all the species into one ###########################

All_Spp_meanAUCcv_tab = rbind(Alobs_meanAUCcv_tab, Bovar_meanAUCcv_tab, Epcal_meanAUCcv_tab,
                              Hyarb_meanAUCcv_tab, Peagg_meanAUCcv_tab, Perid_meanAUCcv_tab)

#Turn distances into numeric
All_Spp_meanAUCcv_tab$Max_disp_dist_m <- as.numeric(All_Spp_meanAUCcv_tab$Max_disp_dist_m)

#Do linear model of Mean_AUC_cv~Max_disp_dist_m to generate trendline

#Plot
ggplot(All_Spp_meanAUCcv_tab, aes(x=Max_disp_dist_m, y = Mean_AUC_cv, color = Species, group = Species)) +
  geom_point(size = 3) + geom_line(size = 1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Mean AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))
# scale_y_continuous(breaks = c(0.1:1))

  #se -> std error -> 95% confidence interval
ggplot(All_Spp_meanAUCcv_tab, aes(x=Max_disp_dist_m, y = Mean_AUC_cv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Mean AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))



################################################################################
#### Plot analyses with components #############################################

### Import table with component values (created by hand from results generated here) ###
Table_edges_components <- read.csv("Table_edges_components.csv")
head(Table_edges_components)
Table_edges_components <- setnames(Table_edges_components, "X", "Max_disp_dist_m")
Table_edges_components <- setnames(Table_edges_components, "Species", "Gender_Species")
plot(Number_of_components~Max_disp_dist_m, data = Table_edges_components)


#### Join All_Spp_meanAUCcv_tab and Table_edges_components ####
#inner join= join by common variable names
All_Spp_meanAUCcv_Components <- merge(Table_edges_components, All_Spp_meanAUCcv_tab, by = "model")
head(All_Spp_meanAUCcv_Components)
All_Spp_meanAUCcv_Components$Max_disp_dist_m.y <- NULL
All_Spp_meanAUCcv_Components <- setnames(All_Spp_meanAUCcv_Components, "Max_disp_dist_m.x", "Max_disp_dist_m")

ggplot(All_Spp_meanAUCcv_Components, aes(x=Max_disp_dist_m, y = Number_of_components, color = Species, group = Species)) +
  geom_point(size = 4) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Number of components")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))
# theme(legend.key.size = unit(10,"point"))

ggplot(All_Spp_meanAUCcv_Components, aes(x= log(Number_of_components), y =Mean_AUC_cv , color = Species, group = Species)) +
  geom_smooth(method = "lm")+
  geom_point(size = 3) + #geom_line(size = 1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "log(Number of components)", y= "Mean AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

plot(All_Spp_meanAUCcv_Components)


#####################################################################
#### Plot relationship number of components vs number of patches #### 

# scale_color_manual(values = c("red", "green", "blue"))
#scale_color_manual(values = c("brown2", "darkorange2", "khaki4","greenyellow", "olivedrab","green4","turquoise", "deepskyblue2", "mediumpurple1", "magenta", "deeppink1"))+

All_Spp_meanAUCcv_Components$Max_Disp_Dist_factor <- as.factor(All_Spp_meanAUCcv_Components$Max_disp_dist_m)
levels(All_Spp_meanAUCcv_Components$Max_Disp_Dist_factor)

ggplot(All_Spp_meanAUCcv_Components, aes(x= Number_of_patches, y = Number_of_components , color = Max_Disp_Dist_factor, group = Max_Disp_Dist_factor)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) + #geom_line(size = 1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Number of patches", y= "Number of components")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))



###################################################################################################
###### Get data of importance by Variable, join to the main df "All_Spp_meanAUCcv_Components" #####
##### Do plots also involving variabble importance #### ########################

#Load the default individual variable importance tabs
#HSI
Alobs_defaultDispDist_HSI_tab <- read.csv("Alobs_defaultDispDist_HSI_tab.csv")
Bovar_defaultDispDist_HSI_tab <- read.csv("Bovar_defaultDispDist_HSI_tab.csv")
Epcal_defaultDispDist_HSI_tab <- read.csv("Epcal_defaultDispDist_HSI_tab.csv")
Hyarb_defaultDispDist_HSI_tab <- read.csv("Hyarb_defaultDispDist_HSI_tab.csv")
Peagg_defaultDispDist_HSI_tab <- read.csv("Peagg_defaultDispDist_HSI_tab.csv")
Perid_defaultDispDist_HSI_tab <- read.csv("Perid_defaultDispDist_HSI_tab.csv")

#EgoSize = "3rd. ord. neigh."
Alobs_defaultDispDist_EgoSize_tab <- read.csv("Alobs_defaultDispDist_EgoSize_tab.csv")
Bovar_defaultDispDist_EgoSize_tab <- read.csv("Bovar_defaultDispDist_EgoSize_tab.csv")
Epcal_defaultDispDist_EgoSize_tab <- read.csv("Epcal_defaultDispDist_EgoSize_tab.csv")
Hyarb_defaultDispDist_EgoSize_tab <- read.csv("Hyarb_defaultDispDist_EgoSize_tab.csv")
Peagg_defaultDispDist_EgoSize_tab <- read.csv("Peagg_defaultDispDist_EgoSize_tab.csv")
Perid_defaultDispDist_EgoSize_tab <- read.csv("Perid_defaultDispDist_EgoSize_tab.csv")

#habAv
Alobs_defaultDispDist_habAv_tab <- read.csv("Alobs_defaultDispDist_habAv_tab.csv")
Bovar_defaultDispDist_habAv_tab <- read.csv("Bovar_defaultDispDist_habAv_tab.csv")
Epcal_defaultDispDist_habAv_tab <- read.csv("Epcal_defaultDispDist_habAv_tab.csv")
Hyarb_defaultDispDist_habAv_tab <- read.csv("Hyarb_defaultDispDist_habAv_tab.csv")
Peagg_defaultDispDist_habAv_tab <- read.csv("Peagg_defaultDispDist_habAv_tab.csv")
Perid_defaultDispDist_habAv_tab <- read.csv("Perid_defaultDispDist_habAv_tab.csv")


#### Do the dataframes with the mean importance values ###### 

### HSI #########################################################################

## Hyarb

Hyarb_defaultDispDist_HSI_mean_importance = data.frame(mean(Hyarb_defaultDispDist_HSI_tab$imp))
Hyarb_defaultDispDist_HSI_mean_importance <- setnames(Hyarb_defaultDispDist_HSI_mean_importance, "mean.Hyarb_defaultDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Hyarb_defaultDispDist_HSI_mean_importance$model = "Hyarb_defaultDispDist"

Hyarb_300mDispDist_HSI_mean_importance = data.frame(mean(Hyarb_300mDispDist_HSI_tab$imp))
Hyarb_300mDispDist_HSI_mean_importance <- setnames(Hyarb_300mDispDist_HSI_mean_importance, "mean.Hyarb_300mDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Hyarb_300mDispDist_HSI_mean_importance$model = "Hyarb_300mDispDist"

Hyarb_1KmDispDist_HSI_mean_importance = data.frame(mean(Hyarb_1KmDispDist_HSI_tab$imp))
Hyarb_1KmDispDist_HSI_mean_importance <- setnames(Hyarb_1KmDispDist_HSI_mean_importance, "mean.Hyarb_1KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Hyarb_1KmDispDist_HSI_mean_importance$model = "Hyarb_1KmDispDist"

Hyarb_2KmDispDist_HSI_mean_importance = data.frame(mean(Hyarb_2KmDispDist_HSI_tab$imp))
Hyarb_2KmDispDist_HSI_mean_importance <- setnames(Hyarb_2KmDispDist_HSI_mean_importance, "mean.Hyarb_2KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Hyarb_2KmDispDist_HSI_mean_importance$model = "Hyarb_2KmDispDist"

Hyarb_4KmDispDist_HSI_mean_importance = data.frame(mean(Hyarb_4KmDispDist_HSI_tab$imp))
Hyarb_4KmDispDist_HSI_mean_importance <- setnames(Hyarb_4KmDispDist_HSI_mean_importance, "mean.Hyarb_4KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Hyarb_4KmDispDist_HSI_mean_importance$model = "Hyarb_4KmDispDist"

Hyarb_6KmDispDist_HSI_mean_importance = data.frame(mean(Hyarb_6KmDispDist_HSI_tab$imp))
Hyarb_6KmDispDist_HSI_mean_importance <- setnames(Hyarb_6KmDispDist_HSI_mean_importance, "mean.Hyarb_6KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Hyarb_6KmDispDist_HSI_mean_importance$model = "Hyarb_6KmDispDist"

Hyarb_8KmDispDist_HSI_mean_importance = data.frame(mean(Hyarb_8KmDispDist_HSI_tab$imp))
Hyarb_8KmDispDist_HSI_mean_importance <- setnames(Hyarb_8KmDispDist_HSI_mean_importance, "mean.Hyarb_8KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Hyarb_8KmDispDist_HSI_mean_importance$model = "Hyarb_8KmDispDist"

Hyarb_10KmDispDist_HSI_mean_importance = data.frame(mean(Hyarb_10KmDispDist_HSI_tab$imp))
Hyarb_10KmDispDist_HSI_mean_importance <- setnames(Hyarb_10KmDispDist_HSI_mean_importance, "mean.Hyarb_10KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Hyarb_10KmDispDist_HSI_mean_importance$model = "Hyarb_10KmDispDist"

##join them
Hyarb_HSI_mean_importance_tab = rbind(Hyarb_300mDispDist_HSI_mean_importance, Hyarb_1KmDispDist_HSI_mean_importance, Hyarb_defaultDispDist_HSI_mean_importance,
                                      Hyarb_2KmDispDist_HSI_mean_importance, Hyarb_4KmDispDist_HSI_mean_importance, Hyarb_6KmDispDist_HSI_mean_importance,
                                      Hyarb_8KmDispDist_HSI_mean_importance, Hyarb_10KmDispDist_HSI_mean_importance)


## Alobs

Alobs_defaultDispDist_HSI_mean_importance = data.frame(mean(Alobs_defaultDispDist_HSI_tab$imp))
Alobs_defaultDispDist_HSI_mean_importance <- setnames(Alobs_defaultDispDist_HSI_mean_importance, "mean.Alobs_defaultDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Alobs_defaultDispDist_HSI_mean_importance$model = "Alobs_defaultDispDist"

Alobs_300mDispDist_HSI_mean_importance = data.frame(mean(Alobs_300mDispDist_HSI_tab$imp))
Alobs_300mDispDist_HSI_mean_importance <- setnames(Alobs_300mDispDist_HSI_mean_importance, "mean.Alobs_300mDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Alobs_300mDispDist_HSI_mean_importance$model = "Alobs_300mDispDist"

Alobs_1KmDispDist_HSI_mean_importance = data.frame(mean(Alobs_1KmDispDist_HSI_tab$imp))
Alobs_1KmDispDist_HSI_mean_importance <- setnames(Alobs_1KmDispDist_HSI_mean_importance, "mean.Alobs_1KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Alobs_1KmDispDist_HSI_mean_importance$model = "Alobs_1KmDispDist"

Alobs_2KmDispDist_HSI_mean_importance = data.frame(mean(Alobs_2KmDispDist_HSI_tab$imp))
Alobs_2KmDispDist_HSI_mean_importance <- setnames(Alobs_2KmDispDist_HSI_mean_importance, "mean.Alobs_2KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Alobs_2KmDispDist_HSI_mean_importance$model = "Alobs_2KmDispDist"

Alobs_4KmDispDist_HSI_mean_importance = data.frame(mean(Alobs_4KmDispDist_HSI_tab$imp))
Alobs_4KmDispDist_HSI_mean_importance <- setnames(Alobs_4KmDispDist_HSI_mean_importance, "mean.Alobs_4KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Alobs_4KmDispDist_HSI_mean_importance$model = "Alobs_4KmDispDist"

Alobs_6KmDispDist_HSI_mean_importance = data.frame(mean(Alobs_6KmDispDist_HSI_tab$imp))
Alobs_6KmDispDist_HSI_mean_importance <- setnames(Alobs_6KmDispDist_HSI_mean_importance, "mean.Alobs_6KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Alobs_6KmDispDist_HSI_mean_importance$model = "Alobs_6KmDispDist"

Alobs_8KmDispDist_HSI_mean_importance = data.frame(mean(Alobs_8KmDispDist_HSI_tab$imp))
Alobs_8KmDispDist_HSI_mean_importance <- setnames(Alobs_8KmDispDist_HSI_mean_importance, "mean.Alobs_8KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Alobs_8KmDispDist_HSI_mean_importance$model = "Alobs_8KmDispDist"

Alobs_10KmDispDist_HSI_mean_importance = data.frame(mean(Alobs_10KmDispDist_HSI_tab$imp))
Alobs_10KmDispDist_HSI_mean_importance <- setnames(Alobs_10KmDispDist_HSI_mean_importance, "mean.Alobs_10KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Alobs_10KmDispDist_HSI_mean_importance$model = "Alobs_10KmDispDist"

##join them
Alobs_HSI_mean_importance_tab = rbind(Alobs_300mDispDist_HSI_mean_importance, Alobs_1KmDispDist_HSI_mean_importance, Alobs_defaultDispDist_HSI_mean_importance,
                                      Alobs_2KmDispDist_HSI_mean_importance, Alobs_4KmDispDist_HSI_mean_importance, Alobs_6KmDispDist_HSI_mean_importance,
                                      Alobs_8KmDispDist_HSI_mean_importance, Alobs_10KmDispDist_HSI_mean_importance)


## Bovar

Bovar_defaultDispDist_HSI_mean_importance = data.frame(mean(Bovar_defaultDispDist_HSI_tab$imp))
Bovar_defaultDispDist_HSI_mean_importance <- setnames(Bovar_defaultDispDist_HSI_mean_importance, "mean.Bovar_defaultDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Bovar_defaultDispDist_HSI_mean_importance$model = "Bovar_defaultDispDist"

Bovar_300mDispDist_HSI_mean_importance = data.frame(mean(Bovar_300mDispDist_HSI_tab$imp))
Bovar_300mDispDist_HSI_mean_importance <- setnames(Bovar_300mDispDist_HSI_mean_importance, "mean.Bovar_300mDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Bovar_300mDispDist_HSI_mean_importance$model = "Bovar_300mDispDist"

Bovar_1KmDispDist_HSI_mean_importance = data.frame(mean(Bovar_1KmDispDist_HSI_tab$imp))
Bovar_1KmDispDist_HSI_mean_importance <- setnames(Bovar_1KmDispDist_HSI_mean_importance, "mean.Bovar_1KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Bovar_1KmDispDist_HSI_mean_importance$model = "Bovar_1KmDispDist"

Bovar_2KmDispDist_HSI_mean_importance = data.frame(mean(Bovar_2KmDispDist_HSI_tab$imp))
Bovar_2KmDispDist_HSI_mean_importance <- setnames(Bovar_2KmDispDist_HSI_mean_importance, "mean.Bovar_2KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Bovar_2KmDispDist_HSI_mean_importance$model = "Bovar_2KmDispDist"

#Default for Bovar = 4Km
# Bovar_4KmDispDist_HSI_mean_importance = data.frame(mean(Bovar_4KmDispDist_HSI_tab$imp))
# Bovar_4KmDispDist_HSI_mean_importance <- setnames(Bovar_4KmDispDist_HSI_mean_importance, "mean.Bovar_4KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
# #Set identifier labels
# Bovar_4KmDispDist_HSI_mean_importance$model = "Bovar_4KmDispDist"

Bovar_6KmDispDist_HSI_mean_importance = data.frame(mean(Bovar_6KmDispDist_HSI_tab$imp))
Bovar_6KmDispDist_HSI_mean_importance <- setnames(Bovar_6KmDispDist_HSI_mean_importance, "mean.Bovar_6KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Bovar_6KmDispDist_HSI_mean_importance$model = "Bovar_6KmDispDist"

Bovar_8KmDispDist_HSI_mean_importance = data.frame(mean(Bovar_8KmDispDist_HSI_tab$imp))
Bovar_8KmDispDist_HSI_mean_importance <- setnames(Bovar_8KmDispDist_HSI_mean_importance, "mean.Bovar_8KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Bovar_8KmDispDist_HSI_mean_importance$model = "Bovar_8KmDispDist"

Bovar_10KmDispDist_HSI_mean_importance = data.frame(mean(Bovar_10KmDispDist_HSI_tab$imp))
Bovar_10KmDispDist_HSI_mean_importance <- setnames(Bovar_10KmDispDist_HSI_mean_importance, "mean.Bovar_10KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Bovar_10KmDispDist_HSI_mean_importance$model = "Bovar_10KmDispDist"

##join them
Bovar_HSI_mean_importance_tab = rbind(Bovar_300mDispDist_HSI_mean_importance, Bovar_1KmDispDist_HSI_mean_importance, Bovar_defaultDispDist_HSI_mean_importance,
                                      Bovar_2KmDispDist_HSI_mean_importance, Bovar_6KmDispDist_HSI_mean_importance,
                                      Bovar_8KmDispDist_HSI_mean_importance, Bovar_10KmDispDist_HSI_mean_importance)


## Epcal

Epcal_defaultDispDist_HSI_mean_importance = data.frame(mean(Epcal_defaultDispDist_HSI_tab$imp))
Epcal_defaultDispDist_HSI_mean_importance <- setnames(Epcal_defaultDispDist_HSI_mean_importance, "mean.Epcal_defaultDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Epcal_defaultDispDist_HSI_mean_importance$model = "Epcal_defaultDispDist"

Epcal_300mDispDist_HSI_mean_importance = data.frame(mean(Epcal_300mDispDist_HSI_tab$imp))
Epcal_300mDispDist_HSI_mean_importance <- setnames(Epcal_300mDispDist_HSI_mean_importance, "mean.Epcal_300mDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Epcal_300mDispDist_HSI_mean_importance$model = "Epcal_300mDispDist"

Epcal_1KmDispDist_HSI_mean_importance = data.frame(mean(Epcal_1KmDispDist_HSI_tab$imp))
Epcal_1KmDispDist_HSI_mean_importance <- setnames(Epcal_1KmDispDist_HSI_mean_importance, "mean.Epcal_1KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Epcal_1KmDispDist_HSI_mean_importance$model = "Epcal_1KmDispDist"

Epcal_2KmDispDist_HSI_mean_importance = data.frame(mean(Epcal_2KmDispDist_HSI_tab$imp))
Epcal_2KmDispDist_HSI_mean_importance <- setnames(Epcal_2KmDispDist_HSI_mean_importance, "mean.Epcal_2KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Epcal_2KmDispDist_HSI_mean_importance$model = "Epcal_2KmDispDist"

Epcal_4KmDispDist_HSI_mean_importance = data.frame(mean(Epcal_4KmDispDist_HSI_tab$imp))
Epcal_4KmDispDist_HSI_mean_importance <- setnames(Epcal_4KmDispDist_HSI_mean_importance, "mean.Epcal_4KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Epcal_4KmDispDist_HSI_mean_importance$model = "Epcal_4KmDispDist"

Epcal_6KmDispDist_HSI_mean_importance = data.frame(mean(Epcal_6KmDispDist_HSI_tab$imp))
Epcal_6KmDispDist_HSI_mean_importance <- setnames(Epcal_6KmDispDist_HSI_mean_importance, "mean.Epcal_6KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Epcal_6KmDispDist_HSI_mean_importance$model = "Epcal_6KmDispDist"

Epcal_8KmDispDist_HSI_mean_importance = data.frame(mean(Epcal_8KmDispDist_HSI_tab$imp))
Epcal_8KmDispDist_HSI_mean_importance <- setnames(Epcal_8KmDispDist_HSI_mean_importance, "mean.Epcal_8KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Epcal_8KmDispDist_HSI_mean_importance$model = "Epcal_8KmDispDist"

Epcal_10KmDispDist_HSI_mean_importance = data.frame(mean(Epcal_10KmDispDist_HSI_tab$imp))
Epcal_10KmDispDist_HSI_mean_importance <- setnames(Epcal_10KmDispDist_HSI_mean_importance, "mean.Epcal_10KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Epcal_10KmDispDist_HSI_mean_importance$model = "Epcal_10KmDispDist"

##join them
Epcal_HSI_mean_importance_tab = rbind(Epcal_300mDispDist_HSI_mean_importance, Epcal_1KmDispDist_HSI_mean_importance, Epcal_defaultDispDist_HSI_mean_importance,
                                      Epcal_2KmDispDist_HSI_mean_importance, Epcal_4KmDispDist_HSI_mean_importance, Epcal_6KmDispDist_HSI_mean_importance,
                                      Epcal_8KmDispDist_HSI_mean_importance, Epcal_10KmDispDist_HSI_mean_importance)


## Peagg

Peagg_defaultDispDist_HSI_mean_importance = data.frame(mean(Peagg_defaultDispDist_HSI_tab$imp))
Peagg_defaultDispDist_HSI_mean_importance <- setnames(Peagg_defaultDispDist_HSI_mean_importance, "mean.Peagg_defaultDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Peagg_defaultDispDist_HSI_mean_importance$model = "Peagg_defaultDispDist"

Peagg_300mDispDist_HSI_mean_importance = data.frame(mean(Peagg_300mDispDist_HSI_tab$imp))
Peagg_300mDispDist_HSI_mean_importance <- setnames(Peagg_300mDispDist_HSI_mean_importance, "mean.Peagg_300mDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Peagg_300mDispDist_HSI_mean_importance$model = "Peagg_300mDispDist"

Peagg_1KmDispDist_HSI_mean_importance = data.frame(mean(Peagg_1KmDispDist_HSI_tab$imp))
Peagg_1KmDispDist_HSI_mean_importance <- setnames(Peagg_1KmDispDist_HSI_mean_importance, "mean.Peagg_1KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Peagg_1KmDispDist_HSI_mean_importance$model = "Peagg_1KmDispDist"

Peagg_2KmDispDist_HSI_mean_importance = data.frame(mean(Peagg_2KmDispDist_HSI_tab$imp))
Peagg_2KmDispDist_HSI_mean_importance <- setnames(Peagg_2KmDispDist_HSI_mean_importance, "mean.Peagg_2KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Peagg_2KmDispDist_HSI_mean_importance$model = "Peagg_2KmDispDist"

Peagg_4KmDispDist_HSI_mean_importance = data.frame(mean(Peagg_4KmDispDist_HSI_tab$imp))
Peagg_4KmDispDist_HSI_mean_importance <- setnames(Peagg_4KmDispDist_HSI_mean_importance, "mean.Peagg_4KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Peagg_4KmDispDist_HSI_mean_importance$model = "Peagg_4KmDispDist"

Peagg_6KmDispDist_HSI_mean_importance = data.frame(mean(Peagg_6KmDispDist_HSI_tab$imp))
Peagg_6KmDispDist_HSI_mean_importance <- setnames(Peagg_6KmDispDist_HSI_mean_importance, "mean.Peagg_6KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Peagg_6KmDispDist_HSI_mean_importance$model = "Peagg_6KmDispDist"

Peagg_8KmDispDist_HSI_mean_importance = data.frame(mean(Peagg_8KmDispDist_HSI_tab$imp))
Peagg_8KmDispDist_HSI_mean_importance <- setnames(Peagg_8KmDispDist_HSI_mean_importance, "mean.Peagg_8KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Peagg_8KmDispDist_HSI_mean_importance$model = "Peagg_8KmDispDist"

Peagg_10KmDispDist_HSI_mean_importance = data.frame(mean(Peagg_10KmDispDist_HSI_tab$imp))
Peagg_10KmDispDist_HSI_mean_importance <- setnames(Peagg_10KmDispDist_HSI_mean_importance, "mean.Peagg_10KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Peagg_10KmDispDist_HSI_mean_importance$model = "Peagg_10KmDispDist"

##join them
Peagg_HSI_mean_importance_tab = rbind(Peagg_300mDispDist_HSI_mean_importance, Peagg_1KmDispDist_HSI_mean_importance, Peagg_defaultDispDist_HSI_mean_importance,
                                      Peagg_2KmDispDist_HSI_mean_importance, Peagg_4KmDispDist_HSI_mean_importance, Peagg_6KmDispDist_HSI_mean_importance,
                                      Peagg_8KmDispDist_HSI_mean_importance, Peagg_10KmDispDist_HSI_mean_importance)


## Perid

Perid_defaultDispDist_HSI_mean_importance = data.frame(mean(Perid_defaultDispDist_HSI_tab$imp))
Perid_defaultDispDist_HSI_mean_importance <- setnames(Perid_defaultDispDist_HSI_mean_importance, "mean.Perid_defaultDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Perid_defaultDispDist_HSI_mean_importance$model = "Perid_defaultDispDist"

Perid_300mDispDist_HSI_mean_importance = data.frame(mean(Perid_300mDispDist_HSI_tab$imp))
Perid_300mDispDist_HSI_mean_importance <- setnames(Perid_300mDispDist_HSI_mean_importance, "mean.Perid_300mDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Perid_300mDispDist_HSI_mean_importance$model = "Perid_300mDispDist"

Perid_1KmDispDist_HSI_mean_importance = data.frame(mean(Perid_1KmDispDist_HSI_tab$imp))
Perid_1KmDispDist_HSI_mean_importance <- setnames(Perid_1KmDispDist_HSI_mean_importance, "mean.Perid_1KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Perid_1KmDispDist_HSI_mean_importance$model = "Perid_1KmDispDist"

Perid_2KmDispDist_HSI_mean_importance = data.frame(mean(Perid_2KmDispDist_HSI_tab$imp))
Perid_2KmDispDist_HSI_mean_importance <- setnames(Perid_2KmDispDist_HSI_mean_importance, "mean.Perid_2KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Perid_2KmDispDist_HSI_mean_importance$model = "Perid_2KmDispDist"

Perid_4KmDispDist_HSI_mean_importance = data.frame(mean(Perid_4KmDispDist_HSI_tab$imp))
Perid_4KmDispDist_HSI_mean_importance <- setnames(Perid_4KmDispDist_HSI_mean_importance, "mean.Perid_4KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Perid_4KmDispDist_HSI_mean_importance$model = "Perid_4KmDispDist"

Perid_6KmDispDist_HSI_mean_importance = data.frame(mean(Perid_6KmDispDist_HSI_tab$imp))
Perid_6KmDispDist_HSI_mean_importance <- setnames(Perid_6KmDispDist_HSI_mean_importance, "mean.Perid_6KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Perid_6KmDispDist_HSI_mean_importance$model = "Perid_6KmDispDist"

Perid_8KmDispDist_HSI_mean_importance = data.frame(mean(Perid_8KmDispDist_HSI_tab$imp))
Perid_8KmDispDist_HSI_mean_importance <- setnames(Perid_8KmDispDist_HSI_mean_importance, "mean.Perid_8KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Perid_8KmDispDist_HSI_mean_importance$model = "Perid_8KmDispDist"

Perid_10KmDispDist_HSI_mean_importance = data.frame(mean(Perid_10KmDispDist_HSI_tab$imp))
Perid_10KmDispDist_HSI_mean_importance <- setnames(Perid_10KmDispDist_HSI_mean_importance, "mean.Perid_10KmDispDist_HSI_tab.imp.","HSI_Mean_importance")
#Set identifier labels
Perid_10KmDispDist_HSI_mean_importance$model = "Perid_10KmDispDist"

##join them
Perid_HSI_mean_importance_tab = rbind(Perid_300mDispDist_HSI_mean_importance, Perid_1KmDispDist_HSI_mean_importance, Perid_defaultDispDist_HSI_mean_importance,
                                      Perid_2KmDispDist_HSI_mean_importance, Perid_4KmDispDist_HSI_mean_importance, Perid_6KmDispDist_HSI_mean_importance,
                                      Perid_8KmDispDist_HSI_mean_importance, Perid_10KmDispDist_HSI_mean_importance)


#### Join imp_tabs of all the species into one #################################

All_Spp_HSI_mean_importance_tab = rbind(Alobs_HSI_mean_importance_tab, Bovar_HSI_mean_importance_tab, Epcal_HSI_mean_importance_tab,
                                        Hyarb_HSI_mean_importance_tab, Peagg_HSI_mean_importance_tab, Perid_HSI_mean_importance_tab)

#### Join to main df with AUC cv & components ####

All_Spp_meanAUCcv_Components <- merge(All_Spp_meanAUCcv_Components, All_Spp_HSI_mean_importance_tab, by = "model")
head(All_Spp_meanAUCcv_Components)
All_Spp_meanAUCcv_Components$Max_disp_dist_m.y <- NULL
All_Spp_meanAUCcv_Components <- setnames(All_Spp_meanAUCcv_Components, "Max_disp_dist_m.x", "Max_disp_dist_m")

### Plot importance~distance
ggplot(All_Spp_meanAUCcv_Components, aes(x=Max_disp_dist_m, y = HSI_Mean_importance, color = Species, group = Species)) +
  geom_point(size = 3) + geom_line(size = 1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Mean importance of habitat suitability index")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

## Plot AUC-cv ~ importance
ggplot(All_Spp_meanAUCcv_Components, aes(x=HSI_Mean_importance, y = Mean_AUC_cv, color = Species, group = Species)) +
  geom_point() + geom_line()



###############################################################################################################################
#### Get 3rd Order Neighborhood (EgoSize) for all dists & species into main df and do multiplot of it vs Distance #############

# Hyarb
Hyarb_defaultDispDist_EgoSize_mean_importance = data.frame(mean(Hyarb_defaultDispDist_EgoSize_tab$imp))
Hyarb_defaultDispDist_EgoSize_mean_importance <- setnames(Hyarb_defaultDispDist_EgoSize_mean_importance, "mean.Hyarb_defaultDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Hyarb_defaultDispDist_EgoSize_mean_importance$model = "Hyarb_defaultDispDist"

Hyarb_300mDispDist_EgoSize_mean_importance = data.frame(mean(Hyarb_300mDispDist_EgoSize_tab$imp))
Hyarb_300mDispDist_EgoSize_mean_importance <- setnames(Hyarb_300mDispDist_EgoSize_mean_importance, "mean.Hyarb_300mDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Hyarb_300mDispDist_EgoSize_mean_importance$model = "Hyarb_300mDispDist"

Hyarb_1KmDispDist_EgoSize_mean_importance = data.frame(mean(Hyarb_1KmDispDist_EgoSize_tab$imp))
Hyarb_1KmDispDist_EgoSize_mean_importance <- setnames(Hyarb_1KmDispDist_EgoSize_mean_importance, "mean.Hyarb_1KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Hyarb_1KmDispDist_EgoSize_mean_importance$model = "Hyarb_1KmDispDist"

Hyarb_2KmDispDist_EgoSize_mean_importance = data.frame(mean(Hyarb_2KmDispDist_EgoSize_tab$imp))
Hyarb_2KmDispDist_EgoSize_mean_importance <- setnames(Hyarb_2KmDispDist_EgoSize_mean_importance, "mean.Hyarb_2KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Hyarb_2KmDispDist_EgoSize_mean_importance$model = "Hyarb_2KmDispDist"

Hyarb_4KmDispDist_EgoSize_mean_importance = data.frame(mean(Hyarb_4KmDispDist_EgoSize_tab$imp))
Hyarb_4KmDispDist_EgoSize_mean_importance <- setnames(Hyarb_4KmDispDist_EgoSize_mean_importance, "mean.Hyarb_4KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Hyarb_4KmDispDist_EgoSize_mean_importance$model = "Hyarb_4KmDispDist"

Hyarb_6KmDispDist_EgoSize_mean_importance = data.frame(mean(Hyarb_6KmDispDist_EgoSize_tab$imp))
Hyarb_6KmDispDist_EgoSize_mean_importance <- setnames(Hyarb_6KmDispDist_EgoSize_mean_importance, "mean.Hyarb_6KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Hyarb_6KmDispDist_EgoSize_mean_importance$model = "Hyarb_6KmDispDist"

Hyarb_8KmDispDist_EgoSize_mean_importance = data.frame(mean(Hyarb_8KmDispDist_EgoSize_tab$imp))
Hyarb_8KmDispDist_EgoSize_mean_importance <- setnames(Hyarb_8KmDispDist_EgoSize_mean_importance, "mean.Hyarb_8KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Hyarb_8KmDispDist_EgoSize_mean_importance$model = "Hyarb_8KmDispDist"

Hyarb_10KmDispDist_EgoSize_mean_importance = data.frame(mean(Hyarb_10KmDispDist_EgoSize_tab$imp))
Hyarb_10KmDispDist_EgoSize_mean_importance <- setnames(Hyarb_10KmDispDist_EgoSize_mean_importance, "mean.Hyarb_10KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Hyarb_10KmDispDist_EgoSize_mean_importance$model = "Hyarb_10KmDispDist"

##join them
Hyarb_EgoSize_mean_importance_tab = rbind(Hyarb_300mDispDist_EgoSize_mean_importance, Hyarb_1KmDispDist_EgoSize_mean_importance, Hyarb_defaultDispDist_EgoSize_mean_importance,
                                          Hyarb_2KmDispDist_EgoSize_mean_importance, Hyarb_4KmDispDist_EgoSize_mean_importance, Hyarb_6KmDispDist_EgoSize_mean_importance,
                                          Hyarb_8KmDispDist_EgoSize_mean_importance, Hyarb_10KmDispDist_EgoSize_mean_importance)


# Alobs
Alobs_defaultDispDist_EgoSize_mean_importance = data.frame(mean(Alobs_defaultDispDist_EgoSize_tab$imp))
Alobs_defaultDispDist_EgoSize_mean_importance <- setnames(Alobs_defaultDispDist_EgoSize_mean_importance, "mean.Alobs_defaultDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Alobs_defaultDispDist_EgoSize_mean_importance$model = "Alobs_defaultDispDist"

Alobs_300mDispDist_EgoSize_mean_importance = data.frame(mean(Alobs_300mDispDist_EgoSize_tab$imp))
Alobs_300mDispDist_EgoSize_mean_importance <- setnames(Alobs_300mDispDist_EgoSize_mean_importance, "mean.Alobs_300mDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Alobs_300mDispDist_EgoSize_mean_importance$model = "Alobs_300mDispDist"

Alobs_1KmDispDist_EgoSize_mean_importance = data.frame(mean(Alobs_1KmDispDist_EgoSize_tab$imp))
Alobs_1KmDispDist_EgoSize_mean_importance <- setnames(Alobs_1KmDispDist_EgoSize_mean_importance, "mean.Alobs_1KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Alobs_1KmDispDist_EgoSize_mean_importance$model = "Alobs_1KmDispDist"

Alobs_2KmDispDist_EgoSize_mean_importance = data.frame(mean(Alobs_2KmDispDist_EgoSize_tab$imp))
Alobs_2KmDispDist_EgoSize_mean_importance <- setnames(Alobs_2KmDispDist_EgoSize_mean_importance, "mean.Alobs_2KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Alobs_2KmDispDist_EgoSize_mean_importance$model = "Alobs_2KmDispDist"

Alobs_4KmDispDist_EgoSize_mean_importance = data.frame(mean(Alobs_4KmDispDist_EgoSize_tab$imp))
Alobs_4KmDispDist_EgoSize_mean_importance <- setnames(Alobs_4KmDispDist_EgoSize_mean_importance, "mean.Alobs_4KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Alobs_4KmDispDist_EgoSize_mean_importance$model = "Alobs_4KmDispDist"

Alobs_6KmDispDist_EgoSize_mean_importance = data.frame(mean(Alobs_6KmDispDist_EgoSize_tab$imp))
Alobs_6KmDispDist_EgoSize_mean_importance <- setnames(Alobs_6KmDispDist_EgoSize_mean_importance, "mean.Alobs_6KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Alobs_6KmDispDist_EgoSize_mean_importance$model = "Alobs_6KmDispDist"

Alobs_8KmDispDist_EgoSize_mean_importance = data.frame(mean(Alobs_8KmDispDist_EgoSize_tab$imp))
Alobs_8KmDispDist_EgoSize_mean_importance <- setnames(Alobs_8KmDispDist_EgoSize_mean_importance, "mean.Alobs_8KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Alobs_8KmDispDist_EgoSize_mean_importance$model = "Alobs_8KmDispDist"

Alobs_10KmDispDist_EgoSize_mean_importance = data.frame(mean(Alobs_10KmDispDist_EgoSize_tab$imp))
Alobs_10KmDispDist_EgoSize_mean_importance <- setnames(Alobs_10KmDispDist_EgoSize_mean_importance, "mean.Alobs_10KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Alobs_10KmDispDist_EgoSize_mean_importance$model = "Alobs_10KmDispDist"

##join them
Alobs_EgoSize_mean_importance_tab = rbind(Alobs_300mDispDist_EgoSize_mean_importance, Alobs_1KmDispDist_EgoSize_mean_importance, Alobs_defaultDispDist_EgoSize_mean_importance,
                                          Alobs_2KmDispDist_EgoSize_mean_importance, Alobs_4KmDispDist_EgoSize_mean_importance, Alobs_6KmDispDist_EgoSize_mean_importance,
                                          Alobs_8KmDispDist_EgoSize_mean_importance, Alobs_10KmDispDist_EgoSize_mean_importance)


# Bovar
Bovar_defaultDispDist_EgoSize_mean_importance = data.frame(mean(Bovar_defaultDispDist_EgoSize_tab$imp))
Bovar_defaultDispDist_EgoSize_mean_importance <- setnames(Bovar_defaultDispDist_EgoSize_mean_importance, "mean.Bovar_defaultDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Bovar_defaultDispDist_EgoSize_mean_importance$model = "Bovar_defaultDispDist"

Bovar_300mDispDist_EgoSize_mean_importance = data.frame(mean(Bovar_300mDispDist_EgoSize_tab$imp))
Bovar_300mDispDist_EgoSize_mean_importance <- setnames(Bovar_300mDispDist_EgoSize_mean_importance, "mean.Bovar_300mDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Bovar_300mDispDist_EgoSize_mean_importance$model = "Bovar_300mDispDist"

Bovar_1KmDispDist_EgoSize_mean_importance = data.frame(mean(Bovar_1KmDispDist_EgoSize_tab$imp))
Bovar_1KmDispDist_EgoSize_mean_importance <- setnames(Bovar_1KmDispDist_EgoSize_mean_importance, "mean.Bovar_1KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Bovar_1KmDispDist_EgoSize_mean_importance$model = "Bovar_1KmDispDist"

Bovar_2KmDispDist_EgoSize_mean_importance = data.frame(mean(Bovar_2KmDispDist_EgoSize_tab$imp))
Bovar_2KmDispDist_EgoSize_mean_importance <- setnames(Bovar_2KmDispDist_EgoSize_mean_importance, "mean.Bovar_2KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Bovar_2KmDispDist_EgoSize_mean_importance$model = "Bovar_2KmDispDist"

#Default for Bovar = 4Km
# Bovar_4KmDispDist_EgoSize_mean_importance = data.frame(mean(Bovar_4KmDispDist_EgoSize_tab$imp))
# Bovar_4KmDispDist_EgoSize_mean_importance <- setnames(Bovar_4KmDispDist_EgoSize_mean_importance, "mean.Bovar_4KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
# #Set identifier labels
# Bovar_4KmDispDist_EgoSize_mean_importance$model = "Bovar_4KmDispDist"

Bovar_6KmDispDist_EgoSize_mean_importance = data.frame(mean(Bovar_6KmDispDist_EgoSize_tab$imp))
Bovar_6KmDispDist_EgoSize_mean_importance <- setnames(Bovar_6KmDispDist_EgoSize_mean_importance, "mean.Bovar_6KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Bovar_6KmDispDist_EgoSize_mean_importance$model = "Bovar_6KmDispDist"

Bovar_8KmDispDist_EgoSize_mean_importance = data.frame(mean(Bovar_8KmDispDist_EgoSize_tab$imp))
Bovar_8KmDispDist_EgoSize_mean_importance <- setnames(Bovar_8KmDispDist_EgoSize_mean_importance, "mean.Bovar_8KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Bovar_8KmDispDist_EgoSize_mean_importance$model = "Bovar_8KmDispDist"

Bovar_10KmDispDist_EgoSize_mean_importance = data.frame(mean(Bovar_10KmDispDist_EgoSize_tab$imp))
Bovar_10KmDispDist_EgoSize_mean_importance <- setnames(Bovar_10KmDispDist_EgoSize_mean_importance, "mean.Bovar_10KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Bovar_10KmDispDist_EgoSize_mean_importance$model = "Bovar_10KmDispDist"

##join them
Bovar_EgoSize_mean_importance_tab = rbind(Bovar_300mDispDist_EgoSize_mean_importance, Bovar_1KmDispDist_EgoSize_mean_importance, Bovar_defaultDispDist_EgoSize_mean_importance,
                                          Bovar_2KmDispDist_EgoSize_mean_importance, Bovar_6KmDispDist_EgoSize_mean_importance,
                                          Bovar_8KmDispDist_EgoSize_mean_importance, Bovar_10KmDispDist_EgoSize_mean_importance)


# Epcal
Epcal_defaultDispDist_EgoSize_mean_importance = data.frame(mean(Epcal_defaultDispDist_EgoSize_tab$imp))
Epcal_defaultDispDist_EgoSize_mean_importance <- setnames(Epcal_defaultDispDist_EgoSize_mean_importance, "mean.Epcal_defaultDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Epcal_defaultDispDist_EgoSize_mean_importance$model = "Epcal_defaultDispDist"

Epcal_300mDispDist_EgoSize_mean_importance = data.frame(mean(Epcal_300mDispDist_EgoSize_tab$imp))
Epcal_300mDispDist_EgoSize_mean_importance <- setnames(Epcal_300mDispDist_EgoSize_mean_importance, "mean.Epcal_300mDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Epcal_300mDispDist_EgoSize_mean_importance$model = "Epcal_300mDispDist"

Epcal_1KmDispDist_EgoSize_mean_importance = data.frame(mean(Epcal_1KmDispDist_EgoSize_tab$imp))
Epcal_1KmDispDist_EgoSize_mean_importance <- setnames(Epcal_1KmDispDist_EgoSize_mean_importance, "mean.Epcal_1KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Epcal_1KmDispDist_EgoSize_mean_importance$model = "Epcal_1KmDispDist"

Epcal_2KmDispDist_EgoSize_mean_importance = data.frame(mean(Epcal_2KmDispDist_EgoSize_tab$imp))
Epcal_2KmDispDist_EgoSize_mean_importance <- setnames(Epcal_2KmDispDist_EgoSize_mean_importance, "mean.Epcal_2KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Epcal_2KmDispDist_EgoSize_mean_importance$model = "Epcal_2KmDispDist"

Epcal_4KmDispDist_EgoSize_mean_importance = data.frame(mean(Epcal_4KmDispDist_EgoSize_tab$imp))
Epcal_4KmDispDist_EgoSize_mean_importance <- setnames(Epcal_4KmDispDist_EgoSize_mean_importance, "mean.Epcal_4KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Epcal_4KmDispDist_EgoSize_mean_importance$model = "Epcal_4KmDispDist"

Epcal_6KmDispDist_EgoSize_mean_importance = data.frame(mean(Epcal_6KmDispDist_EgoSize_tab$imp))
Epcal_6KmDispDist_EgoSize_mean_importance <- setnames(Epcal_6KmDispDist_EgoSize_mean_importance, "mean.Epcal_6KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Epcal_6KmDispDist_EgoSize_mean_importance$model = "Epcal_6KmDispDist"

Epcal_8KmDispDist_EgoSize_mean_importance = data.frame(mean(Epcal_8KmDispDist_EgoSize_tab$imp))
Epcal_8KmDispDist_EgoSize_mean_importance <- setnames(Epcal_8KmDispDist_EgoSize_mean_importance, "mean.Epcal_8KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Epcal_8KmDispDist_EgoSize_mean_importance$model = "Epcal_8KmDispDist"

Epcal_10KmDispDist_EgoSize_mean_importance = data.frame(mean(Epcal_10KmDispDist_EgoSize_tab$imp))
Epcal_10KmDispDist_EgoSize_mean_importance <- setnames(Epcal_10KmDispDist_EgoSize_mean_importance, "mean.Epcal_10KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Epcal_10KmDispDist_EgoSize_mean_importance$model = "Epcal_10KmDispDist"

##join them
Epcal_EgoSize_mean_importance_tab = rbind(Epcal_300mDispDist_EgoSize_mean_importance, Epcal_1KmDispDist_EgoSize_mean_importance, Epcal_defaultDispDist_EgoSize_mean_importance,
                                          Epcal_2KmDispDist_EgoSize_mean_importance, Epcal_4KmDispDist_EgoSize_mean_importance, Epcal_6KmDispDist_EgoSize_mean_importance,
                                          Epcal_8KmDispDist_EgoSize_mean_importance, Epcal_10KmDispDist_EgoSize_mean_importance)


# Peagg
Peagg_defaultDispDist_EgoSize_mean_importance = data.frame(mean(Peagg_defaultDispDist_EgoSize_tab$imp))
Peagg_defaultDispDist_EgoSize_mean_importance <- setnames(Peagg_defaultDispDist_EgoSize_mean_importance, "mean.Peagg_defaultDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Peagg_defaultDispDist_EgoSize_mean_importance$model = "Peagg_defaultDispDist"

Peagg_300mDispDist_EgoSize_mean_importance = data.frame(mean(Peagg_300mDispDist_EgoSize_tab$imp))
Peagg_300mDispDist_EgoSize_mean_importance <- setnames(Peagg_300mDispDist_EgoSize_mean_importance, "mean.Peagg_300mDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Peagg_300mDispDist_EgoSize_mean_importance$model = "Peagg_300mDispDist"

Peagg_1KmDispDist_EgoSize_mean_importance = data.frame(mean(Peagg_1KmDispDist_EgoSize_tab$imp))
Peagg_1KmDispDist_EgoSize_mean_importance <- setnames(Peagg_1KmDispDist_EgoSize_mean_importance, "mean.Peagg_1KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Peagg_1KmDispDist_EgoSize_mean_importance$model = "Peagg_1KmDispDist"

Peagg_2KmDispDist_EgoSize_mean_importance = data.frame(mean(Peagg_2KmDispDist_EgoSize_tab$imp))
Peagg_2KmDispDist_EgoSize_mean_importance <- setnames(Peagg_2KmDispDist_EgoSize_mean_importance, "mean.Peagg_2KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Peagg_2KmDispDist_EgoSize_mean_importance$model = "Peagg_2KmDispDist"

Peagg_4KmDispDist_EgoSize_mean_importance = data.frame(mean(Peagg_4KmDispDist_EgoSize_tab$imp))
Peagg_4KmDispDist_EgoSize_mean_importance <- setnames(Peagg_4KmDispDist_EgoSize_mean_importance, "mean.Peagg_4KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Peagg_4KmDispDist_EgoSize_mean_importance$model = "Peagg_4KmDispDist"

Peagg_6KmDispDist_EgoSize_mean_importance = data.frame(mean(Peagg_6KmDispDist_EgoSize_tab$imp))
Peagg_6KmDispDist_EgoSize_mean_importance <- setnames(Peagg_6KmDispDist_EgoSize_mean_importance, "mean.Peagg_6KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Peagg_6KmDispDist_EgoSize_mean_importance$model = "Peagg_6KmDispDist"

Peagg_8KmDispDist_EgoSize_mean_importance = data.frame(mean(Peagg_8KmDispDist_EgoSize_tab$imp))
Peagg_8KmDispDist_EgoSize_mean_importance <- setnames(Peagg_8KmDispDist_EgoSize_mean_importance, "mean.Peagg_8KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Peagg_8KmDispDist_EgoSize_mean_importance$model = "Peagg_8KmDispDist"

Peagg_10KmDispDist_EgoSize_mean_importance = data.frame(mean(Peagg_10KmDispDist_EgoSize_tab$imp))
Peagg_10KmDispDist_EgoSize_mean_importance <- setnames(Peagg_10KmDispDist_EgoSize_mean_importance, "mean.Peagg_10KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Peagg_10KmDispDist_EgoSize_mean_importance$model = "Peagg_10KmDispDist"

##join them
Peagg_EgoSize_mean_importance_tab = rbind(Peagg_300mDispDist_EgoSize_mean_importance, Peagg_1KmDispDist_EgoSize_mean_importance, Peagg_defaultDispDist_EgoSize_mean_importance,
                                          Peagg_2KmDispDist_EgoSize_mean_importance, Peagg_4KmDispDist_EgoSize_mean_importance, Peagg_6KmDispDist_EgoSize_mean_importance,
                                          Peagg_8KmDispDist_EgoSize_mean_importance, Peagg_10KmDispDist_EgoSize_mean_importance)


# Perid
Perid_defaultDispDist_EgoSize_mean_importance = data.frame(mean(Perid_defaultDispDist_EgoSize_tab$imp))
Perid_defaultDispDist_EgoSize_mean_importance <- setnames(Perid_defaultDispDist_EgoSize_mean_importance, "mean.Perid_defaultDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Perid_defaultDispDist_EgoSize_mean_importance$model = "Perid_defaultDispDist"

Perid_300mDispDist_EgoSize_mean_importance = data.frame(mean(Perid_300mDispDist_EgoSize_tab$imp))
Perid_300mDispDist_EgoSize_mean_importance <- setnames(Perid_300mDispDist_EgoSize_mean_importance, "mean.Perid_300mDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Perid_300mDispDist_EgoSize_mean_importance$model = "Perid_300mDispDist"

Perid_1KmDispDist_EgoSize_mean_importance = data.frame(mean(Perid_1KmDispDist_EgoSize_tab$imp))
Perid_1KmDispDist_EgoSize_mean_importance <- setnames(Perid_1KmDispDist_EgoSize_mean_importance, "mean.Perid_1KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Perid_1KmDispDist_EgoSize_mean_importance$model = "Perid_1KmDispDist"

Perid_2KmDispDist_EgoSize_mean_importance = data.frame(mean(Perid_2KmDispDist_EgoSize_tab$imp))
Perid_2KmDispDist_EgoSize_mean_importance <- setnames(Perid_2KmDispDist_EgoSize_mean_importance, "mean.Perid_2KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Perid_2KmDispDist_EgoSize_mean_importance$model = "Perid_2KmDispDist"

Perid_4KmDispDist_EgoSize_mean_importance = data.frame(mean(Perid_4KmDispDist_EgoSize_tab$imp))
Perid_4KmDispDist_EgoSize_mean_importance <- setnames(Perid_4KmDispDist_EgoSize_mean_importance, "mean.Perid_4KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Perid_4KmDispDist_EgoSize_mean_importance$model = "Perid_4KmDispDist"

Perid_6KmDispDist_EgoSize_mean_importance = data.frame(mean(Perid_6KmDispDist_EgoSize_tab$imp))
Perid_6KmDispDist_EgoSize_mean_importance <- setnames(Perid_6KmDispDist_EgoSize_mean_importance, "mean.Perid_6KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Perid_6KmDispDist_EgoSize_mean_importance$model = "Perid_6KmDispDist"

Perid_8KmDispDist_EgoSize_mean_importance = data.frame(mean(Perid_8KmDispDist_EgoSize_tab$imp))
Perid_8KmDispDist_EgoSize_mean_importance <- setnames(Perid_8KmDispDist_EgoSize_mean_importance, "mean.Perid_8KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Perid_8KmDispDist_EgoSize_mean_importance$model = "Perid_8KmDispDist"

Perid_10KmDispDist_EgoSize_mean_importance = data.frame(mean(Perid_10KmDispDist_EgoSize_tab$imp))
Perid_10KmDispDist_EgoSize_mean_importance <- setnames(Perid_10KmDispDist_EgoSize_mean_importance, "mean.Perid_10KmDispDist_EgoSize_tab.imp.","EgoSize_Mean_importance")
#Set identifier labels
Perid_10KmDispDist_EgoSize_mean_importance$model = "Perid_10KmDispDist"

##join them
Perid_EgoSize_mean_importance_tab = rbind(Perid_300mDispDist_EgoSize_mean_importance, Perid_1KmDispDist_EgoSize_mean_importance, Perid_defaultDispDist_EgoSize_mean_importance,
                                          Perid_2KmDispDist_EgoSize_mean_importance, Perid_4KmDispDist_EgoSize_mean_importance, Perid_6KmDispDist_EgoSize_mean_importance,
                                          Perid_8KmDispDist_EgoSize_mean_importance, Perid_10KmDispDist_EgoSize_mean_importance)


#### Join imp_tabs of all the species into one ####
All_Spp_EgoSize_mean_importance_tab = rbind(Alobs_EgoSize_mean_importance_tab, Bovar_EgoSize_mean_importance_tab, Epcal_EgoSize_mean_importance_tab,
                                            Hyarb_EgoSize_mean_importance_tab, Peagg_EgoSize_mean_importance_tab, Perid_EgoSize_mean_importance_tab)
All_Spp_EgoSize_mean_importance_tab <- setnames(All_Spp_EgoSize_mean_importance_tab, "EgoSize_Mean_importance", "3rd_Order_Neigh_Mean_Importance")


#### Join to main df with AUC cv & components ####
All_Spp_meanAUCcv_Components <- merge(All_Spp_meanAUCcv_Components, All_Spp_EgoSize_mean_importance_tab, by = "model")
head(All_Spp_meanAUCcv_Components)

# All_Spp_meanAUCcv_Components$EgoSize_Mean_importance <- NULL
# All_Spp_meanAUCcv_Components <- setnames(All_Spp_meanAUCcv_Components, "3rd_Order_Neigh_Mean_Importance", "Third_Order_Neigh_Mean_Importance")

### Plot importance~distance
ggplot(All_Spp_meanAUCcv_Components, aes(x=Max_disp_dist_m, y = Third_Order_Neigh_Mean_Importance, color = Species, group = Species)) +
  geom_point(size = 3) + geom_line(size = 1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Mean importance of 3rd order neighborhood")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

### Plot AUC-cv~importance
ggplot(All_Spp_meanAUCcv_Components, aes(x=Third_Order_Neigh_Mean_Importance, y = Mean_AUC_cv, color = Species, group = Species)) +
  geom_point() + geom_line()



###############################################################################################################################
#### Get habAv for all dists & species into main df and do multiplot of it vs Distance #############

# Hyarb
Hyarb_defaultDispDist_habAv_mean_importance = data.frame(mean(Hyarb_defaultDispDist_habAv_tab$imp))
Hyarb_defaultDispDist_habAv_mean_importance <- setnames(Hyarb_defaultDispDist_habAv_mean_importance, "mean.Hyarb_defaultDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Hyarb_defaultDispDist_habAv_mean_importance$model = "Hyarb_defaultDispDist"

Hyarb_300mDispDist_habAv_mean_importance = data.frame(mean(Hyarb_300mDispDist_habAv_tab$imp))
Hyarb_300mDispDist_habAv_mean_importance <- setnames(Hyarb_300mDispDist_habAv_mean_importance, "mean.Hyarb_300mDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Hyarb_300mDispDist_habAv_mean_importance$model = "Hyarb_300mDispDist"

Hyarb_1KmDispDist_habAv_mean_importance = data.frame(mean(Hyarb_1KmDispDist_habAv_tab$imp))
Hyarb_1KmDispDist_habAv_mean_importance <- setnames(Hyarb_1KmDispDist_habAv_mean_importance, "mean.Hyarb_1KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Hyarb_1KmDispDist_habAv_mean_importance$model = "Hyarb_1KmDispDist"

Hyarb_2KmDispDist_habAv_mean_importance = data.frame(mean(Hyarb_2KmDispDist_habAv_tab$imp))
Hyarb_2KmDispDist_habAv_mean_importance <- setnames(Hyarb_2KmDispDist_habAv_mean_importance, "mean.Hyarb_2KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Hyarb_2KmDispDist_habAv_mean_importance$model = "Hyarb_2KmDispDist"

Hyarb_4KmDispDist_habAv_mean_importance = data.frame(mean(Hyarb_4KmDispDist_habAv_tab$imp))
Hyarb_4KmDispDist_habAv_mean_importance <- setnames(Hyarb_4KmDispDist_habAv_mean_importance, "mean.Hyarb_4KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Hyarb_4KmDispDist_habAv_mean_importance$model = "Hyarb_4KmDispDist"

Hyarb_6KmDispDist_habAv_mean_importance = data.frame(mean(Hyarb_6KmDispDist_habAv_tab$imp))
Hyarb_6KmDispDist_habAv_mean_importance <- setnames(Hyarb_6KmDispDist_habAv_mean_importance, "mean.Hyarb_6KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Hyarb_6KmDispDist_habAv_mean_importance$model = "Hyarb_6KmDispDist"

Hyarb_8KmDispDist_habAv_mean_importance = data.frame(mean(Hyarb_8KmDispDist_habAv_tab$imp))
Hyarb_8KmDispDist_habAv_mean_importance <- setnames(Hyarb_8KmDispDist_habAv_mean_importance, "mean.Hyarb_8KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Hyarb_8KmDispDist_habAv_mean_importance$model = "Hyarb_8KmDispDist"

Hyarb_10KmDispDist_habAv_mean_importance = data.frame(mean(Hyarb_10KmDispDist_habAv_tab$imp))
Hyarb_10KmDispDist_habAv_mean_importance <- setnames(Hyarb_10KmDispDist_habAv_mean_importance, "mean.Hyarb_10KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Hyarb_10KmDispDist_habAv_mean_importance$model = "Hyarb_10KmDispDist"

##join them
Hyarb_habAv_mean_importance_tab = rbind(Hyarb_300mDispDist_habAv_mean_importance, Hyarb_1KmDispDist_habAv_mean_importance, Hyarb_defaultDispDist_habAv_mean_importance,
                                        Hyarb_2KmDispDist_habAv_mean_importance, Hyarb_4KmDispDist_habAv_mean_importance, Hyarb_6KmDispDist_habAv_mean_importance,
                                        Hyarb_8KmDispDist_habAv_mean_importance, Hyarb_10KmDispDist_habAv_mean_importance)


# Alobs
Alobs_defaultDispDist_habAv_mean_importance = data.frame(mean(Alobs_defaultDispDist_habAv_tab$imp))
Alobs_defaultDispDist_habAv_mean_importance <- setnames(Alobs_defaultDispDist_habAv_mean_importance, "mean.Alobs_defaultDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Alobs_defaultDispDist_habAv_mean_importance$model = "Alobs_defaultDispDist"

Alobs_300mDispDist_habAv_mean_importance = data.frame(mean(Alobs_300mDispDist_habAv_tab$imp))
Alobs_300mDispDist_habAv_mean_importance <- setnames(Alobs_300mDispDist_habAv_mean_importance, "mean.Alobs_300mDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Alobs_300mDispDist_habAv_mean_importance$model = "Alobs_300mDispDist"

Alobs_1KmDispDist_habAv_mean_importance = data.frame(mean(Alobs_1KmDispDist_habAv_tab$imp))
Alobs_1KmDispDist_habAv_mean_importance <- setnames(Alobs_1KmDispDist_habAv_mean_importance, "mean.Alobs_1KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Alobs_1KmDispDist_habAv_mean_importance$model = "Alobs_1KmDispDist"

Alobs_2KmDispDist_habAv_mean_importance = data.frame(mean(Alobs_2KmDispDist_habAv_tab$imp))
Alobs_2KmDispDist_habAv_mean_importance <- setnames(Alobs_2KmDispDist_habAv_mean_importance, "mean.Alobs_2KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Alobs_2KmDispDist_habAv_mean_importance$model = "Alobs_2KmDispDist"

Alobs_4KmDispDist_habAv_mean_importance = data.frame(mean(Alobs_4KmDispDist_habAv_tab$imp))
Alobs_4KmDispDist_habAv_mean_importance <- setnames(Alobs_4KmDispDist_habAv_mean_importance, "mean.Alobs_4KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Alobs_4KmDispDist_habAv_mean_importance$model = "Alobs_4KmDispDist"

Alobs_6KmDispDist_habAv_mean_importance = data.frame(mean(Alobs_6KmDispDist_habAv_tab$imp))
Alobs_6KmDispDist_habAv_mean_importance <- setnames(Alobs_6KmDispDist_habAv_mean_importance, "mean.Alobs_6KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Alobs_6KmDispDist_habAv_mean_importance$model = "Alobs_6KmDispDist"

Alobs_8KmDispDist_habAv_mean_importance = data.frame(mean(Alobs_8KmDispDist_habAv_tab$imp))
Alobs_8KmDispDist_habAv_mean_importance <- setnames(Alobs_8KmDispDist_habAv_mean_importance, "mean.Alobs_8KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Alobs_8KmDispDist_habAv_mean_importance$model = "Alobs_8KmDispDist"

Alobs_10KmDispDist_habAv_mean_importance = data.frame(mean(Alobs_10KmDispDist_habAv_tab$imp))
Alobs_10KmDispDist_habAv_mean_importance <- setnames(Alobs_10KmDispDist_habAv_mean_importance, "mean.Alobs_10KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Alobs_10KmDispDist_habAv_mean_importance$model = "Alobs_10KmDispDist"

##join them
Alobs_habAv_mean_importance_tab = rbind(Alobs_300mDispDist_habAv_mean_importance, Alobs_1KmDispDist_habAv_mean_importance, Alobs_defaultDispDist_habAv_mean_importance,
                                        Alobs_2KmDispDist_habAv_mean_importance, Alobs_4KmDispDist_habAv_mean_importance, Alobs_6KmDispDist_habAv_mean_importance,
                                        Alobs_8KmDispDist_habAv_mean_importance, Alobs_10KmDispDist_habAv_mean_importance)


# Bovar
Bovar_defaultDispDist_habAv_mean_importance = data.frame(mean(Bovar_defaultDispDist_habAv_tab$imp))
Bovar_defaultDispDist_habAv_mean_importance <- setnames(Bovar_defaultDispDist_habAv_mean_importance, "mean.Bovar_defaultDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Bovar_defaultDispDist_habAv_mean_importance$model = "Bovar_defaultDispDist"

Bovar_300mDispDist_habAv_mean_importance = data.frame(mean(Bovar_300mDispDist_habAv_tab$imp))
Bovar_300mDispDist_habAv_mean_importance <- setnames(Bovar_300mDispDist_habAv_mean_importance, "mean.Bovar_300mDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Bovar_300mDispDist_habAv_mean_importance$model = "Bovar_300mDispDist"

Bovar_1KmDispDist_habAv_mean_importance = data.frame(mean(Bovar_1KmDispDist_habAv_tab$imp))
Bovar_1KmDispDist_habAv_mean_importance <- setnames(Bovar_1KmDispDist_habAv_mean_importance, "mean.Bovar_1KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Bovar_1KmDispDist_habAv_mean_importance$model = "Bovar_1KmDispDist"

Bovar_2KmDispDist_habAv_mean_importance = data.frame(mean(Bovar_2KmDispDist_habAv_tab$imp))
Bovar_2KmDispDist_habAv_mean_importance <- setnames(Bovar_2KmDispDist_habAv_mean_importance, "mean.Bovar_2KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Bovar_2KmDispDist_habAv_mean_importance$model = "Bovar_2KmDispDist"

#Default for Bovar = 4Km
# Bovar_4KmDispDist_habAv_mean_importance = data.frame(mean(Bovar_4KmDispDist_habAv_tab$imp))
# Bovar_4KmDispDist_habAv_mean_importance <- setnames(Bovar_4KmDispDist_habAv_mean_importance, "mean.Bovar_4KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
# #Set identifier labels
# Bovar_4KmDispDist_habAv_mean_importance$model = "Bovar_4KmDispDist"

Bovar_6KmDispDist_habAv_mean_importance = data.frame(mean(Bovar_6KmDispDist_habAv_tab$imp))
Bovar_6KmDispDist_habAv_mean_importance <- setnames(Bovar_6KmDispDist_habAv_mean_importance, "mean.Bovar_6KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Bovar_6KmDispDist_habAv_mean_importance$model = "Bovar_6KmDispDist"

Bovar_8KmDispDist_habAv_mean_importance = data.frame(mean(Bovar_8KmDispDist_habAv_tab$imp))
Bovar_8KmDispDist_habAv_mean_importance <- setnames(Bovar_8KmDispDist_habAv_mean_importance, "mean.Bovar_8KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Bovar_8KmDispDist_habAv_mean_importance$model = "Bovar_8KmDispDist"

Bovar_10KmDispDist_habAv_mean_importance = data.frame(mean(Bovar_10KmDispDist_habAv_tab$imp))
Bovar_10KmDispDist_habAv_mean_importance <- setnames(Bovar_10KmDispDist_habAv_mean_importance, "mean.Bovar_10KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Bovar_10KmDispDist_habAv_mean_importance$model = "Bovar_10KmDispDist"

##join them
Bovar_habAv_mean_importance_tab = rbind(Bovar_300mDispDist_habAv_mean_importance, Bovar_1KmDispDist_habAv_mean_importance, Bovar_defaultDispDist_habAv_mean_importance,
                                        Bovar_2KmDispDist_habAv_mean_importance, Bovar_6KmDispDist_habAv_mean_importance,
                                        Bovar_8KmDispDist_habAv_mean_importance, Bovar_10KmDispDist_habAv_mean_importance)


# Epcal
Epcal_defaultDispDist_habAv_mean_importance = data.frame(mean(Epcal_defaultDispDist_habAv_tab$imp))
Epcal_defaultDispDist_habAv_mean_importance <- setnames(Epcal_defaultDispDist_habAv_mean_importance, "mean.Epcal_defaultDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Epcal_defaultDispDist_habAv_mean_importance$model = "Epcal_defaultDispDist"

Epcal_300mDispDist_habAv_mean_importance = data.frame(mean(Epcal_300mDispDist_habAv_tab$imp))
Epcal_300mDispDist_habAv_mean_importance <- setnames(Epcal_300mDispDist_habAv_mean_importance, "mean.Epcal_300mDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Epcal_300mDispDist_habAv_mean_importance$model = "Epcal_300mDispDist"

Epcal_1KmDispDist_habAv_mean_importance = data.frame(mean(Epcal_1KmDispDist_habAv_tab$imp))
Epcal_1KmDispDist_habAv_mean_importance <- setnames(Epcal_1KmDispDist_habAv_mean_importance, "mean.Epcal_1KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Epcal_1KmDispDist_habAv_mean_importance$model = "Epcal_1KmDispDist"

Epcal_2KmDispDist_habAv_mean_importance = data.frame(mean(Epcal_2KmDispDist_habAv_tab$imp))
Epcal_2KmDispDist_habAv_mean_importance <- setnames(Epcal_2KmDispDist_habAv_mean_importance, "mean.Epcal_2KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Epcal_2KmDispDist_habAv_mean_importance$model = "Epcal_2KmDispDist"

Epcal_4KmDispDist_habAv_mean_importance = data.frame(mean(Epcal_4KmDispDist_habAv_tab$imp))
Epcal_4KmDispDist_habAv_mean_importance <- setnames(Epcal_4KmDispDist_habAv_mean_importance, "mean.Epcal_4KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Epcal_4KmDispDist_habAv_mean_importance$model = "Epcal_4KmDispDist"

Epcal_6KmDispDist_habAv_mean_importance = data.frame(mean(Epcal_6KmDispDist_habAv_tab$imp))
Epcal_6KmDispDist_habAv_mean_importance <- setnames(Epcal_6KmDispDist_habAv_mean_importance, "mean.Epcal_6KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Epcal_6KmDispDist_habAv_mean_importance$model = "Epcal_6KmDispDist"

Epcal_8KmDispDist_habAv_mean_importance = data.frame(mean(Epcal_8KmDispDist_habAv_tab$imp))
Epcal_8KmDispDist_habAv_mean_importance <- setnames(Epcal_8KmDispDist_habAv_mean_importance, "mean.Epcal_8KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Epcal_8KmDispDist_habAv_mean_importance$model = "Epcal_8KmDispDist"

Epcal_10KmDispDist_habAv_mean_importance = data.frame(mean(Epcal_10KmDispDist_habAv_tab$imp))
Epcal_10KmDispDist_habAv_mean_importance <- setnames(Epcal_10KmDispDist_habAv_mean_importance, "mean.Epcal_10KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Epcal_10KmDispDist_habAv_mean_importance$model = "Epcal_10KmDispDist"

##join them
Epcal_habAv_mean_importance_tab = rbind(Epcal_300mDispDist_habAv_mean_importance, Epcal_1KmDispDist_habAv_mean_importance, Epcal_defaultDispDist_habAv_mean_importance,
                                        Epcal_2KmDispDist_habAv_mean_importance, Epcal_4KmDispDist_habAv_mean_importance, Epcal_6KmDispDist_habAv_mean_importance,
                                        Epcal_8KmDispDist_habAv_mean_importance, Epcal_10KmDispDist_habAv_mean_importance)


# Peagg
Peagg_defaultDispDist_habAv_mean_importance = data.frame(mean(Peagg_defaultDispDist_habAv_tab$imp))
Peagg_defaultDispDist_habAv_mean_importance <- setnames(Peagg_defaultDispDist_habAv_mean_importance, "mean.Peagg_defaultDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Peagg_defaultDispDist_habAv_mean_importance$model = "Peagg_defaultDispDist"

Peagg_300mDispDist_habAv_mean_importance = data.frame(mean(Peagg_300mDispDist_habAv_tab$imp))
Peagg_300mDispDist_habAv_mean_importance <- setnames(Peagg_300mDispDist_habAv_mean_importance, "mean.Peagg_300mDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Peagg_300mDispDist_habAv_mean_importance$model = "Peagg_300mDispDist"

Peagg_1KmDispDist_habAv_mean_importance = data.frame(mean(Peagg_1KmDispDist_habAv_tab$imp))
Peagg_1KmDispDist_habAv_mean_importance <- setnames(Peagg_1KmDispDist_habAv_mean_importance, "mean.Peagg_1KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Peagg_1KmDispDist_habAv_mean_importance$model = "Peagg_1KmDispDist"

Peagg_2KmDispDist_habAv_mean_importance = data.frame(mean(Peagg_2KmDispDist_habAv_tab$imp))
Peagg_2KmDispDist_habAv_mean_importance <- setnames(Peagg_2KmDispDist_habAv_mean_importance, "mean.Peagg_2KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Peagg_2KmDispDist_habAv_mean_importance$model = "Peagg_2KmDispDist"

Peagg_4KmDispDist_habAv_mean_importance = data.frame(mean(Peagg_4KmDispDist_habAv_tab$imp))
Peagg_4KmDispDist_habAv_mean_importance <- setnames(Peagg_4KmDispDist_habAv_mean_importance, "mean.Peagg_4KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Peagg_4KmDispDist_habAv_mean_importance$model = "Peagg_4KmDispDist"

Peagg_6KmDispDist_habAv_mean_importance = data.frame(mean(Peagg_6KmDispDist_habAv_tab$imp))
Peagg_6KmDispDist_habAv_mean_importance <- setnames(Peagg_6KmDispDist_habAv_mean_importance, "mean.Peagg_6KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Peagg_6KmDispDist_habAv_mean_importance$model = "Peagg_6KmDispDist"

Peagg_8KmDispDist_habAv_mean_importance = data.frame(mean(Peagg_8KmDispDist_habAv_tab$imp))
Peagg_8KmDispDist_habAv_mean_importance <- setnames(Peagg_8KmDispDist_habAv_mean_importance, "mean.Peagg_8KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Peagg_8KmDispDist_habAv_mean_importance$model = "Peagg_8KmDispDist"

Peagg_10KmDispDist_habAv_mean_importance = data.frame(mean(Peagg_10KmDispDist_habAv_tab$imp))
Peagg_10KmDispDist_habAv_mean_importance <- setnames(Peagg_10KmDispDist_habAv_mean_importance, "mean.Peagg_10KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Peagg_10KmDispDist_habAv_mean_importance$model = "Peagg_10KmDispDist"

##join them
Peagg_habAv_mean_importance_tab = rbind(Peagg_300mDispDist_habAv_mean_importance, Peagg_1KmDispDist_habAv_mean_importance, Peagg_defaultDispDist_habAv_mean_importance,
                                        Peagg_2KmDispDist_habAv_mean_importance, Peagg_4KmDispDist_habAv_mean_importance, Peagg_6KmDispDist_habAv_mean_importance,
                                        Peagg_8KmDispDist_habAv_mean_importance, Peagg_10KmDispDist_habAv_mean_importance)


# Perid
Perid_defaultDispDist_habAv_mean_importance = data.frame(mean(Perid_defaultDispDist_habAv_tab$imp))
Perid_defaultDispDist_habAv_mean_importance <- setnames(Perid_defaultDispDist_habAv_mean_importance, "mean.Perid_defaultDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Perid_defaultDispDist_habAv_mean_importance$model = "Perid_defaultDispDist"

Perid_300mDispDist_habAv_mean_importance = data.frame(mean(Perid_300mDispDist_habAv_tab$imp))
Perid_300mDispDist_habAv_mean_importance <- setnames(Perid_300mDispDist_habAv_mean_importance, "mean.Perid_300mDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Perid_300mDispDist_habAv_mean_importance$model = "Perid_300mDispDist"

Perid_1KmDispDist_habAv_mean_importance = data.frame(mean(Perid_1KmDispDist_habAv_tab$imp))
Perid_1KmDispDist_habAv_mean_importance <- setnames(Perid_1KmDispDist_habAv_mean_importance, "mean.Perid_1KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Perid_1KmDispDist_habAv_mean_importance$model = "Perid_1KmDispDist"

Perid_2KmDispDist_habAv_mean_importance = data.frame(mean(Perid_2KmDispDist_habAv_tab$imp))
Perid_2KmDispDist_habAv_mean_importance <- setnames(Perid_2KmDispDist_habAv_mean_importance, "mean.Perid_2KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Perid_2KmDispDist_habAv_mean_importance$model = "Perid_2KmDispDist"

Perid_4KmDispDist_habAv_mean_importance = data.frame(mean(Perid_4KmDispDist_habAv_tab$imp))
Perid_4KmDispDist_habAv_mean_importance <- setnames(Perid_4KmDispDist_habAv_mean_importance, "mean.Perid_4KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Perid_4KmDispDist_habAv_mean_importance$model = "Perid_4KmDispDist"

Perid_6KmDispDist_habAv_mean_importance = data.frame(mean(Perid_6KmDispDist_habAv_tab$imp))
Perid_6KmDispDist_habAv_mean_importance <- setnames(Perid_6KmDispDist_habAv_mean_importance, "mean.Perid_6KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Perid_6KmDispDist_habAv_mean_importance$model = "Perid_6KmDispDist"

Perid_8KmDispDist_habAv_mean_importance = data.frame(mean(Perid_8KmDispDist_habAv_tab$imp))
Perid_8KmDispDist_habAv_mean_importance <- setnames(Perid_8KmDispDist_habAv_mean_importance, "mean.Perid_8KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Perid_8KmDispDist_habAv_mean_importance$model = "Perid_8KmDispDist"

Perid_10KmDispDist_habAv_mean_importance = data.frame(mean(Perid_10KmDispDist_habAv_tab$imp))
Perid_10KmDispDist_habAv_mean_importance <- setnames(Perid_10KmDispDist_habAv_mean_importance, "mean.Perid_10KmDispDist_habAv_tab.imp.","habAv_Mean_importance")
#Set identifier labels
Perid_10KmDispDist_habAv_mean_importance$model = "Perid_10KmDispDist"

##join them
Perid_habAv_mean_importance_tab = rbind(Perid_300mDispDist_habAv_mean_importance, Perid_1KmDispDist_habAv_mean_importance, Perid_defaultDispDist_habAv_mean_importance,
                                        Perid_2KmDispDist_habAv_mean_importance, Perid_4KmDispDist_habAv_mean_importance, Perid_6KmDispDist_habAv_mean_importance,
                                        Perid_8KmDispDist_habAv_mean_importance, Perid_10KmDispDist_habAv_mean_importance)


#### Join imp_tabs of all the species into one ####
All_Spp_habAv_mean_importance_tab = rbind(Alobs_habAv_mean_importance_tab, Bovar_habAv_mean_importance_tab, Epcal_habAv_mean_importance_tab,
                                          Hyarb_habAv_mean_importance_tab, Peagg_habAv_mean_importance_tab, Perid_habAv_mean_importance_tab)
# All_Spp_habAv_mean_importance_tab <- setnames(All_Spp_habAv_mean_importance_tab, "habAv_Mean_importance", "3rd_Order_Neigh_Mean_Importance")


#### Join to main df with AUC cv & components ####
All_Spp_meanAUCcv_Components <- merge(All_Spp_meanAUCcv_Components, All_Spp_habAv_mean_importance_tab, by = "model")
head(All_Spp_meanAUCcv_Components)

# All_Spp_meanAUCcv_Components$habAv_Mean_importance <- NULL
# All_Spp_meanAUCcv_Components <- setnames(All_Spp_meanAUCcv_Components, "3rd_Order_Neigh_Mean_Importance", "Third_Order_Neigh_Mean_Importance")

### Plot importance~distance
ggplot(All_Spp_meanAUCcv_Components, aes(x=Max_disp_dist_m, y = habAv_Mean_importance, color = Species, group = Species)) +
  geom_point(size = 3) + geom_line(size = 1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Mean importance of habitat availability")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))


### Plot AUC-cv~importance
ggplot(All_Spp_meanAUCcv_Components, aes(x=habAv_Mean_importance, y = Mean_AUC_cv, color = Species, group = Species)) +
  geom_point() + geom_line()



################################################################################
####### Get number of trees into main df for joint plotting ####################
#### Do invidual df's with only the mean values of each setting and species ####
################################################################################

### Hyarb ######################################################################

#default
Hyarb_defaultDispDist_meanntrees_tab <- data.frame(mean(Hyarb_DefaultDispDist_output_tab$ntrees))
Hyarb_defaultDispDist_meanntrees_tab <- setnames(Hyarb_defaultDispDist_meanntrees_tab, "mean.Hyarb_DefaultDispDist_output_tab.ntrees.",
                                                 "Mean_ntrees")
#Set identifier labels
Hyarb_defaultDispDist_meanntrees_tab$model = "Hyarb_defaultDispDist"

#300m
Hyarb_300mDispDist_meanntrees_tab <- data.frame(mean(Hyarb_300mDispDist_output_tab$ntrees))
Hyarb_300mDispDist_meanntrees_tab <- setnames(Hyarb_300mDispDist_meanntrees_tab, "mean.Hyarb_300mDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Hyarb_300mDispDist_meanntrees_tab$model = "Hyarb_300mDispDist"

#1Km
Hyarb_1KmDispDist_meanntrees_tab <- data.frame(mean(Hyarb_1KmDispDist_output_tab$ntrees))
Hyarb_1KmDispDist_meanntrees_tab <- setnames(Hyarb_1KmDispDist_meanntrees_tab, "mean.Hyarb_1KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Hyarb_1KmDispDist_meanntrees_tab$model = "Hyarb_1KmDispDist"

#2Km
Hyarb_2KmDispDist_meanntrees_tab <- data.frame(mean(Hyarb_2KmDispDist_output_tab$ntrees))
Hyarb_2KmDispDist_meanntrees_tab <- setnames(Hyarb_2KmDispDist_meanntrees_tab, "mean.Hyarb_2KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Hyarb_2KmDispDist_meanntrees_tab$model = "Hyarb_2KmDispDist"

#4Km
Hyarb_4KmDispDist_meanntrees_tab <- data.frame(mean(Hyarb_4KmDispDist_output_tab$ntrees))
Hyarb_4KmDispDist_meanntrees_tab <- setnames(Hyarb_4KmDispDist_meanntrees_tab, "mean.Hyarb_4KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Hyarb_4KmDispDist_meanntrees_tab$model = "Hyarb_4KmDispDist"

#6Km
Hyarb_6KmDispDist_meanntrees_tab <- data.frame(mean(Hyarb_6KmDispDist_output_tab$ntrees))
Hyarb_6KmDispDist_meanntrees_tab <- setnames(Hyarb_6KmDispDist_meanntrees_tab, "mean.Hyarb_6KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Hyarb_6KmDispDist_meanntrees_tab$model = "Hyarb_6KmDispDist"

#8Km
Hyarb_8KmDispDist_meanntrees_tab <- data.frame(mean(Hyarb_8KmDispDist_output_tab$ntrees))
Hyarb_8KmDispDist_meanntrees_tab <- setnames(Hyarb_8KmDispDist_meanntrees_tab, "mean.Hyarb_8KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Hyarb_8KmDispDist_meanntrees_tab$model = "Hyarb_8KmDispDist"

#10Km
Hyarb_10KmDispDist_meanntrees_tab <- data.frame(mean(Hyarb_10KmDispDist_output_tab$ntrees))
Hyarb_10KmDispDist_meanntrees_tab <- setnames(Hyarb_10KmDispDist_meanntrees_tab, "mean.Hyarb_10KmDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Hyarb_10KmDispDist_meanntrees_tab$model = "Hyarb_10KmDispDist"

##join individual model dfs
Hyarb_meanntrees_tab = rbind(Hyarb_300mDispDist_meanntrees_tab, Hyarb_1KmDispDist_meanntrees_tab, Hyarb_2KmDispDist_meanntrees_tab,
                             Hyarb_defaultDispDist_meanntrees_tab, Hyarb_4KmDispDist_meanntrees_tab, Hyarb_6KmDispDist_meanntrees_tab,
                             Hyarb_8KmDispDist_meanntrees_tab, Hyarb_10KmDispDist_meanntrees_tab)


#### Alobs #####################################################################

#default
Alobs_defaultDispDist_meanntrees_tab <- data.frame(mean(Alobs_DefaultDispDist_output_tab$ntrees))
Alobs_defaultDispDist_meanntrees_tab <- setnames(Alobs_defaultDispDist_meanntrees_tab, "mean.Alobs_DefaultDispDist_output_tab.ntrees.",
                                                 "Mean_ntrees")
#Set identifier labels
Alobs_defaultDispDist_meanntrees_tab$model = "Alobs_defaultDispDist"

#300m
Alobs_300mDispDist_meanntrees_tab <- data.frame(mean(Alobs_300mDispDist_output_tab$ntrees))
Alobs_300mDispDist_meanntrees_tab <- setnames(Alobs_300mDispDist_meanntrees_tab, "mean.Alobs_300mDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Alobs_300mDispDist_meanntrees_tab$model = "Alobs_300mDispDist"

#1Km
Alobs_1KmDispDist_meanntrees_tab <- data.frame(mean(Alobs_1KmDispDist_output_tab$ntrees))
Alobs_1KmDispDist_meanntrees_tab <- setnames(Alobs_1KmDispDist_meanntrees_tab, "mean.Alobs_1KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Alobs_1KmDispDist_meanntrees_tab$model = "Alobs_1KmDispDist"

#2Km
Alobs_2KmDispDist_meanntrees_tab <- data.frame(mean(Alobs_2KmDispDist_output_tab$ntrees))
Alobs_2KmDispDist_meanntrees_tab <- setnames(Alobs_2KmDispDist_meanntrees_tab, "mean.Alobs_2KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Alobs_2KmDispDist_meanntrees_tab$model = "Alobs_2KmDispDist"

#4Km
Alobs_4KmDispDist_meanntrees_tab <- data.frame(mean(Alobs_4KmDispDist_output_tab$ntrees))
Alobs_4KmDispDist_meanntrees_tab <- setnames(Alobs_4KmDispDist_meanntrees_tab, "mean.Alobs_4KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Alobs_4KmDispDist_meanntrees_tab$model = "Alobs_4KmDispDist"

#6Km
Alobs_6KmDispDist_meanntrees_tab <- data.frame(mean(Alobs_6KmDispDist_output_tab$ntrees))
Alobs_6KmDispDist_meanntrees_tab <- setnames(Alobs_6KmDispDist_meanntrees_tab, "mean.Alobs_6KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Alobs_6KmDispDist_meanntrees_tab$model = "Alobs_6KmDispDist"

#8Km
Alobs_8KmDispDist_meanntrees_tab <- data.frame(mean(Alobs_8KmDispDist_output_tab$ntrees))
Alobs_8KmDispDist_meanntrees_tab <- setnames(Alobs_8KmDispDist_meanntrees_tab, "mean.Alobs_8KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Alobs_8KmDispDist_meanntrees_tab$model = "Alobs_8KmDispDist"

#10Km
Alobs_10KmDispDist_meanntrees_tab <- data.frame(mean(Alobs_10KmDispDist_output_tab$ntrees))
Alobs_10KmDispDist_meanntrees_tab <- setnames(Alobs_10KmDispDist_meanntrees_tab, "mean.Alobs_10KmDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Alobs_10KmDispDist_meanntrees_tab$model = "Alobs_10KmDispDist"

##join individual model dfs
Alobs_meanntrees_tab = rbind(Alobs_300mDispDist_meanntrees_tab, Alobs_1KmDispDist_meanntrees_tab, Alobs_2KmDispDist_meanntrees_tab,
                             Alobs_defaultDispDist_meanntrees_tab, Alobs_4KmDispDist_meanntrees_tab, Alobs_6KmDispDist_meanntrees_tab,
                             Alobs_8KmDispDist_meanntrees_tab, Alobs_10KmDispDist_meanntrees_tab)


#### Bovar #####################################################################

#default
Bovar_defaultDispDist_meanntrees_tab <- data.frame(mean(Bovar_DefaultDispDist_output_tab$ntrees))
Bovar_defaultDispDist_meanntrees_tab <- setnames(Bovar_defaultDispDist_meanntrees_tab, "mean.Bovar_DefaultDispDist_output_tab.ntrees.",
                                                 "Mean_ntrees")
#Set identifier labels
Bovar_defaultDispDist_meanntrees_tab$model = "Bovar_defaultDispDist"

#300m
Bovar_300mDispDist_meanntrees_tab <- data.frame(mean(Bovar_300mDispDist_output_tab$ntrees))
Bovar_300mDispDist_meanntrees_tab <- setnames(Bovar_300mDispDist_meanntrees_tab, "mean.Bovar_300mDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Bovar_300mDispDist_meanntrees_tab$model = "Bovar_300mDispDist"

#1Km
Bovar_1KmDispDist_meanntrees_tab <- data.frame(mean(Bovar_1KmDispDist_output_tab$ntrees))
Bovar_1KmDispDist_meanntrees_tab <- setnames(Bovar_1KmDispDist_meanntrees_tab, "mean.Bovar_1KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Bovar_1KmDispDist_meanntrees_tab$model = "Bovar_1KmDispDist"

#2Km
Bovar_2KmDispDist_meanntrees_tab <- data.frame(mean(Bovar_2KmDispDist_output_tab$ntrees))
Bovar_2KmDispDist_meanntrees_tab <- setnames(Bovar_2KmDispDist_meanntrees_tab, "mean.Bovar_2KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Bovar_2KmDispDist_meanntrees_tab$model = "Bovar_2KmDispDist"

#4Km
# Bovar_4KmDispDist_meanntrees_tab <- data.frame(mean(Bovar_4KmDispDist_output_tab$ntrees))
# Bovar_4KmDispDist_meanntrees_tab <- setnames(Bovar_4KmDispDist_meanntrees_tab, "mean.Bovar_4KmDispDist_output_tab.ntrees.",
# "Mean_ntrees")
#Set identifier labels
# Bovar_4KmDispDist_meanntrees_tab$model = "Bovar_4KmDispDist"

#6Km
Bovar_6KmDispDist_meanntrees_tab <- data.frame(mean(Bovar_6KmDispDist_output_tab$ntrees))
Bovar_6KmDispDist_meanntrees_tab <- setnames(Bovar_6KmDispDist_meanntrees_tab, "mean.Bovar_6KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Bovar_6KmDispDist_meanntrees_tab$model = "Bovar_6KmDispDist"

#8Km
Bovar_8KmDispDist_meanntrees_tab <- data.frame(mean(Bovar_8KmDispDist_output_tab$ntrees))
Bovar_8KmDispDist_meanntrees_tab <- setnames(Bovar_8KmDispDist_meanntrees_tab, "mean.Bovar_8KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Bovar_8KmDispDist_meanntrees_tab$model = "Bovar_8KmDispDist"

#10Km
Bovar_10KmDispDist_meanntrees_tab <- data.frame(mean(Bovar_10KmDispDist_output_tab$ntrees))
Bovar_10KmDispDist_meanntrees_tab <- setnames(Bovar_10KmDispDist_meanntrees_tab, "mean.Bovar_10KmDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Bovar_10KmDispDist_meanntrees_tab$model = "Bovar_10KmDispDist"

##join individual model dfs
Bovar_meanntrees_tab = rbind(Bovar_300mDispDist_meanntrees_tab, Bovar_1KmDispDist_meanntrees_tab, Bovar_2KmDispDist_meanntrees_tab,
                             Bovar_defaultDispDist_meanntrees_tab, 
                             # Bovar_4KmDispDist_meanntrees_tab, 
                             Bovar_6KmDispDist_meanntrees_tab,
                             Bovar_8KmDispDist_meanntrees_tab, Bovar_10KmDispDist_meanntrees_tab)


#### Epcal #####################################################################

#default
Epcal_defaultDispDist_meanntrees_tab <- data.frame(mean(Epcal_DefaultDispDist_output_tab$ntrees))
Epcal_defaultDispDist_meanntrees_tab <- setnames(Epcal_defaultDispDist_meanntrees_tab, "mean.Epcal_DefaultDispDist_output_tab.ntrees.",
                                                 "Mean_ntrees")
#Set identifier labels
Epcal_defaultDispDist_meanntrees_tab$model = "Epcal_defaultDispDist"

#300m
Epcal_300mDispDist_meanntrees_tab <- data.frame(mean(Epcal_300mDispDist_output_tab$ntrees))
Epcal_300mDispDist_meanntrees_tab <- setnames(Epcal_300mDispDist_meanntrees_tab, "mean.Epcal_300mDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Epcal_300mDispDist_meanntrees_tab$model = "Epcal_300mDispDist"

#1Km
Epcal_1KmDispDist_meanntrees_tab <- data.frame(mean(Epcal_1KmDispDist_output_tab$ntrees))
Epcal_1KmDispDist_meanntrees_tab <- setnames(Epcal_1KmDispDist_meanntrees_tab, "mean.Epcal_1KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Epcal_1KmDispDist_meanntrees_tab$model = "Epcal_1KmDispDist"

#2Km
Epcal_2KmDispDist_meanntrees_tab <- data.frame(mean(Epcal_2KmDispDist_output_tab$ntrees))
Epcal_2KmDispDist_meanntrees_tab <- setnames(Epcal_2KmDispDist_meanntrees_tab, "mean.Epcal_2KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Epcal_2KmDispDist_meanntrees_tab$model = "Epcal_2KmDispDist"

#4Km
Epcal_4KmDispDist_meanntrees_tab <- data.frame(mean(Epcal_4KmDispDist_output_tab$ntrees))
Epcal_4KmDispDist_meanntrees_tab <- setnames(Epcal_4KmDispDist_meanntrees_tab, "mean.Epcal_4KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Epcal_4KmDispDist_meanntrees_tab$model = "Epcal_4KmDispDist"

#6Km
Epcal_6KmDispDist_meanntrees_tab <- data.frame(mean(Epcal_6KmDispDist_output_tab$ntrees))
Epcal_6KmDispDist_meanntrees_tab <- setnames(Epcal_6KmDispDist_meanntrees_tab, "mean.Epcal_6KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Epcal_6KmDispDist_meanntrees_tab$model = "Epcal_6KmDispDist"

#8Km
Epcal_8KmDispDist_meanntrees_tab <- data.frame(mean(Epcal_8KmDispDist_output_tab$ntrees))
Epcal_8KmDispDist_meanntrees_tab <- setnames(Epcal_8KmDispDist_meanntrees_tab, "mean.Epcal_8KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Epcal_8KmDispDist_meanntrees_tab$model = "Epcal_8KmDispDist"

#10Km
Epcal_10KmDispDist_meanntrees_tab <- data.frame(mean(Epcal_10KmDispDist_output_tab$ntrees))
Epcal_10KmDispDist_meanntrees_tab <- setnames(Epcal_10KmDispDist_meanntrees_tab, "mean.Epcal_10KmDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Epcal_10KmDispDist_meanntrees_tab$model = "Epcal_10KmDispDist"

##join individual model dfs
Epcal_meanntrees_tab = rbind(Epcal_300mDispDist_meanntrees_tab, Epcal_1KmDispDist_meanntrees_tab, Epcal_2KmDispDist_meanntrees_tab,
                             Epcal_defaultDispDist_meanntrees_tab, Epcal_4KmDispDist_meanntrees_tab, Epcal_6KmDispDist_meanntrees_tab,
                             Epcal_8KmDispDist_meanntrees_tab, Epcal_10KmDispDist_meanntrees_tab)


#### Peagg #####################################################################

#default
Peagg_defaultDispDist_meanntrees_tab <- data.frame(mean(Peagg_DefaultDispDist_output_tab$ntrees))
Peagg_defaultDispDist_meanntrees_tab <- setnames(Peagg_defaultDispDist_meanntrees_tab, "mean.Peagg_DefaultDispDist_output_tab.ntrees.",
                                                 "Mean_ntrees")
#Set identifier labels
Peagg_defaultDispDist_meanntrees_tab$model = "Peagg_defaultDispDist"

#300m
Peagg_300mDispDist_meanntrees_tab <- data.frame(mean(Peagg_300mDispDist_output_tab$ntrees))
Peagg_300mDispDist_meanntrees_tab <- setnames(Peagg_300mDispDist_meanntrees_tab, "mean.Peagg_300mDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Peagg_300mDispDist_meanntrees_tab$model = "Peagg_300mDispDist"

#1Km
Peagg_1KmDispDist_meanntrees_tab <- data.frame(mean(Peagg_1KmDispDist_output_tab$ntrees))
Peagg_1KmDispDist_meanntrees_tab <- setnames(Peagg_1KmDispDist_meanntrees_tab, "mean.Peagg_1KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Peagg_1KmDispDist_meanntrees_tab$model = "Peagg_1KmDispDist"

#2Km
Peagg_2KmDispDist_meanntrees_tab <- data.frame(mean(Peagg_2KmDispDist_output_tab$ntrees))
Peagg_2KmDispDist_meanntrees_tab <- setnames(Peagg_2KmDispDist_meanntrees_tab, "mean.Peagg_2KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Peagg_2KmDispDist_meanntrees_tab$model = "Peagg_2KmDispDist"

#4Km
Peagg_4KmDispDist_meanntrees_tab <- data.frame(mean(Peagg_4KmDispDist_output_tab$ntrees))
Peagg_4KmDispDist_meanntrees_tab <- setnames(Peagg_4KmDispDist_meanntrees_tab, "mean.Peagg_4KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Peagg_4KmDispDist_meanntrees_tab$model = "Peagg_4KmDispDist"

#6Km
Peagg_6KmDispDist_meanntrees_tab <- data.frame(mean(Peagg_6KmDispDist_output_tab$ntrees))
Peagg_6KmDispDist_meanntrees_tab <- setnames(Peagg_6KmDispDist_meanntrees_tab, "mean.Peagg_6KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Peagg_6KmDispDist_meanntrees_tab$model = "Peagg_6KmDispDist"

#8Km
Peagg_8KmDispDist_meanntrees_tab <- data.frame(mean(Peagg_8KmDispDist_output_tab$ntrees))
Peagg_8KmDispDist_meanntrees_tab <- setnames(Peagg_8KmDispDist_meanntrees_tab, "mean.Peagg_8KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Peagg_8KmDispDist_meanntrees_tab$model = "Peagg_8KmDispDist"

#10Km
Peagg_10KmDispDist_meanntrees_tab <- data.frame(mean(Peagg_10KmDispDist_output_tab$ntrees))
Peagg_10KmDispDist_meanntrees_tab <- setnames(Peagg_10KmDispDist_meanntrees_tab, "mean.Peagg_10KmDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Peagg_10KmDispDist_meanntrees_tab$model = "Peagg_10KmDispDist"

##join individual model dfs
Peagg_meanntrees_tab = rbind(Peagg_300mDispDist_meanntrees_tab, Peagg_1KmDispDist_meanntrees_tab, Peagg_2KmDispDist_meanntrees_tab,
                             Peagg_defaultDispDist_meanntrees_tab, Peagg_4KmDispDist_meanntrees_tab, Peagg_6KmDispDist_meanntrees_tab,
                             Peagg_8KmDispDist_meanntrees_tab, Peagg_10KmDispDist_meanntrees_tab)


#### Perid #####################################################################

#default
Perid_defaultDispDist_meanntrees_tab <- data.frame(mean(Perid_DefaultDispDist_output_tab$ntrees))
Perid_defaultDispDist_meanntrees_tab <- setnames(Perid_defaultDispDist_meanntrees_tab, "mean.Perid_DefaultDispDist_output_tab.ntrees.",
                                                 "Mean_ntrees")
#Set identifier labels
Perid_defaultDispDist_meanntrees_tab$model = "Perid_defaultDispDist"

#300m
Perid_300mDispDist_meanntrees_tab <- data.frame(mean(Perid_300mDispDist_output_tab$ntrees))
Perid_300mDispDist_meanntrees_tab <- setnames(Perid_300mDispDist_meanntrees_tab, "mean.Perid_300mDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Perid_300mDispDist_meanntrees_tab$model = "Perid_300mDispDist"

#1Km
Perid_1KmDispDist_meanntrees_tab <- data.frame(mean(Perid_1KmDispDist_output_tab$ntrees))
Perid_1KmDispDist_meanntrees_tab <- setnames(Perid_1KmDispDist_meanntrees_tab, "mean.Perid_1KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Perid_1KmDispDist_meanntrees_tab$model = "Perid_1KmDispDist"

#2Km
Perid_2KmDispDist_meanntrees_tab <- data.frame(mean(Perid_2KmDispDist_output_tab$ntrees))
Perid_2KmDispDist_meanntrees_tab <- setnames(Perid_2KmDispDist_meanntrees_tab, "mean.Perid_2KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Perid_2KmDispDist_meanntrees_tab$model = "Perid_2KmDispDist"

#4Km
Perid_4KmDispDist_meanntrees_tab <- data.frame(mean(Perid_4KmDispDist_output_tab$ntrees))
Perid_4KmDispDist_meanntrees_tab <- setnames(Perid_4KmDispDist_meanntrees_tab, "mean.Perid_4KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Perid_4KmDispDist_meanntrees_tab$model = "Perid_4KmDispDist"

#6Km
Perid_6KmDispDist_meanntrees_tab <- data.frame(mean(Perid_6KmDispDist_output_tab$ntrees))
Perid_6KmDispDist_meanntrees_tab <- setnames(Perid_6KmDispDist_meanntrees_tab, "mean.Perid_6KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Perid_6KmDispDist_meanntrees_tab$model = "Perid_6KmDispDist"

#8Km
Perid_8KmDispDist_meanntrees_tab <- data.frame(mean(Perid_8KmDispDist_output_tab$ntrees))
Perid_8KmDispDist_meanntrees_tab <- setnames(Perid_8KmDispDist_meanntrees_tab, "mean.Perid_8KmDispDist_output_tab.ntrees.",
                                             "Mean_ntrees")
#Set identifier labels
Perid_8KmDispDist_meanntrees_tab$model = "Perid_8KmDispDist"

#10Km
Perid_10KmDispDist_meanntrees_tab <- data.frame(mean(Perid_10KmDispDist_output_tab$ntrees))
Perid_10KmDispDist_meanntrees_tab <- setnames(Perid_10KmDispDist_meanntrees_tab, "mean.Perid_10KmDispDist_output_tab.ntrees.",
                                              "Mean_ntrees")
#Set identifier labels
Perid_10KmDispDist_meanntrees_tab$model = "Perid_10KmDispDist"

##join individual model dfs
Perid_meanntrees_tab = rbind(Perid_300mDispDist_meanntrees_tab, Perid_1KmDispDist_meanntrees_tab, Perid_2KmDispDist_meanntrees_tab,
                             Perid_defaultDispDist_meanntrees_tab, Perid_4KmDispDist_meanntrees_tab, Perid_6KmDispDist_meanntrees_tab,
                             Perid_8KmDispDist_meanntrees_tab, Perid_10KmDispDist_meanntrees_tab)


#### Join the AUC-cv dfs of all the species into one ###########################
All_Spp_meanntrees_tab = rbind(Alobs_meanntrees_tab, Bovar_meanntrees_tab, Epcal_meanntrees_tab,
                               Hyarb_meanntrees_tab, Peagg_meanntrees_tab, Perid_meanntrees_tab)

#### Join to main df with AUC cv & components ####
All_Spp_meanAUCcv_Components <- merge(All_Spp_meanAUCcv_Components, All_Spp_meanntrees_tab, by = "model")
head(All_Spp_meanAUCcv_Components)

#check if ntrees numeric, otherwise convert, to do linear model


### Plot AUC-cv ~ ntrees 

# All_Spp_meanAUCcv_Components$Max_Disp_Dist_factor <- as.factor(All_Spp_meanAUCcv_Components$Max_disp_dist_m)
# levels(All_Spp_meanAUCcv_Components$Max_Disp_Dist_factor)

ggplot(All_Spp_meanAUCcv_Components, aes(x= Mean_ntrees, y = Mean_AUC_cv , color = Species, group = Species)) +
  geom_smooth(method = "lm")+
  geom_point(size = 3)+ #geom_line(size = 1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Mean number of trees", y= "Mean cross-validated AUC")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))




################################################################################
#### Set dataframes per species including all the runs of each setting #########
#### All runs, not only means ##################################################

#### Get all the dfs into a list #### ###############################
all_dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
# names(dfs)

### Make sublist of default dfs
DefaultDispDist_output_list <- all_dfs[grep("DefaultDispDist_output_tab", names(all_dfs))]
names(DefaultDispDist_output_list)

### Remove the "X" column from default models
DefaultDispDist_output_list <- lapply(DefaultDispDist_output_list, function(x) x[names(x) != "X"])
list2env(DefaultDispDist_output_list, .GlobalEnv)

#### Get list of dfs by species ####
### List of all output_tabs
Output_tabs_list <- all_dfs[grep("output_tab", names(all_dfs))]
names(Output_tabs_list)

#remove repeat dfs (Bovar extras, because of 4Km being its default)
Output_tabs_list <- within(Output_tabs_list, rm(Bovar_DefaultDispDist_BRToutput_tab, Bovar_4KmDispDist_output_tab)) 

### remove 'X' column from default df's
Output_tabs_list <- lapply(Output_tabs_list, function(x) x[names(x) != "X"])
list2env(Output_tabs_list, .GlobalEnv)

#### Add extra column "model" to all the dfs in the Output_tabs_list
#When already with IDs (named list)
# Output_tabs_list <- mapply(cbind, Output_tabs_list, "model"=ID, SIMPLIFY=F)

### Function to turn name into element of new column
#version to pass through loop / lapply
dfnames_to_IDcol <- function(df, nm) {
  df$model <- nm
  return(df)
}


# dfnames_to_IDcol <- function(df) {
#   df$model <- deparse(substitute(df))
#   return(df)
# }

### Collect the names of the dfs into a table to work with them and name the elements of the list
setwd("M:/people/damiano/PhD/Sensitivity_DispDist/Mega_df_allmodels")
write.csv(names(Output_tabs_list), "Output_tabs_list.csv")

#Set the names of the new list corresponding to the original so they can be saved
dfs <- setNames(Output_tabs_list, names(Output_tabs_list))
#do the addition of the new column in the loop
dfs <- lapply(names(dfs), function(x) dfnames_to_IDcol(dfs[[x]], x))
#define again the names with the source df, as in new df the procedure made them being lost
dfs <- setNames(dfs, names(Output_tabs_list))
#store to the global environment the procedures done on the list
list2env(dfs, .GlobalEnv)


### remove 'X' column from default df's
dfs <- lapply(dfs, function(x) x[names(x) != "X"])
list2env(dfs, .GlobalEnv)

#doing it by hand works, but managed to get the loop to work
# Alobs_10KmDispDist_output_tab <- dfnames_to_IDcol(Alobs_10KmDispDist_output_tab)
# Perid_DefaultDispDist_output_tab <- dfnames_to_IDcol(Perid_DefaultDispDist_output_tab)

### Remove "_output_tab"  string part from model IDs
# library(stringr)
dfs <- lapply(dfs, function(x) x %>%
                mutate_at("model", str_replace, "_output_tab", ""))
list2env(dfs, .GlobalEnv)


#### Append the 'species', 'num_components, 'num_edges' and 'num patches'
### Default vs default (capital vs lowercase) discrepancy in dfs
#Do a copy
All_Spp_meanAUCcv_Components_2 <- All_Spp_meanAUCcv_Components
## change all the lowercase for uppercase in 'All_Spp_meanAUCcv_Components'
All_Spp_meanAUCcv_Components_2$model <- gsub('default', 'Default', All_Spp_meanAUCcv_Components_2$model) 

### join (merge) by model ID - Sets the appended columns contents to NA - Find other solution
#test
Perid_DefaultDispDist_output_tab_test <- merge(Perid_DefaultDispDist_output_tab_test, All_Spp_meanAUCcv_Components_2[ , c("Species", "Gender_Species", "Max_Disp_Dist", "Max_disp_dist_m", "Max_Disp_Dist_factor", "Number_of_patches", "Number_of_edges", "Number_of_components", "model")], by = "model", all.x = TRUE)
# merge(table1, table2[, c("pid", "val2")], by="pid")
### It works!
#another method
# table1$val2 <- table2$val2[match(table1$pid, table2$pid)]

### Get into the loop /lapply
dfs <- lapply(dfs, function(x) 
  merge(x, All_Spp_meanAUCcv_Components_2[ , c("Species", "Gender_Species", "Max_Disp_Dist", "Max_disp_dist_m", "Max_Disp_Dist_factor", "Number_of_patches", "Number_of_edges", "Number_of_components", "model")], by = "model", all.x = TRUE))
list2env(dfs, .GlobalEnv)


#### Do sublists by species ####
Alobs_Output_tabs_list <- dfs[grep("Alobs", names(dfs))]
names(Alobs_Output_tabs_list)
Bovar_Output_tabs_list <- dfs[grep("Bovar", names(dfs))]
Epcal_Output_tabs_list <- dfs[grep("Epcal", names(dfs))]
Hyarb_Output_tabs_list <- dfs[grep("Hyarb", names(dfs))]
Peagg_Output_tabs_list <- dfs[grep("Peagg", names(dfs))]
Perid_Output_tabs_list <- dfs[grep("Perid", names(dfs))]

# list2env(dfs, .GlobalEnv)

##### Get megadf's per species, containing the scores for all the models, not just the means
#### Bind the dfs into one for each species
#Use do.call to pass all the elements of the list and rbind them in one go
Alobs_values_tab <- do.call('rbind', Alobs_Output_tabs_list)
Bovar_values_tab <- do.call('rbind', Bovar_Output_tabs_list)
Epcal_values_tab <- do.call('rbind', Epcal_Output_tabs_list)
Hyarb_values_tab <- do.call('rbind', Hyarb_Output_tabs_list)
Peagg_values_tab <- do.call('rbind', Peagg_Output_tabs_list)
Perid_values_tab <- do.call('rbind', Perid_Output_tabs_list)

#### Do linear models of relations between variables (including all the models and not only means)
### AUC-cv ~ MaxDispDist
lm_AUCcv_MaxDispDist_Alobs <- lm(AUC_cv~Max_disp_dist_m, data = Alobs_values_tab)
lm_AUCcv_MaxDispDist_Bovar <- lm(AUC_cv~Max_disp_dist_m, data = Bovar_values_tab)
lm_AUCcv_MaxDispDist_Epcal <- lm(AUC_cv~Max_disp_dist_m, data = Epcal_values_tab)
lm_AUCcv_MaxDispDist_Hyarb <- lm(AUC_cv~Max_disp_dist_m, data = Hyarb_values_tab)
lm_AUCcv_MaxDispDist_Peagg <- lm(AUC_cv~Max_disp_dist_m, data = Peagg_values_tab)
lm_AUCcv_MaxDispDist_Perid <- lm(AUC_cv~Max_disp_dist_m, data = Perid_values_tab)

#### Do Linear models of AUC-cv ~ ntrees
lm_AUCcv_ntrees_Alobs <- lm(AUC_cv~ntrees, data = Alobs_values_tab)
lm_AUCcv_ntrees_Bovar <- lm(AUC_cv~ntrees, data = Bovar_values_tab)
lm_AUCcv_ntrees_Epcal <- lm(AUC_cv~ntrees, data = Epcal_values_tab)
lm_AUCcv_ntrees_Hyarb <- lm(AUC_cv~ntrees, data = Hyarb_values_tab)
lm_AUCcv_ntrees_Peagg <- lm(AUC_cv~ntrees, data = Peagg_values_tab)
lm_AUCcv_ntrees_Perid <- lm(AUC_cv~ntrees, data = Perid_values_tab)

### between AUC-cv ~ number of components
lm_AUCcv_Numcomponents_Alobs <- lm(AUC_cv~Number_of_components, data = Alobs_values_tab)
lm_AUCcv_Numcomponents_Bovar <- lm(AUC_cv~Number_of_components, data = Bovar_values_tab)
lm_AUCcv_Numcomponents_Epcal <- lm(AUC_cv~Number_of_components, data = Epcal_values_tab)
lm_AUCcv_Numcomponents_Hyarb <- lm(AUC_cv~Number_of_components, data = Hyarb_values_tab)
lm_AUCcv_Numcomponents_Peagg <- lm(AUC_cv~Number_of_components, data = Peagg_values_tab)
lm_AUCcv_Numcomponents_Perid <- lm(AUC_cv~Number_of_components, data = Perid_values_tab)

### between ntrees ~ max disp dist?

### plot again histogram of each species AUC-cv distr
hist(Alobs_values_tab$AUC_cv, xlab = "AUC-cv")
hist(Bovar_values_tab$AUC_cv, xlab = "AUC-cv")
hist(Epcal_values_tab$AUC_cv, xlab = "AUC-cv")
hist(Hyarb_values_tab$AUC_cv, xlab = "AUC-cv")
hist(Peagg_values_tab$AUC_cv, xlab = "AUC-cv")
hist(Perid_values_tab$AUC_cv, xlab = "AUC-cv")

#See the distribution by setting
#In two species
lapply(Alobs_Output_tabs_list, function(x) hist(x$AUC_cv))
lapply(Bovar_Output_tabs_list, function(x) hist(x$AUC_cv, main = x$model))
#distrs are quite varied

#### Get medians to plot relations ####
#AUC-cv
median_inList <- function(df) {
  df$median_AUCcv <- median(df$AUC_cv)
  return(df)
}

dfs <- lapply(dfs, function(df) median_inList(df))
list2env(dfs, .GlobalEnv)

#### Get mean AUC_cv too
mean_inList <- function(df) {
  df$mean_AUCcv <- mean(df$AUC_cv)
  return(df)
}

dfs <- lapply(dfs, function(df) mean_inList(df))
list2env(dfs, .GlobalEnv)


#### Same for number of trees
#median
median_inList_nt <- function(df) {
  df$median_ntrees <- median(df$ntrees)
  return(df)
}

dfs <- lapply(dfs, function(df) median_inList_nt(df))
list2env(dfs, .GlobalEnv)

#mean
mean_inList_nt <- function(df) {
  df$mean_ntrees <- mean(df$ntrees)
  return(df)
}

dfs <- lapply(dfs, function(df) mean_inList_nt(df))
list2env(dfs, .GlobalEnv)


#### Do again plots of AUCcv, ntrees, etc, with all the models

### Get all the dfs into the MacroDF for plotting
MacroDF_allSpp_DispDists <- do.call('rbind', dfs)


#Plot point to point
#Median  AUC-cv ~ MaxDispDist
ggplot(MacroDF_allSpp_DispDists, aes(x=Max_disp_dist_m, y = median_AUCcv, color = Species, group = Species)) +
  geom_point(size = 3) + geom_line(size = 1)+
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Median AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

# ggplot(All_Spp_meanAUCcv_tab, aes(x=Max_disp_dist_m, y = Mean_AUC_cv, color = Species, group = Species)) +
#   geom_point(size = 3) + geom_line(size = 1)+
#   theme_bw() +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
#   labs(x = "Maximum dispersal distance (m)", y= "Mean AUC-cv")+
#   theme(axis.text = element_text(size = 15))+
#   theme(axis.title = element_text(size = 20))+
#   theme(legend.text = element_text(size = 20))+
#   theme(legend.title = element_text(size = 20))
# scale_y_continuous(breaks = c(0.1:1))


#Plot w/ fitted trendline

#Mean  AUC-cv ~ MaxDispDist without #se -> std error -> 95% confidence interval

ggplot(MacroDF_allSpp_DispDists, aes(x=Max_disp_dist_m, y = mean_AUCcv, color = Species, group = Species)) +
  geom_smooth(method = "lm")+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Mean AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

#Median  AUC-cv ~ MaxDispDist
ggplot(MacroDF_allSpp_DispDists, aes(x=Max_disp_dist_m, y = median_AUCcv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "Median AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))


#scatterplot with ALL models (no median or mean)
ggplot(MacroDF_allSpp_DispDists, aes(x=Max_disp_dist_m, y = AUC_cv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Maximum dispersal distance (m)", y= "AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))


#### Plot again Num_components ~ num patches (with the df that includes all models) 
ggplot(MacroDF_allSpp_DispDists, aes(x=Number_of_patches, y = Number_of_components, color = Max_Disp_Dist_factor, group = Max_Disp_Dist_factor)) +
  geom_smooth(method = "lm", se = F)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Number of patches", y= "Number of components")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

#### Plot AUC-cv ~ num trees
### Mean w/o std error shade
ggplot(MacroDF_allSpp_DispDists, aes(x= mean_ntrees, y =  mean_AUCcv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Mean number of trees", y= "Mean AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

### Median
ggplot(MacroDF_allSpp_DispDists, aes(x= median_ntrees, y =  median_AUCcv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Median number of trees", y= "Median AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

### ALL
ggplot(MacroDF_allSpp_DispDists, aes(x= ntrees, y =  AUC_cv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Number of trees", y= "AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))


### check num trees ~ num patches
ggplot(MacroDF_allSpp_DispDists, aes(x= Number_of_patches, y =  mean_ntrees, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "Number of patches", y= "mean Number of trees")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

#### lm AUC-cv ~ num components

#### Plot AUC-cv ~ num components
### Mean AUC-cv w/o std error shade
ggplot(MacroDF_allSpp_DispDists, aes(x= log(Number_of_components), y =  mean_AUCcv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "log(Number of components)", y= "Mean AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

### Median AUC-cv
ggplot(MacroDF_allSpp_DispDists, aes(x= log(Number_of_components), y =  median_AUCcv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "log(Number of components)", y= "Median AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

## ALL models AUC-cv
ggplot(MacroDF_allSpp_DispDists, aes(x= log(Number_of_components), y =  AUC_cv, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "log(Number of components)", y= "AUC-cv")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))


#### Variable importance ####
### HSI
ggplot(MacroDF_allSpp_DispDists, aes(x= Max_disp_dist_m, y =  HSI_imp, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 80)) + 
  labs(x = "Max. disp. dist. (m)", y= "Importance of HSI")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

ggplot(MacroDF_allSpp_DispDists, aes(x= Max_disp_dist_m, y =  HSI_imp, color = Species, group = Species)) +
  # geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +  geom_line(size = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 80)) +
  labs(x = "Max. disp. dist. (m)", y= "Importance of HSI")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

###EgoSize (3rd order neighborhood)
ggplot(MacroDF_allSpp_DispDists, aes(x= Max_disp_dist_m, y =  EgoSize_imp, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 80)) +
  labs(x = "Max. disp. dist. (m)", y= "Importance of 3rd order neighborhood")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))

### HabAv 
ggplot(MacroDF_allSpp_DispDists, aes(x= Max_disp_dist_m, y =  habAv_imp, color = Species, group = Species)) +
  geom_smooth(method = "lm", se = FALSE)+
  geom_point(size = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  coord_cartesian(xlim = c(0, 10000), ylim = c(0, 80)) +
  labs(x = "Max. disp. dist. (m)", y= "Importance of habitat availability")+
  theme(axis.text = element_text(size = 15))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 20))+
  theme(legend.title = element_text(size = 20))


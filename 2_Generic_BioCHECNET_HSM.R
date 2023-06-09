###############################################################################################################
#### Script for Ensemble Habitat Suitability Modelling ########################################################
### With previous selection of variables and correlation checking #############################################
#### For "any" species ########################################################################################
#### Author: Damian O. Ortiz-Rodríguez, Antoine Guisan, Maarten J. van Strien #################################
#### Article: "Sensitivity of habitat network models to changes in maximum dispersal distance" ################
###############################################################################################################

library(sp)
library(raster)
library(rgdal)
library(biomod2)
library(tools)
library(usdm)
library(ecospat)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)
library(nnet)
library(forcats)
library(sna)
library(dismo)
library(ROCR)
library(ecodist)

##############################
#### Prepare inputs for HSM ####

#### Import rasters ####
setwd("E:/PhD/Automated_Run_BioCHECNET/GeoData/HSM_Predictors")
rasts = list.files(pattern = "\\.tif$")

#Test for unmatching properties like coord. ref. syst (CRS), extent, rowcol
rasts_var = vector()

for(rast in rasts)
{
  rstr = raster(rast)
  print(rast)
  print(crs(rstr))
  # assign(file_path_sans_ext(rast), rstr)
  # rasts_var= append(rasts_var,file_path_sans_ext(rast))
}

##### stack all predictors
PredictorStack <- raster::stack(rasts)

#Write R object to actual .grd raster stack
writeRaster(PredictorStack, filename=("PredictorStack2"), bandorder='BIL', suffix='numbers', overwrite = T)

#############################################
#When done once, for all the other species the previous step can be skipped, just load the saved stack

### Load the environmental predictor variables raster stack  ###

# setwd("C:/Users/damiano/Documents/PhD/Automated_Run_BioCHECNET/GeoData/HSM_Predictors")
PredictorStack <- stack('PredictorStack2.grd')
names(PredictorStack)


### load species data ###
setwd("C:/Users/damiano/Documents/PhD/Additional_species_runs/HSM") ### <- Actual workspace intended for most things

#needs 0's and 1's next to each coordinate pair
Gender_species <- read.csv("./Presences_tables/Gender_species_npr.csv")
head(Gender_species)

# resp.name, the name of studied species (Header of the binary presence column)
#_sa -> marker for specific study area
Ge_spe <- 'Gender_species_sa'
# The presence (or presence/absence) data, the resp.var
Ge_spePres <- as.numeric(Gender_species[,Ge_spe])
# the XY coordinates of species data
Ge_speXY <- Gender_species[,c("CX_HA_C","CY_HA_C")]


#####  Check for spatial autocorrelation
#Make an 'Ecospat-Biomod' compatible dataframe
Gender_species_redux <- data.frame(Gender_species$CX_HA_C, Gender_species$CY_HA_C, Gender_species$Gender_species_sa)

Gender_species_redux <- setnames(Gender_species_redux, "Gender_species.CX_HA_C", "CX_HA_C")
Gender_species_redux <- setnames(Gender_species_redux, "Gender_species.CY_HA_C", "CY_HA_C")
Gender_species_redux <- setnames(Gender_species_redux, "Gender_species.Gender_species_sa", "Gender_species_sa")

#extract
Gender_species_xy <- data.frame(Gender_species_redux$CX_HA_C, Gender_species_redux$CY_HA_C)

Gender_species_xy <- setnames(Gender_species_xy, "Gender_species_redux.CX_HA_C", "CX_HA_C")
Gender_species_xy <- setnames(Gender_species_xy, "Gender_species_redux.CY_HA_C", "CY_HA_C")


Gender_species_points <- extract(PredictorStack, Gender_species_xy, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE, 
                                      fun=NULL, na.rm=TRUE, 1, 20, df=TRUE, factors=FALSE)
#Add the X and Y to Gender_species_points
Gender_species_points$CX_HA_C <- Gender_species_xy$CX_HA_C
Gender_species_points$CY_HA_C <- Gender_species_xy$CY_HA_C

#Mantel correlogram

Gender_species_correlogram <- ecospat.mantel.correlogram (dfvar=Gender_species_points, colxy=22:23, n=100, colvar=2:21, max=1000, nclass=10, nperm=100)
#max, nclass, nperm: from Ecospat vignette


#############################
#### biomod2 modelling ######

#### Set data settings ####
Ge_spe_Data <- BIOMOD_FormatingData(resp.var = Ge_spePres,
                                    expl.var = PredictorStack,
                                    resp.xy = Ge_speXY,
                                    resp.name = Ge_spe,
                                    eval.resp.var=NULL, 
                                    eval.expl.var = NULL,
                                    eval.resp.xy = NULL,
                                    PA.nb.rep = 1,
                                    PA.nb.absences = 10000,
                                    PA.strategy = 'random',
                                    na.rm = TRUE)

#check if data are correctly formatted
Ge_spe_Data
plot(Ge_spe_Data)


#### Building models ####  
#keep default options
Ge_spe_ModOptions <-BIOMOD_ModelingOptions()

#In case Maxent doesn't work
#library("rJava")
#MAXENT.Phillips = list( path_to_maxent.jar = getwd())
#system.file("java", package="dismo")

### Modelling ###
# Computing the models
Ge_spe_ModelOut <- BIOMOD_Modeling(
  Ge_spe_Data,
  models = c('GLM','RF','MAXENT.Phillips'),
  models.options = Ge_spe_ModOptions,
  NbRunEval=3,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(Ge_spe,"Auto1",sep=""))


#When this step is over, have a look at some outputs :
# modeling summary
Ge_spe_ModelOut

#### Models evaluations ####

# get all models evaluation
Ge_spe_ModelEval <- get_evaluations(Ge_spe_ModelOut)
# print the dimnames of this object
dimnames(Ge_spe_ModelEval)

#print the ROC scores of all selected models
Ge_spe_ModelEval["ROC","Testing.data",,,]


#print relative importance of the explanatory variables
get_variables_importance(Ge_spe_ModelOut)

################################################

#### Ensemble modeling ####

Ge_spe_EnsembleMod1 <- BIOMOD_EnsembleModeling(
  modeling.output = Ge_spe_ModelOut,
  chosen.models = 'all', 
  em.by='all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )


# print summary
Ge_spe_EnsembleMod1

# get evaluation scores
get_evaluations(Ge_spe_EnsembleMod1)




#### 4 Projection ####

# projection over study area under current conditions
Ge_spe_Auto1_Proj <- BIOMOD_Projection(
  modeling.output = Ge_spe_ModelOut,
  new.env = PredictorStack ,
  proj.name = 'Ge_spe_Auto1st',
  selected.models = 'all',
  binary.meth = 'ROC',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of created object
Ge_spe_Auto1_Proj

# files created on hard drive
list.files("Gender.species.sa/proj_Ge_spe_Auto1st/")


# make some plots sub-selected by str.grep argument
plot(Ge_spe_Auto1_Proj, str.grep = 'MAXENT.Phillips') 
plot(Ge_spe_Auto1_Proj, str.grep = 'GLM')
plot(Ge_spe_Auto1_Proj, str.grep = 'RF')


#get the projected map

Ge_spe_Auto1_Projectedmap <- get_predictions(Ge_spe_Auto1_Proj)
Ge_spe_Auto1_Projectedmap

plot(Ge_spe_Auto1_Projectedmap)

#Visualize maps
getwd()
Auto1st_Gespe_HSM <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.grd")
#It's a raster with multiple bands, one for each projection made by each model run
plot(Auto1st_Gespe_HSM)

#To check the names and order of the layers in the stack
names(Auto1st_Gespe_HSM)

#Get individual layers to check if they are all proper maps
#load the individual bands to a raster object
#use .gri, not .grd

Auto1st_Gespe_HSM_Band1 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 1)
Auto1st_Gespe_HSM_Band2 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 2)
Auto1st_Gespe_HSM_Band3 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 3)
Auto1st_Gespe_HSM_Band4 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 4)
Auto1st_Gespe_HSM_Band5 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 5)
Auto1st_Gespe_HSM_Band6 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 6)
Auto1st_Gespe_HSM_Band7 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 7)
Auto1st_Gespe_HSM_Band8 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 8)
Auto1st_Gespe_HSM_Band9 = raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa.gri", band = 9)

plot(Auto1st_Gespe_HSM_Band1, main= "PA1_RUN1_GLM")
plot(Auto1st_Gespe_HSM_Band2, main= "PA1_RUN1_RF")
plot(Auto1st_Gespe_HSM_Band3, main= "PA1_RUN1_MAXENT.Phillips")
plot(Auto1st_Gespe_HSM_Band4, main= "PA1_RUN2_GLM")
plot(Auto1st_Gespe_HSM_Band5, main= "PA1_RUN2_RF")
plot(Auto1st_Gespe_HSM_Band6, main= "PA1_RUN2_MAXENT.Phillips")
plot(Auto1st_Gespe_HSM_Band7, main= "PA1_RUN3_GLM")
plot(Auto1st_Gespe_HSM_Band8, main= "PA1_RUN3_RF")
plot(Auto1st_Gespe_HSM_Band9, main= "PA1_RUN3_MAXENT.Phillips")

# #Same as above but binarized
ROCBin_Auto1st_Gespe_HSM <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa_ROCbin.grd")
plot(ROCBin_Auto1st_Gespe_HSM)

# #### Convert out.rasters to GeoTIFF to be able to visualize in ArcGIS ####
writeRaster(Auto1st_Gespe_HSM, filename="Auto1st_Gespe_HSM_stack", format="GTiff")
writeRaster(ROCBin_Auto1st_Gespe_HSM, filename="ROCBin_Auto1st_Gespe_HSM_stack", format="GTiff")


#### 5 Ensemble_Forecasting ####
Ge_spe_Auto1_EF <- BIOMOD_EnsembleForecasting(
  EM.output = Ge_spe_EnsembleMod1,
  projection.output = Ge_spe_Auto1_Proj,
  binary.meth="ROC") 

Ge_spe_Auto1_EF
plot(Ge_spe_Auto1_EF)

#Visualize ensemble maps
#(binary transformation of single models based on ROC optimized threshold)

#mean by ROC = mean ensemble model
Auto1st_Gespe_mean_EM <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/individual_projections/Gender.species.sa_EMmeanByROC_mergedAlgo_mergedRun_mergedData.grd")
#ca by ROC = committee averaging ensemble model (binary transformation of single models based on ROC optimized threshold 
Auto1st_Gespe_ca_EM <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/individual_projections/Gender.species.sa_EMcaByROC_mergedAlgo_mergedRun_mergedData.grd")
#cv by ROC = coefficient of variation - lower score, better model
Auto1st_Gespe_cv_EM <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/individual_projections/Gender.species.sa_EMcvByROC_mergedAlgo_mergedRun_mergedData.grd")
#median by ROC = median ensemble model
Auto1st_Gespe_median_EM <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/individual_projections/Gender.species.sa_EMmedianByROC_mergedAlgo_mergedRun_mergedData.grd")
#wmean by ROC = weighted mean ensemble models (weighted by single models ROC score)
Auto1st_Gespe_wmean_EM <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/individual_projections/Gender.species.sa_EMwmeanByROC_mergedAlgo_mergedRun_mergedData.grd")

#This one is a stack of all the previous 5 plus the ciInf and ciSup   
  #ciInf = Confidence interval Inferior, ciSup = Confidence interval Superior
Auto1st_Gespe_ensemble <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa_ensemble.grd")
#ROC Binarization
ROCBin_Auto1st_Gespe_ensemble <-raster("./Gender.species.sa/proj_Ge_spe_Auto1st/proj_Ge_spe_Auto1st_Gender.species.sa_ensemble_ROCbin.grd")

#Plot to check visually bands
plot(Auto1st_Gespe_ca_EM, main= "comittee_averaging")
plot(Auto1st_Gespe_mean_EM, main= "mean ensemble model")
plot(Auto1st_Gespe_median_EM, main= "median  ensemble model")
plot(Auto1st_Gespe_wmean_EM, main= "weighted mean ensemble model")
plot(Auto1st_Gespe_cv_EM, main= "coefficient of variation")

plot(Auto1st_Gespe_ensemble, main= "ensemble") #band 1 is 'mean'
plot(ROCBin_Auto1st_Gespe_ensemble, main= "ROC Bin")


# #### Convert out.rasters to GeoTIFF to be able to visualize in ArcGIS ####
#band 1 is 'mean'
writeRaster(Auto1st_Gespe_ensemble, filename="Auto1st_Gespe_ensemble_1", band = 1, format="GTiff")
writeRaster(ROCBin_Auto1st_Gespe_ensemble, filename="ROCBin_Auto1st_Gespe_ensemble_1", band = 1, format="GTiff")


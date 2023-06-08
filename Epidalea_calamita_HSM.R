###############################################################################################################
#### Script for Ensemble Habitat Suitability Modelling ########################################################
### With previous selection of variables and correlation checking #############################################
#### For Epidalea calamita ####################################################################################
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

setwd("C:/Users/damiano/Documents/PhD/Automated_Run_BioCHECNET/GeoData/HSM_Predictors")

### load the environmental variables raster stack 
#PredictorStack 
PredictorStack <- stack('PredictorStack2.grd')
names(PredictorStack)

### load species data ###
setwd("C:/Users/damiano/Documents/PhD/Additional_species_runs/HSM") ### <- Actual workspace intended for most things

#needs 0's and 1's next to each coordinate pair
#Gender_species <- read.csv("./Presences_tables/Gender_species_npr.csv")
Epidalea_calamita <- read.csv("./Presences_tables/Epidalea_calamita_npr.csv")
head(Epidalea_calamita)

# resp.name, the name of studied species (Header of the binary presence column)
#_sa -> marker for specific study area
Ep_cal <- 'Epidalea_calamita_sa'
# The presence (or presence/absence) data, the resp.var
Ep_calPres <- as.numeric(Epidalea_calamita[,Ep_cal])
# the XY coordinates of species data
Ep_calXY <- Epidalea_calamita[,c("CX_HA_C","CY_HA_C")]


#####  Check for spatial autocorrelation
#Make an 'Ecospat-Biomod' compatible dataframe
Epidalea_calamita_redux <- data.frame(Epidalea_calamita$CX_HA_C, Epidalea_calamita$CY_HA_C, Epidalea_calamita$Epidalea_calamita_sa)

Epidalea_calamita_redux <- setnames(Epidalea_calamita_redux, "Epidalea_calamita.CX_HA_C", "CX_HA_C")
Epidalea_calamita_redux <- setnames(Epidalea_calamita_redux, "Epidalea_calamita.CY_HA_C", "CY_HA_C")
Epidalea_calamita_redux <- setnames(Epidalea_calamita_redux, "Epidalea_calamita.Epidalea_calamita_sa", "Epidalea_calamita_sa")

# extract
Epidalea_calamita_xy <- data.frame(Epidalea_calamita_redux$CX_HA_C, Epidalea_calamita_redux$CY_HA_C)

Epidalea_calamita_xy <- setnames(Epidalea_calamita_xy, "Epidalea_calamita_redux.CX_HA_C", "CX_HA_C")
Epidalea_calamita_xy <- setnames(Epidalea_calamita_xy, "Epidalea_calamita_redux.CY_HA_C", "CY_HA_C")


Epidalea_calamita_points <- extract(PredictorStack, Epidalea_calamita_xy, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE, 
                                 fun=NULL, na.rm=TRUE, 1, 20, df=TRUE, factors=FALSE)
#Add the X and Y to Gender_species_points
Epidalea_calamita_points$CX_HA_C <- Epidalea_calamita_xy$CX_HA_C
Epidalea_calamita_points$CY_HA_C <- Epidalea_calamita_xy$CY_HA_C

#Mantel correlogram

Epidalea_calamita_correlogram <- ecospat.mantel.correlogram (dfvar=Epidalea_calamita_points, colxy=22:23, n=100, colvar=2:21, max=1000, nclass=10, nperm=100)
#max, nclass, nperm: from Ecospat vignette


#############################
#### biomod2 modelling ######

#### Set data settings ####
Ep_cal_Data <- BIOMOD_FormatingData(resp.var = Ep_calPres,
                                    expl.var = PredictorStack,
                                    resp.xy = Ep_calXY,
                                    resp.name = Ep_cal,
                                    eval.resp.var=NULL, 
                                    eval.expl.var = NULL,
                                    eval.resp.xy = NULL,
                                    PA.nb.rep = 1,
                                    PA.nb.absences = 10000,
                                    PA.strategy = 'random',
                                    na.rm = TRUE)

#check if data are correctly formatted
Ep_cal_Data
plot(Ep_cal_Data)

#### Building models ####  
#keep default options
Ep_cal_ModOptions <-BIOMOD_ModelingOptions()

#In case Maxent doesn't work
#library("rJava")
#MAXENT.Phillips = list( path_to_maxent.jar = getwd())
#system.file("java", package="dismo")

### Modelling ###
# Computing the models
Ep_cal_ModelOut <- BIOMOD_Modeling(
  Ep_cal_Data,
  models = c('GLM','RF','MAXENT.Phillips'),
  models.options = Ep_cal_ModOptions,
  NbRunEval=3,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(Ep_cal,"Auto1",sep=""))


#When this step is over, have a look at some outputs :
# modeling summary
Ep_cal_ModelOut

#### Models evaluations ####

# get all models evaluation
Ep_cal_ModelEval <- get_evaluations(Ep_cal_ModelOut)
# print the dimnames of this object
dimnames(Ep_cal_ModelEval)

#print the ROC scores of all selected models
Ep_cal_ModelEval["ROC","Testing.data",,,]


#print relative importance of the explanatory variables
get_variables_importance(Ep_cal_ModelOut)

################################################

#### Ensemble modeling ####

Ep_cal_EnsembleMod1 <- BIOMOD_EnsembleModeling(
  modeling.output = Ep_cal_ModelOut,
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
Ep_cal_EnsembleMod1

# get evaluation scores
get_evaluations(Ep_cal_EnsembleMod1)




#### 4 Projection ####

# projection over study area under current conditions
Ep_cal_Auto1_Proj <- BIOMOD_Projection(
  modeling.output = Ep_cal_ModelOut,
  new.env = PredictorStack ,
  proj.name = 'Ep_cal_Auto1st',
  selected.models = 'all',
  binary.meth = 'ROC',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of created object
Ep_cal_Auto1_Proj

# files created on hard drive
list.files("Epidalea.calamita.sa/proj_Ep_cal_Auto1st/")



# make some plots sub-selected by str.grep argument
plot(Ep_cal_Auto1_Proj, str.grep = 'MAXENT.Phillips') #Change method
plot(Ep_cal_Auto1_Proj, str.grep = 'GLM')
plot(Ep_cal_Auto1_Proj, str.grep = 'RF')


#get the projected map

Ep_cal_Auto1_Projectedmap <- get_predictions(Ep_cal_Auto1_Proj)
Ep_cal_Auto1_Projectedmap

plot(Ep_cal_Auto1_Projectedmap)

#Visualize maps
getwd()
Auto1st_Epcal_HSM <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.grd")
#It's a raster with multiple bands, one for each projection made by each model run
plot(Auto1st_Epcal_HSM)

#To check the names and order of the layers in the stack
names(Auto1st_Epcal_HSM)

#Get individual layers to check if they are all proper maps
#load the individual bands to a raster object
#use .gri, not .grd

Auto1st_Epcal_HSM_Band1 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 1)
Auto1st_Epcal_HSM_Band2 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 2)
Auto1st_Epcal_HSM_Band3 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 3)
Auto1st_Epcal_HSM_Band4 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 4)
Auto1st_Epcal_HSM_Band5 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 5)
Auto1st_Epcal_HSM_Band6 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 6)
Auto1st_Epcal_HSM_Band7 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 7)
Auto1st_Epcal_HSM_Band8 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 8)
Auto1st_Epcal_HSM_Band9 = raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa.gri", band = 9)

plot(Auto1st_Epcal_HSM_Band1, main= "PA1_RUN1_GLM")
plot(Auto1st_Epcal_HSM_Band2, main= "PA1_RUN1_RF")
plot(Auto1st_Epcal_HSM_Band3, main= "PA1_RUN1_MAXENT.Phillips")
plot(Auto1st_Epcal_HSM_Band4, main= "PA1_RUN2_GLM")
plot(Auto1st_Epcal_HSM_Band5, main= "PA1_RUN2_RF")
plot(Auto1st_Epcal_HSM_Band6, main= "PA1_RUN2_MAXENT.Phillips")
plot(Auto1st_Epcal_HSM_Band7, main= "PA1_RUN3_GLM")
plot(Auto1st_Epcal_HSM_Band8, main= "PA1_RUN3_RF")
plot(Auto1st_Epcal_HSM_Band9, main= "PA1_RUN3_MAXENT.Phillips")

# #Same as above but binarized
ROCBin_Auto1st_Epcal_HSM <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa_ROCbin.grd")
plot(ROCBin_Auto1st_Epcal_HSM)

# #### Convert out.rasters to GeoTIFF to be able to visualize in ArcGIS ####
writeRaster(Auto1st_Epcal_HSM, filename="Auto1st_Epcal_HSM_stack", format="GTiff")
writeRaster(ROCBin_Auto1st_Epcal_HSM, filename="ROCBin_Auto1st_Epcal_HSM_stack", format="GTiff")


#### 5 Ensemble_Forecasting ####
Ep_cal_Auto1_EF <- BIOMOD_EnsembleForecasting(
  EM.output = Ep_cal_EnsembleMod1,
  projection.output = Ep_cal_Auto1_Proj,
  binary.meth="ROC")

Ep_cal_Auto1_EF
plot(Ep_cal_Auto1_EF)

#Visualize ensemble maps
#(binary transformation of single models based on ROC optimized threshold)

#mean by ROC = mean ensemble model
Auto1st_Epcal_mean_EM <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/individual_projections/Epidalea.calamita.sa_EMmeanByROC_mergedAlgo_mergedRun_mergedData.grd")
#ca by ROC = committee averaging ensemble model (binary transformation of single models based on ROC optimized threshold 
Auto1st_Epcal_ca_EM <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/individual_projections/Epidalea.calamita.sa_EMcaByROC_mergedAlgo_mergedRun_mergedData.grd")
#cv by ROC = coefficient of variation - lower score, better model
Auto1st_Epcal_cv_EM <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/individual_projections/Epidalea.calamita.sa_EMcvByROC_mergedAlgo_mergedRun_mergedData.grd")
#median by ROC = median ensemble model
Auto1st_Epcal_median_EM <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/individual_projections/Epidalea.calamita.sa_EMmedianByROC_mergedAlgo_mergedRun_mergedData.grd")
#wmean by ROC = weighted mean ensemble models (weighted by single models ROC score)
Auto1st_Epcal_wmean_EM <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/individual_projections/Epidalea.calamita.sa_EMwmeanByROC_mergedAlgo_mergedRun_mergedData.grd")

#This one is a stack of all the previous 5 plus the ciInf and ciSup   
  #ciInf = Confidence interval Inferior, ciSup = Confidence interval Superior
Auto1st_Epcal_ensemble <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa_ensemble.grd")
#ROC Binarization
ROCBin_Auto1st_Epcal_ensemble <-raster("./Epidalea.calamita.sa/proj_Ep_cal_Auto1st/proj_Ep_cal_Auto1st_Epidalea.calamita.sa_ensemble_ROCbin.grd")

#Plot to check visually bands
plot(Auto1st_Epcal_ca_EM, main= "comittee_averaging")
plot(Auto1st_Epcal_mean_EM, main= "mean ensemble model")
plot(Auto1st_Epcal_median_EM, main= "median  ensemble model")
plot(Auto1st_Epcal_wmean_EM, main= "weighted mean ensemble model")
plot(Auto1st_Epcal_cv_EM, main= "coefficient of variation")

plot(Auto1st_Epcal_ensemble, main= "ensemble") #band 1 is 'mean'
plot(ROCBin_Auto1st_Epcal_ensemble, main= "ROC Bin")


# #### Convert out.rasters to GeoTIFF to be able to visualize in ArcGIS ####
#band 1 is 'mean'
writeRaster(Auto1st_Epcal_ensemble, filename="Auto1st_Epcal_ensemble_1", band = 1, format="GTiff") 
writeRaster(ROCBin_Auto1st_Epcal_ensemble, filename="ROCBin_Auto1st_Epcal_ensemble_1", band = 1, format="GTiff") 

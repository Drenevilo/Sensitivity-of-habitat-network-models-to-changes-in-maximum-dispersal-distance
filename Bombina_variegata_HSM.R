###############################################################################################################
#### Script for Ensemble Habitat Suitability Modelling ########################################################
### With previous selection of variables and correlation checking #############################################
#### For Bombina_variegata ####################################################################################
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
Bombina_variegata <- read.csv("./Presences_tables/Bombina_variegata_npr.csv")
head(Bombina_variegata)

# resp.name, the name of studied species (Header of the binary presence column)
#_sa -> marker for specific study area
Bo_var <- 'Bombina_variegata_sa'
# The presence (or presence/absence) data, the resp.var
Bo_varPres <- as.numeric(Bombina_variegata[,Bo_var])
# the XY coordinates of species data
Bo_varXY <- Bombina_variegata[,c("CX_HA_C","CY_HA_C")]


#####  Check for spatial autocorrelation
#Make an 'Ecospat-Biomod' compatible dataframe
Bombina_variegata_redux <- data.frame(Bombina_variegata$CX_HA_C, Bombina_variegata$CY_HA_C, Bombina_variegata$Bombina_variegata_sa)

Bombina_variegata_redux <- setnames(Bombina_variegata_redux, "Bombina_variegata.CX_HA_C", "CX_HA_C")
Bombina_variegata_redux <- setnames(Bombina_variegata_redux, "Bombina_variegata.CY_HA_C", "CY_HA_C")
Bombina_variegata_redux <- setnames(Bombina_variegata_redux, "Bombina_variegata.Bombina_variegata_sa", "Bombina_variegata_sa")

# extract
Bombina_variegata_xy <- data.frame(Bombina_variegata_redux$CX_HA_C, Bombina_variegata_redux$CY_HA_C)

Bombina_variegata_xy <- setnames(Bombina_variegata_xy, "Bombina_variegata_redux.CX_HA_C", "CX_HA_C")
Bombina_variegata_xy <- setnames(Bombina_variegata_xy, "Bombina_variegata_redux.CY_HA_C", "CY_HA_C")


Bombina_variegata_points <- extract(PredictorStack, Bombina_variegata_xy, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE, 
                                 fun=NULL, na.rm=TRUE, 1, 20, df=TRUE, factors=FALSE)
#Add the X and Y to Gender_species_points
Bombina_variegata_points$CX_HA_C <- Bombina_variegata_xy$CX_HA_C
Bombina_variegata_points$CY_HA_C <- Bombina_variegata_xy$CY_HA_C

#Mantel correlogram

Bombina_variegata_correlogram <- ecospat.mantel.correlogram (dfvar=Bombina_variegata_points, colxy=22:23, n=100, colvar=2:21, max=1000, nclass=10, nperm=100)
#max, nclass, nperm: from Ecospat vignette


#############################
#### biomod2 modelling ######

## Follow same settings as for Hyla arborea

#### Set data settings ####
Bo_var_Data <- BIOMOD_FormatingData(resp.var = Bo_varPres,
                                    expl.var = PredictorStack,
                                    resp.xy = Bo_varXY,
                                    resp.name = Bo_var,
                                    eval.resp.var=NULL, 
                                    eval.expl.var = NULL,
                                    eval.resp.xy = NULL,
                                    PA.nb.rep = 1,
                                    PA.nb.absences = 10000,
                                    PA.strategy = 'random',
                                    na.rm = TRUE)

#check if data are correctly formatted
Bo_var_Data
plot(Bo_var_Data)

#### Building models ####  
#keep default options
Bo_var_ModOptions <-BIOMOD_ModelingOptions()

#In case Maxent doesn't work
#library("rJava")
#MAXENT.Phillips = list( path_to_maxent.jar = getwd())
#system.file("java", package="dismo")

### Modelling ###
# Computing the models
Bo_var_ModelOut <- BIOMOD_Modeling(
  Bo_var_Data,
  models = c('GLM','RF','MAXENT.Phillips'),
  models.options = Bo_var_ModOptions,
  NbRunEval=3,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(Bo_var,"Auto1",sep=""))


#When this step is over, have a look at some outputs :
# modeling summary
Bo_var_ModelOut

#### Models evaluations ####

# get all models evaluation
Bo_var_ModelEval <- get_evaluations(Bo_var_ModelOut)
# print the dimnames of this object
dimnames(Bo_var_ModelEval)

#print the ROC scores of all selected models
Bo_var_ModelEval["ROC","Testing.data",,,]


#print relative importance of the explanatory variables
get_variables_importance(Bo_var_ModelOut)

################################################

#### Ensemble modeling ####

Bo_var_EnsembleMod1 <- BIOMOD_EnsembleModeling(
  modeling.output = Bo_var_ModelOut,
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
Bo_var_EnsembleMod1

# get evaluation scores
get_evaluations(Bo_var_EnsembleMod1)




#### 4 Projection ####

# projection over study area under current conditions
Bo_var_Auto1_Proj <- BIOMOD_Projection(
  modeling.output = Bo_var_ModelOut,
  new.env = PredictorStack ,
  proj.name = 'Bo_var_Auto1st',
  selected.models = 'all',
  binary.meth = 'ROC',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of created oject
Bo_var_Auto1_Proj

# files created on hard drive
list.files("Bombina.variegata.sa/proj_Bo_var_Auto1st/")



# make some plots sub-selected by str.grep argument
plot(Bo_var_Auto1_Proj, str.grep = 'MAXENT.Phillips')
plot(Bo_var_Auto1_Proj, str.grep = 'GLM')
plot(Bo_var_Auto1_Proj, str.grep = 'RF')


#get the projected map

Bo_var_Auto1_Projectedmap <- get_predictions(Bo_var_Auto1_Proj)
Bo_var_Auto1_Projectedmap

plot(Bo_var_Auto1_Projectedmap)

#Visualize maps
getwd()
Auto1st_Bovar_HSM <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.grd")
#It's a raster with multiple bands, one for each projection made by each model run
plot(Auto1st_Bovar_HSM)

#To check the names and order of the layers in the stack
names(Auto1st_Bovar_HSM)

#Get individual layers to check if they are all proper maps
#load the individual bands to a raster object
#use .gri, not .grd!

Auto1st_Bovar_HSM_Band1 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 1)
Auto1st_Bovar_HSM_Band2 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 2)
Auto1st_Bovar_HSM_Band3 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 3)
Auto1st_Bovar_HSM_Band4 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 4)
Auto1st_Bovar_HSM_Band5 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 5)
Auto1st_Bovar_HSM_Band6 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 6)
Auto1st_Bovar_HSM_Band7 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 7)
Auto1st_Bovar_HSM_Band8 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 8)
Auto1st_Bovar_HSM_Band9 = raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa.gri", band = 9)

plot(Auto1st_Bovar_HSM_Band1, main= "PA1_RUN1_GLM")
plot(Auto1st_Bovar_HSM_Band2, main= "PA1_RUN1_RF")
plot(Auto1st_Bovar_HSM_Band3, main= "PA1_RUN1_MAXENT.Phillips")
plot(Auto1st_Bovar_HSM_Band4, main= "PA1_RUN2_GLM")
plot(Auto1st_Bovar_HSM_Band5, main= "PA1_RUN2_RF")
plot(Auto1st_Bovar_HSM_Band6, main= "PA1_RUN2_MAXENT.Phillips")
plot(Auto1st_Bovar_HSM_Band7, main= "PA1_RUN3_GLM")
plot(Auto1st_Bovar_HSM_Band8, main= "PA1_RUN3_RF")
plot(Auto1st_Bovar_HSM_Band9, main= "PA1_RUN3_MAXENT.Phillips")

# #Same as above but binarized
ROCBin_Auto1st_Bovar_HSM <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa_ROCbin.grd")
plot(ROCBin_Auto1st_Bovar_HSM)

# #### Convert out.rasters to GeoTIFF to be able to visualize in ArcGIS ####
writeRaster(Auto1st_Bovar_HSM, filename="Auto1st_Bovar_HSM_stack", format="GTiff")
writeRaster(ROCBin_Auto1st_Bovar_HSM, filename="ROCBin_Auto1st_Bovar_HSM_stack", format="GTiff")



#### 5 Ensemble_Forecasting ####
Bo_var_Auto1_EF <- BIOMOD_EnsembleForecasting(
  EM.output = Bo_var_EnsembleMod1,
  projection.output = Bo_var_Auto1_Proj,
  binary.meth="ROC")

Bo_var_Auto1_EF
plot(Bo_var_Auto1_EF)

#Visualize ensemble maps
#(binary transformation of single models based on ROC optimized threshold)

#mean by ROC = mean ensemble model
Auto1st_Bovar_mean_EM <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/individual_projections/Bombina.variegata.sa_EMmeanByROC_mergedAlgo_mergedRun_mergedData.grd")
#ca by ROC = committee averaging ensemble model (binary transformation of single models based on ROC optimized threshold 
Auto1st_Bovar_ca_EM <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/individual_projections/Bombina.variegata.sa_EMcaByROC_mergedAlgo_mergedRun_mergedData.grd")
#cv by ROC = coefficient of variation - lower score, better model
Auto1st_Bovar_cv_EM <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/individual_projections/Bombina.variegata.sa_EMcvByROC_mergedAlgo_mergedRun_mergedData.grd")
#median by ROC = median ensemble model
Auto1st_Bovar_median_EM <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/individual_projections/Bombina.variegata.sa_EMmedianByROC_mergedAlgo_mergedRun_mergedData.grd")
#wmean by ROC = weighted mean ensemble models (weighted by single models ROC score)
Auto1st_Bovar_wmean_EM <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/individual_projections/Bombina.variegata.sa_EMwmeanByROC_mergedAlgo_mergedRun_mergedData.grd")

#This is a stack of all the previous 5 plus the ciInf and ciSup
  #ciInf -> Confidence interval Inferior, ciSup = Confidence interval Superior
Auto1st_Bovar_ensemble <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa_ensemble.grd")
#ROC Binarization
ROCBin_Auto1st_Bovar_ensemble <-raster("./Bombina.variegata.sa/proj_Bo_var_Auto1st/proj_Bo_var_Auto1st_Bombina.variegata.sa_ensemble_ROCbin.grd")

#Plot to check visually bands
plot(Auto1st_Bovar_ca_EM, main= "comittee_averaging")
plot(Auto1st_Bovar_mean_EM, main= "mean ensemble model")
plot(Auto1st_Bovar_median_EM, main= "median  ensemble model")
plot(Auto1st_Bovar_wmean_EM, main= "weighted mean ensemble model")
plot(Auto1st_Bovar_cv_EM, main= "coefficient of variation")

plot(Auto1st_Bovar_ensemble, main= "ensemble") #band 1 is 'mean'
plot(ROCBin_Auto1st_Bovar_ensemble, main= "ROC Bin")


# #### Convert out.rasters to GeoTIFF to be able to visualize in ArcGIS ####
#band 1 is 'mean'
writeRaster(Auto1st_Bovar_ensemble, filename="Auto1st_Bovar_ensemble_1", band = 1, format="GTiff")
writeRaster(ROCBin_Auto1st_Bovar_ensemble, filename="ROCBin_Auto1st_Bovar_ensemble_1", band = 1, format="GTiff")

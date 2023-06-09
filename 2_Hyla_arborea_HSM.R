###############################################################################################################
#### Script for Ensemble Habitat Suitability Modelling ########################################################
### With previous selection of variables and correlation checking #############################################
#### For Hyla_arborea #########################################################################################
#### Author: Damian O. Ortiz-Rodríguez, Antoine Guisan, Rolf Holderegger, Maarten J. van Strien ###############
#### Article: "Sensitivity of habitat network models to changes in maximum dispersal distance" ################
#### Original version developed for article "Predicting species occurrences with habitat network models" ######
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

##### stack of all predictors
PredictorStack <- raster::stack(rasts)

#Write R object to actual .grd raster stack
writeRaster(PredictorStack, filename=("PredictorStack2"), bandorder='BIL', suffix='numbers', overwrite = T)

#############################################
#When done once, for all the other species the previous step can be skipped, just load the saved stack

### load the environmental predictor variables raster stack ###

PredictorStack <- stack('PredictorStack2.grd')
names(PredictorStack)


setwd("E:/PhD/Automated_Run_BioCHECNET/NoPseudorep/HSM")

#### load species data ####
#needs 0's and 1's next to each coordinate pair
Hyla_arborea <- read.csv("./Presences_tables/Hyla_arborea_npr.csv")
head(Hyla_arborea)

# resp.name, the name of focal species (Header of the binary presence column) 
Hy_arb <- 'Hyla_arborea_sa'
# The presence (or presence/absence) data, the resp.var, in my case, that's the column 'Hyla.arborea_sa'
Hy_arbPres <- as.numeric(Hyla_arborea[,Hy_arb])
# the XY coordinates of species data
Hy_arbXY <- Hyla_arborea[,c("CX_HA_C","CY_HA_C")]


######  Check for spatial autocorrelation
#Make an 'Ecospat-Biomod' compatible dataframe
Hyla_arborea_redux <- data.frame(Hyla_arborea$CX_HA_C, Hyla_arborea$CY_HA_C, Hyla_arborea$Hyla_arborea_sa)

Hyla_arborea_redux <- setnames(Hyla_arborea_redux, "Hyla_arborea.CX_HA_C", "CX_HA_C")
Hyla_arborea_redux <- setnames(Hyla_arborea_redux, "Hyla_arborea.CY_HA_C", "CY_HA_C")
Hyla_arborea_redux <- setnames(Hyla_arborea_redux, "Hyla_arborea.Hyla_arborea_sa", "Hyla_arborea_sa")

#extract
Hyla_arborea_xy <- data.frame(Hyla_arborea_redux$CX_HA_C, Hyla_arborea_redux$CY_HA_C)

Hyla_arborea_xy <- setnames(Hyla_arborea_xy, "Hyla_arborea_redux.CX_HA_C", "CX_HA_C")
Hyla_arborea_xy <- setnames(Hyla_arborea_xy, "Hyla_arborea_redux.CY_HA_C", "CY_HA_C")

Hyla_arborea_points <- extract(PredictorStack, Hyla_arborea_xy, method='simple', buffer=NULL, small=FALSE, cellnumbers=FALSE, 
                               fun=NULL, na.rm=TRUE, 1, 20, df=TRUE, factors=FALSE)
#Add the X and Y to Gender_species_points
Hyla_arborea_points$CX_HA_C <- Hyla_arborea_xy$CX_HA_C
Hyla_arborea_points$CY_HA_C <- Hyla_arborea_xy$CY_HA_C

#Mantel correlogram
Hyla_arborea_correlogram <- ecospat.mantel.correlogram (dfvar=Hyla_arborea_points, colxy=22:23, n=100, colvar=2:21, max=1000, nclass=10, nperm=100)
#max, nclass, nperm: from Ecospat vignette


#############################
#### biomod2 modelling ######

#### Set data settings ####

Hy_arb_Data <- BIOMOD_FormatingData(resp.var = Hy_arbPres,
                                    expl.var = PredictorStack,
                                    resp.xy = Hy_arbXY,
                                    resp.name = Hy_arb,
                                    eval.resp.var=NULL, 
                                    eval.expl.var = NULL,
                                    eval.resp.xy = NULL,
                                    PA.nb.rep = 1,
                                    PA.nb.absences = 10000,
                                    PA.strategy = 'random',
                                    na.rm = TRUE)

#check if data are correctly formatted
Hy_arb_Data
plot(Hy_arb_Data)



#### Building models ####  
#keep Biomod2 default options
Hy_arb_ModOptions <-BIOMOD_ModelingOptions()

#In case Maxent isn't working
#library("rJava")
#MAXENT.Phillips = list( path_to_maxent.jar = getwd())
#system.file("java", package="dismo")

### Modelling ###
# Computing the models
Hy_arb_ModelOut <- BIOMOD_Modeling(
  Hy_arb_Data,
  models = c('GLM','RF','MAXENT.Phillips'),
  models.options = Hy_arb_ModOptions,
  NbRunEval=3,
  DataSplit=80,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(Hy_arb,"Auto3",sep=""))


#When this step is over, have a look at some outputs :
# modeling summary
Hy_arb_ModelOut


#### Models evaluations ####

# get all models evaluation
Hy_arb_ModelEval <- get_evaluations(Hy_arb_ModelOut)
# print the dimnames of this object
dimnames(Hy_arb_ModelEval)

#print the ROC scores of all selected models
Hy_arb_ModelEval["ROC","Testing.data",,,]


#print Relative importance of the explanatory variables
get_variables_importance(Hy_arb_ModelOut)

################################################

#### Ensemble modeling ####

Hy_arb_EnsembleMod3 <- BIOMOD_EnsembleModeling(
  modeling.output = Hy_arb_ModelOut,
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
Hy_arb_EnsembleMod3

# get evaluation scores
get_evaluations(Hy_arb_EnsembleMod3)


#### 4 Projection ####

# projection over study area under current conditions
Hy_arb_Auto3_Proj <- BIOMOD_Projection(
  modeling.output = Hy_arb_ModelOut,
  new.env = PredictorStack ,
  proj.name = 'Hy_arb_Auto3rd',
  selected.models = 'all',
  binary.meth = 'ROC',
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

# summary of created object
Hy_arb_Auto3_Proj

# files created on hard drive
list.files("Hyla.arborea.sa/proj_Hy_arb_Auto3rd/")


# make some plots sub-selected by str.grep argument
plot(Hy_arb_Auto3_Proj, str.grep = 'MAXENT.Phillips')
plot(Hy_arb_Auto3_Proj, str.grep = 'GLM')
plot(Hy_arb_Auto3_Proj, str.grep = 'RF')


#get the projected map

Hy_arb_Auto3_Projectedmap <- get_predictions(Hy_arb_Auto3_Proj)
Hy_arb_Auto3_Projectedmap
plot(Hy_arb_Auto3_Projectedmap)

#Visualize maps
getwd()
Auto3rd_Hyarb_HSM <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.grd")
#It's a raster with multiple bands, one for each projection made by each model run
plot(Auto3rd_Hyarb_HSM)

#To check the names and order of the layers in the stack
names(Auto3rd_Hyarb_HSM)

#Get individual layers to check if they are all proper maps
#load the individual bands to a raster object
#use .gri, not .grd

Auto3rd_Hyarb_HSM_Band1 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 1)
Auto3rd_Hyarb_HSM_Band2 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 2)
Auto3rd_Hyarb_HSM_Band3 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 3)
Auto3rd_Hyarb_HSM_Band4 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 4)
Auto3rd_Hyarb_HSM_Band5 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 5)
Auto3rd_Hyarb_HSM_Band6 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 6)
Auto3rd_Hyarb_HSM_Band7 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 7)
Auto3rd_Hyarb_HSM_Band8 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 8)
Auto3rd_Hyarb_HSM_Band9 = raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa.gri", band = 9)

plot(Auto3rd_Hyarb_HSM_Band1, main= "PA1_RUN1_GLM")
plot(Auto3rd_Hyarb_HSM_Band2, main= "PA1_RUN1_RF")
plot(Auto3rd_Hyarb_HSM_Band3, main= "PA1_RUN1_MAXENT.Phillips")
plot(Auto3rd_Hyarb_HSM_Band4, main= "PA1_RUN2_GLM")
plot(Auto3rd_Hyarb_HSM_Band5, main= "PA1_RUN2_RF")
plot(Auto3rd_Hyarb_HSM_Band6, main= "PA1_RUN2_MAXENT.Phillips")
plot(Auto3rd_Hyarb_HSM_Band7, main= "PA1_RUN3_GLM")
plot(Auto3rd_Hyarb_HSM_Band8, main= "PA1_RUN3_RF")
plot(Auto3rd_Hyarb_HSM_Band9, main= "PA1_RUN3_MAXENT.Phillips")

#Same as above but binarized
ROCBin_Auto3rd_Hyarb_HSM <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa_ROCbin.grd")
plot(ROCBin_Auto3rd_Hyarb_HSM)


# #### Convert out.rasters to GeoTIFF to be able to visualize in ArcGIS ####
writeRaster(Auto3rd_Hyarb_HSM, filename="Auto3rd_Hyarb_HSM_stack", format="GTiff")
writeRaster(ROCBin_Auto3rd_Hyarb_HSM, filename="ROCBin_Auto3rd_Hyarb_HSM_stack", format="GTiff")



#### 5 Ensemble_Forecasting ####
Hy_arb_Auto3_EF <- BIOMOD_EnsembleForecasting(
  EM.output = Hy_arb_EnsembleMod3,
  projection.output = Hy_arb_Auto3_Proj,
  binary.meth="ROC") 

Hy_arb_Auto3_EF
plot(Hy_arb_Auto3_EF)

#Visualize ensemble maps
#(binary transformation of single models based on ROC optimized threshold)

#mean by ROC = mean ensemble model
Auto3rd_Hyarb_mean_EM <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/individual_projections/Hyla.arborea.sa_EMmeanByROC_mergedAlgo_mergedRun_mergedData.grd")
#ca by ROC = committee averaging ensemble model 
Auto3rd_Hyarb_ca_EM <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/individual_projections/Hyla.arborea.sa_EMcaByROC_mergedAlgo_mergedRun_mergedData.grd")
#cv by ROC = coefficient of variation - lower score, better model
Auto3rd_Hyarb_cv_EM <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/individual_projections/Hyla.arborea.sa_EMcvByROC_mergedAlgo_mergedRun_mergedData.grd")
#median by ROC = median ensemble model
Auto3rd_Hyarb_median_EM <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/individual_projections/Hyla.arborea.sa_EMmedianByROC_mergedAlgo_mergedRun_mergedData.grd")
#wmean by ROC = weighted mean ensemble models (weighted by single models ROC score)
Auto3rd_Hyarb_wmean_EM <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/individual_projections/Hyla.arborea.sa_EMwmeanByROC_mergedAlgo_mergedRun_mergedData.grd")

#This one is a stack of all the previous 5 plus the ciInf and ciSup   
  #ciInf = Confidence interval Inferior, ciSup = Confidence interval Superior
Auto3rd_Hyarb_ensemble <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa_ensemble.grd")
#ROC Binarization
ROCBin_Auto3rd_Hyarb_ensemble <-raster("./Hyla.arborea.sa/proj_Hy_arb_Auto3rd/proj_Hy_arb_Auto3rd_Hyla.arborea.sa_ensemble_ROCbin.grd")

#Plot to check visually bands
plot(Auto3rd_Hyarb_mean_EM, main= "mean ensemble model")
plot(Auto3rd_Hyarb_ca_EM, main= "comittee_averaging")
plot(Auto3rd_Hyarb_median_EM, main= "median  ensemble model")
plot(Auto3rd_Hyarb_wmean_EM, main= "weighted mean ensemble model")
plot(Auto3rd_Hyarb_cv_EM, main= "coefficient of variation")

plot(Auto3rd_Hyarb_ensemble, main= "ensemble") #band 1 is 'mean'
plot(ROCBin_Auto3rd_Hyarb_ensemble, main= "ROC Bin")

# #### Convert out.rasters to GeoTIFF to be able to visualize in ArcGIS ####
#band 1 is 'mean'
writeRaster(Auto3rd_Hyarb_ensemble, filename="Auto3rd_Hyarb_ensemble_1", band = 1, format="GTiff")
writeRaster(ROCBin_Auto3rd_Hyarb_ensemble, filename="ROCBin_Auto3rd_Hyarb_ensemble_1", band = 1, format="GTiff") 



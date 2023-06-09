#### Script for: ##################################################################################################
#### Definition of Absences #### #### Calculation of predictors #### Boosted Regression Trees (BRT's) Modelling ###
#### For Hyla arborea #### ########################################################################################
###################################################################################################################
#### Author: Damian O. Ortiz-Rodríguez, Antoine Guisan, Maarten J. van Strien #####################################
#### Article: "Sensitivity of habitat network models to changes in maximum dispersal distance" ####################
#### Originally developed for "Predicting species occurrences with habitat network models"#########################
###################################################################################################################
#For application on other species/cases, run until decider plot and from the result replace '6Vt' with the Vt value that corresponds to Vh>=1


library(ggplot2)
library(sp)
library(rgdal)
library(tools)
library(igraph)
library(plyr)
library(dplyr)
library(data.table)
library(nnet)
library(forcats)
library(dismo)
library(ROCR)
library(cutpointr)


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


#### Import 'Macrotable' with all the patches and their presence or 'absence' status ####
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Presence_absence')

Macrotable_Hy_arb <- read.csv("Macrotable_Hyarb.csv")
head(Macrotable_Hy_arb)
names(Macrotable_Hy_arb)
#Names_Macrotable <- names(Macrotable_auto)

#Drop OID  
Macrotable_Hy_arb$OID <- NULL

#Reorder columns, to be consistent with species alphabetical order (after Hyla arborea)
Macrotable_Hy_arb <- Macrotable_Hy_arb[,order(colnames(Macrotable_Hy_arb))] 
head(Macrotable_Hy_arb)


#########Get binary values from the counts ####
#ifelse(df>0, 1, 0 ) 
#converts it to a matrix
Pres_abs_allpatches_Hy_arb <- ifelse(Macrotable_Hy_arb >0, 1, 0 )
#recover data.frame type
Pres_abs_allpatches_Hy_arb <- as.data.frame(Pres_abs_allpatches_Hy_arb, row.names = NULL, optional = FALSE)
names(Pres_abs_allpatches_Hy_arb)
head(Pres_abs_allpatches_Hy_arb)
#Drop binarized value & count
Pres_abs_allpatches_Hy_arb$VALUE <- NULL
Pres_abs_allpatches_Hy_arb$COUNT <- NULL

#Recovering the identifier columns
Pres_abs_allpatches_Hy_arb$Value <- Macrotable_Hy_arb$VALUE
Pres_abs_allpatches_Hy_arb$Count <- Macrotable_Hy_arb$COUNT
names(Pres_abs_allpatches_Hy_arb)

#Get fields VALUE & COUNT (patch area in hectares) to the beginning
Pres_abs_allpatches_Hy_arb <- moveMe(Pres_abs_allpatches_Hy_arb, c("Value", "Count"), "first")
#head(Pres_abs_allpatches_Hy_arb)
names(Pres_abs_allpatches_Hy_arb)


#### Get V_sp - Observations of Gender_species ####
#(sum(1 or 0 each year) #Rows from 43 to 52 (43:52) for Hy_arb 
Pres_abs_allpatches_Hy_arb$V_Hy_arb <- rowSums(Pres_abs_allpatches_Hy_arb[grep("Hy_arb", names(Pres_abs_allpatches_Hy_arb))])
names(Pres_abs_allpatches_Hy_arb)


#### Calculate observations of any amphibian per year ####

for(year in c("06","07","08","09","10","11","12","13","14","15")){
  Pres_abs_allpatches_Hy_arb[,paste("V_", year, sep = "")]<- rowSums(Pres_abs_allpatches_Hy_arb[grep(year, names(Pres_abs_allpatches_Hy_arb))])
}

Pres_abs_allpatches_Hy_arb$V_06 <- ifelse(Pres_abs_allpatches_Hy_arb$V_06>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_07 <- ifelse(Pres_abs_allpatches_Hy_arb$V_07>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_08 <- ifelse(Pres_abs_allpatches_Hy_arb$V_08>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_09 <- ifelse(Pres_abs_allpatches_Hy_arb$V_09>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_10 <- ifelse(Pres_abs_allpatches_Hy_arb$V_10>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_11 <- ifelse(Pres_abs_allpatches_Hy_arb$V_11>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_12 <- ifelse(Pres_abs_allpatches_Hy_arb$V_12>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_13 <- ifelse(Pres_abs_allpatches_Hy_arb$V_13>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_14 <- ifelse(Pres_abs_allpatches_Hy_arb$V_14>0, 1, 0)
Pres_abs_allpatches_Hy_arb$V_15 <- ifelse(Pres_abs_allpatches_Hy_arb$V_15>0, 1, 0)


### Get V_t - Total observations of amphibians
#By summing the vectors of spp per year
Pres_abs_allpatches_Hy_arb$V_t <- rowSums(Pres_abs_allpatches_Hy_arb[,c("V_06","V_07","V_08","V_09","V_10","V_11","V_12","V_13","V_14","V_15")])


################

#Plot Vt vs Vs (Vh)
boxplot(Pres_abs_allpatches_Hy_arb$V_Hy_arb ~ Pres_abs_allpatches_Hy_arb$V_t, xlab= "Vt", ylab = "Vs")
#plot bars indicating median 

write.csv(Pres_abs_allpatches_Hy_arb, file = "Pres_abs_Vs_Vt_Hy_arb.csv") 

#Total number of patches
length(Pres_abs_allpatches_Hy_arb$V_t)


################ Occurence-state (Pres_abs) decider Plots #########################

#Pres_abs_allpatches2 is the one modified to make plots, excludes all the NA 
#Pres_abs_allpatches (unchanged, with all the 0's) will be used to incorporate the graph properties
#Get only the records that are visits (Vt >= 1)
Pres_abs_allpatches_Hy_arb2 = Pres_abs_allpatches_Hy_arb[Pres_abs_allpatches_Hy_arb$V_t != 0,] 

#Bar graph with proportions, averages and text --> The main one to decide thresholds!!!

#Data for the plot #Taking mean Vh
tab = ddply(Pres_abs_allpatches_Hy_arb2, .(V_t), summarize,  mean=mean(V_Hy_arb), countH = length(which(V_Hy_arb>0)), countT = length(V_Hy_arb), lab = paste(length(which(V_Hy_arb>0)),"/",length(V_Hy_arb),sep = " "), perc = length(which(V_Hy_arb>0))/length(V_Hy_arb)*100)

## Plot to decide threshold ##
my_y_title <- "Mean Vs" 

ggplot(tab, aes(factor(V_t)))+
  geom_col(aes(y=mean), fill = "transparent", colour="black", size=2)+
  scale_y_continuous(breaks = c(1:10))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black",size=2), axis.text.x= element_text(size=12, face="bold"), axis.text.y = element_text(size=12, face="bold"), axis.ticks = element_line(size = 2))+
  labs(y = my_y_title, x = "Vt")


#### Count and identify absences 

#For value of mean_Vt corresponding to a Vs which looks over 1 in the plot, 
# check exact value to make surre it's 1 or over.
# If yes, then it's the Vt *threshold value* for absences
# E.g.#Vt=6 -> Vh=1.8125 #Threshold value
# (See Ortiz-Rodriguez et al 2019)
mean(Pres_abs_allpatches_Hy_arb$V_Hy_arb[Pres_abs_allpatches_Hy_arb$V_t==6])

#Number of absences at threshold value:
sum(Pres_abs_allpatches_Hy_arb$V_t>5 & Pres_abs_allpatches_Hy_arb$V_Hy_arb<1) 

which(Pres_abs_allpatches_Hy_arb$V_t >= 6 & Pres_abs_allpatches_Hy_arb$V_Hy_arb ==0)


#### rename(data, d = b) ####
Pres_abs_allpatches_Hy_arb <- setnames(Pres_abs_allpatches_Hy_arb, "Value", "PatchID")
Pres_abs_allpatches_Hy_arb <- setnames(Pres_abs_allpatches_Hy_arb, "Count", "Patch_Area")


#### Add HSI ####
HSI = read.csv("C:/Users/damiano/Documents/PhD/Automated_Run_BioCHECNET/Corrected/HSI_habPatchCode_NoPseudorep.csv")
Pres_abs_allpatches_Hy_arb$HSI <- HSI$"MEAN"



######## Importing network ############

setwd('C:/Users/damiano/Documents/PhD/Automated_Run_BioCHECNET/Corrected')

#Import table of patches w/network properties
Hyarb_graph <- read_graph("uniform.graphml", format = "graphml") 

#Look at vertices and edges of a graph
E(Hyarb_graph)
V(Hyarb_graph)


#### get metrics as attributes and into a data.frame ####
Hyarb_graph_metrics <- data.frame(
  PatchID = V(Hyarb_graph)$name,
  deg=degree(Hyarb_graph),
  strength = strength(Hyarb_graph, vids = V(Hyarb_graph), mode = "all", loops = FALSE, weights = E(Hyarb_graph)$weight),
  EgoSize = ego_size(Hyarb_graph, order = 3, nodes = V(Hyarb_graph))
  
)
### EgoSize = Third-order neighborhood


#### join the tables by attribute
#inner join= join by common variable names
patches_Hyarb_topo_attributes <- merge(Pres_abs_allpatches_Hy_arb, Hyarb_graph_metrics, by = "PatchID")


#### Add Habitat Availability predictor ####
Hyarb_habAv <- scan("uniform_habAv.txt")
patches_Hyarb_topo_attributes$habAv <- Hyarb_habAv
head(patches_Hyarb_topo_attributes$habAv)


#### Measure betweenness centrality without weight included ####
#Re-import graphs, change name of 'weight' attribute
#setwd('/Networks')
Hyarb_graph2 <- read_graph("uniform.graphml", format = "graphml")

#rename weight attribute
#Actually copying into a new attr. and deleting original attribute
E(Hyarb_graph2)$Alter_wght <- E(Hyarb_graph2)$weight
Hyarb_graph2 <- remove.edge.attribute(Hyarb_graph2, "weight")

#Calculate explicitely unweighted b_c
#Set dataframe w/only PatchID & unw_b_c
Hyarb_unweighted_b_c <- data.frame(
  PatchID = V(Hyarb_graph2)$name,
  b_c=betweenness(Hyarb_graph2, weights = NULL))

#rename unw_b_c
Hyarb_unweighted_b_c <- rename(Hyarb_unweighted_b_c,  "unw_b_c" = "b_c")

#join the tables by attribute
patches_Hyarb_topo_attributes <- merge(x = patches_Hyarb_topo_attributes, y = Hyarb_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(patches_Hyarb_topo_attributes)
head(patches_Hyarb_topo_attributes)



#### Define occurrence_state (pres_abs) ####
### Classify the dataframe in subsets that classify all the patches in presences, absences or questionmarks/NoIdea (NA) 

#subset of questionmarks, #All the records, except the ones that comply with the following commands
patches_Hyarb_topo_attributes[,"pres_abs"] = NA 
#Subset of likely absences (defined by pres_abs decider plot)
patches_Hyarb_topo_attributes[which(patches_Hyarb_topo_attributes$V_t >= 6 & patches_Hyarb_topo_attributes$V_Hy_arb ==0),"pres_abs"] = 0
#Subset of confirmed presences
patches_Hyarb_topo_attributes[which(patches_Hyarb_topo_attributes$V_Hy_arb > 0),"pres_abs"] = 1

#### Display the 'occupancy state' (0/1/NA) of every patch as defined above (pres_abs) ####
patches_Hyarb_topo_attributes$pres_abs

#Get total amount of values for presence and for absence
table(patches_Hyarb_topo_attributes$pres_abs)


#### Do a subset of the columns that only has the predictors and response ####
Hyarb_stattest <- subset(patches_Hyarb_topo_attributes, select=c("PatchID", "Patch_Area", "V_Hy_arb", "V_t", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Hyarb_stattest)
head(Hyarb_stattest)

#change directory to the one with BRT things 
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/BRTs')
#write file
write.csv(Hyarb_stattest, file = "Hyarb_stattest_t6.csv")
#t6 -> pres_abs threshold=6


##### Do the 'only pres_abs (no-NoData) df's #########################
#Take out NA's, set threshold of absences, change name
Hyarb_6Vt_woNA = Hyarb_stattest
Hyarb_6Vt_woNA[Hyarb_6Vt_woNA$V_Hy_arb==0 & Hyarb_6Vt_woNA$V_t >= 6, 'pres_abs'] = 0
Hyarb_6Vt_woNA = Hyarb_6Vt_woNA[!is.na(Hyarb_6Vt_woNA$pres_abs),]

#check values of unw_b_c
mean(Hyarb_6Vt_woNA$unw_b_c)



#### Perform BRT's ############################################################
#Hyperparameters for all models
lr = 0.001
tc = 5
bf = 0.75

n_repeats = 100

##### Hyarb (Full model with all the predictors) ###########################################

# Create the output table that contains all the values of each of the runs
Hyarb_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb
# Create vectors to record performance measures
Hyarb_AUC_cv_vec = vector() #Cross-validated AUC
Hyarb_AUC_vec = vector() #Training AUC

Hyarb_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Hyarb_HSI_vec = vector() 
Hyarb_EgoSize_vec = vector() 
Hyarb_strength_vec  = vector()
Hyarb_deg_vec  = vector()
Hyarb_habAv_vec  = vector()
Hyarb_unw_b_c_vec  = vector()
Hyarb_Patch_Area_vec  = vector()

#Prediction of occurrence state for all patches in the present
Hyarb_predict_mat = matrix(nrow = length(Hyarb_stattest$PatchID), ncol = n_repeats, byrow = FALSE)


### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Hyarb = gbm.step(data=Hyarb_6Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Hyarb)$var, imp = summary(gbm_mod_Hyarb)$rel.inf)
  
  #Write var_imp results to the vectors
  Hyarb_HSI_vec = append(Hyarb_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Hyarb_EgoSize_vec = append(Hyarb_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Hyarb_strength_vec  = append(Hyarb_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Hyarb_deg_vec  = append(Hyarb_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Hyarb_habAv_vec  = append(Hyarb_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Hyarb_unw_b_c_vec  = append(Hyarb_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Hyarb_Patch_Area_vec  = append(Hyarb_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Hyarb_AUC_cv_vec = append(Hyarb_AUC_cv_vec, gbm_mod_Hyarb$cv.statistics$discrimination.mean)
  Hyarb_AUC_vec = append(Hyarb_AUC_vec, gbm_mod_Hyarb$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Hyarb$n.trees
  Hyarb_nt_vec = append(Hyarb_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Hyarb_predict_mat[,i] = predict(gbm_mod_Hyarb, Hyarb_stattest, gbm_mod_Hyarb$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}


Hyarb_output_tab[,"AUC_train"] = Hyarb_AUC_vec
Hyarb_output_tab[,"AUC_cv"] = Hyarb_AUC_cv_vec
Hyarb_output_tab[,"ntrees"] = Hyarb_nt_vec

#Var. importance columns
Hyarb_output_tab[,"HSI_imp"] = Hyarb_HSI_vec
Hyarb_output_tab[,"EgoSize_imp"] = Hyarb_EgoSize_vec
Hyarb_output_tab[,"strength_imp"] = Hyarb_strength_vec
Hyarb_output_tab[,"deg_imp"] = Hyarb_deg_vec
Hyarb_output_tab[,"habAv_imp"] = Hyarb_habAv_vec
Hyarb_output_tab[,"unw_b_c_imp"] = Hyarb_unw_b_c_vec
Hyarb_output_tab[,"Patch_Area_imp"] = Hyarb_Patch_Area_vec


### Make a dataframe for plotting overlaying histrograms in R
Hyarb_HSI_tab = data.frame(imp = Hyarb_HSI_vec)
Hyarb_EgoSize_tab = data.frame(imp = Hyarb_EgoSize_vec)
Hyarb_strength_tab = data.frame(imp = Hyarb_strength_vec)
Hyarb_deg_tab = data.frame(imp = Hyarb_deg_vec)
Hyarb_habAv_tab = data.frame(imp = Hyarb_habAv_vec)
Hyarb_unw_b_c_tab = data.frame(imp = Hyarb_unw_b_c_vec)
Hyarb_Patch_Area_tab = data.frame(imp = Hyarb_Patch_Area_vec)

Hyarb_HSI_tab$variable = "HSI"
Hyarb_EgoSize_tab$variable = "3rd. ord. neigh."
Hyarb_strength_tab$variable = "Strength"
Hyarb_deg_tab$variable = "Degree"
Hyarb_habAv_tab$variable = "Hab. Av."
Hyarb_unw_b_c_tab$variable = "B.C."
Hyarb_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Hyarb_AUC_cv_tab = data.frame(value = Hyarb_AUC_cv_vec)
Hyarb_AUC_train_tab = data.frame(value = Hyarb_AUC_vec)
#Make label of model for plot
Hyarb_AUC_cv_tab$model = "Hyarb"
Hyarb_AUC_train_tab$model = "Hyarb"

#combine pred. vars. into new data frame 
Hyarb_var_imp_tab = rbind(Hyarb_HSI_tab,Hyarb_EgoSize_tab,Hyarb_strength_tab,Hyarb_habAv_tab,Hyarb_deg_tab,Hyarb_unw_b_c_tab,Hyarb_Patch_Area_tab)

ggplot(Hyarb_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Hyarb_var_imp_tab$imp~Hyarb_var_imp_tab$variable, ylab= "Var. Importance", main = "Hyarb")

#Get mean var. importance of all of the vars. 
mean(Hyarb_HSI_vec)
mean(Hyarb_EgoSize_vec)
mean(Hyarb_strength_vec)
mean(Hyarb_deg_vec)
mean(Hyarb_habAv_vec)
mean(Hyarb_unw_b_c_vec)
mean(Hyarb_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Hyarb_output_tab$AUC_cv)
summary(Hyarb_output_tab$AUC_train)



############ NoTopo (Model with only non-topological predictors included) #################################

# Create the output table that contains all the values of each of the runs
noTopo_output_tab = data.frame(run_nr = c(1:n_repeats)) #Hyarb
# Create vectors to record performance measures
noTopo_AUC_cv_vec = vector() #Cross-validated AUC
noTopo_AUC_vec = vector() # Training AUC
noTopo_nt_vec = vector() # number of trees
# noTopo_predict_vec = vector()

#Create vectors with all the values of a certain predictor along all the runs
noTopo_HSI_vec = vector() 
noTopo_Patch_Area_vec  = vector()

### Loop that goes exactly for 100 iterations, to get distributions 

for(i in c(1:n_repeats)){

  
  #Perform gbm with a fixed number of trees and no cross-validation
  gbm_mod_noTopo = gbm.step(data=Hyarb_6Vt_woNA, gbm.x = c('Patch_Area','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_noTopo)$var, imp = summary(gbm_mod_noTopo)$rel.inf)
  
  #Write var_imp results to the vectors
  noTopo_HSI_vec = append(noTopo_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  noTopo_Patch_Area_vec  = append(noTopo_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the CV-AUC of this model to a vector.
  noTopo_AUC_cv_vec = append(noTopo_AUC_cv_vec, gbm_mod_noTopo$cv.statistics$discrimination.mean)
  noTopo_AUC_vec = append(noTopo_AUC_vec, gbm_mod_noTopo$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_noTopo$n.trees
  noTopo_nt_vec = append(noTopo_nt_vec, nt)
  print(nt)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
  
}

noTopo_output_tab[,"AUC_train"] = noTopo_AUC_vec
noTopo_output_tab[,"AUC_cv"] = noTopo_AUC_cv_vec
noTopo_output_tab[,"ntrees"] = noTopo_nt_vec
#Var. importance columns
noTopo_output_tab[,"HSI_imp"] = noTopo_HSI_vec
noTopo_output_tab[,"Patch_Area_imp"] = noTopo_Patch_Area_vec


### Make a dataframe to make plots for comparison
noTopo_HSI_tab = data.frame(imp = noTopo_HSI_vec)
noTopo_Patch_Area_tab = data.frame(imp = noTopo_Patch_Area_vec)

noTopo_HSI_tab$variable = "HSI"
noTopo_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to compare performance between models
noTopo_AUC_cv_tab = data.frame(value = noTopo_AUC_cv_vec)
noTopo_AUC_train_tab = data.frame(value = noTopo_AUC_vec)

#Make label of model for plot
noTopo_AUC_cv_tab$model = "noTopo"
noTopo_AUC_train_tab$model = "noTopo"

#combine pred. vars. into new data frame 
noTopo_var_imp_tab = rbind(noTopo_HSI_tab,noTopo_Patch_Area_tab)

ggplot(noTopo_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(noTopo_var_imp_tab$imp~noTopo_var_imp_tab$variable, ylab= "Var. Importance", main = "noTopo_Hyarb")

#Get mean var. importance of all of the vars. 
mean(noTopo_HSI_vec)
mean(noTopo_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(noTopo_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_cv)



##### Compare between runs ###########################

#cv AUC
#mean
AUC_cv_tab = rbind(Hyarb_AUC_cv_tab, noTopo_AUC_cv_tab)

#Change order of factors to display noTopo at the edge
AUC_cv_tab$model<- as.factor(AUC_cv_tab$model)
levels(AUC_cv_tab$model)
AUC_cv_tab$model<-factor(AUC_cv_tab$model, levels=c("Hyarb", "noTopo"))
#print out
levels(AUC_cv_tab$model)

#Plot
ggplot(AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_cv_tab$value~AUC_cv_tab$model, ylab= "Cross-validated AUC")


#Training AUC
AUC_train_tab = rbind(Hyarb_AUC_train_tab, noTopo_AUC_train_tab)

ggplot(AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_train_tab$value~AUC_train_tab$model, ylab= "Training AUC")


summary(Hyarb_output_tab$AUC_cv)
summary(noTopo_output_tab$AUC_cv)

summary(Hyarb_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_train)

summary(Hyarb_output_tab$ntrees)
summary(noTopo_output_tab$ntrees)


###Get avg_#trees over the AUC>=.75 threshold
Hyarb_mean_acceptable_ntrees <- mean(Hyarb_output_tab$ntrees[Hyarb_output_tab$AUC_cv>=0.75])
Hyarb_mean_acceptable_ntrees
noTopo_mean_acceptable_ntrees <- mean(noTopo_output_tab$ntrees[noTopo_output_tab$AUC_cv>=0.75])
noTopo_mean_acceptable_ntrees
Total_mean_acceptable_ntrees <- mean(Hyarb_mean_acceptable_ntrees,noTopo_mean_acceptable_ntrees)
Total_mean_acceptable_ntrees



#### Get sample models for response curves (tagged with the kind of network they are) ####
#gbm.mod (No tag) = noTopo model (as it was the last one to be run)

#NoTopo
test = gbm.plot(gbm_mod_noTopo,return.grid=TRUE)
str(test)
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)
#Hyarb
gbm_mod_Hyarb = gbm.step(data=Hyarb_6Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','habAv','EgoSize','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

#plot
gbm.plot(gbm_mod_Hyarb, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,3), write.title=FALSE, cex.axis = 2,  cex.lab=2)
#NoTopo
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)


### specific stats of sample models:
summary(gbm_mod_Hyarb)
gbm_mod_Hyarb$self.statistics$discrimination
#noTopo
summary(gbm_mod_noTopo)
gbm_mod_noTopo$self.statistics$discrimination

write.csv(Hyarb_6Vt_woNA, file = "BinPredRun_Hyarb_6Vt_woNA_b_c_unw.csv")



#################################################################################################
#### Get discrete and rescaled predictions ####
#################################################################################################
#### Add continuous predictions to dataframe w/ all the predictors & occurrence-state ####
#Do new dataframe to keep the original unchanged (only add underscore between gender & sp.)
# Hy_arb_stattest <- Hyarb_stattest
# 
# 
# Hyarb_predict_df = data.frame(Hyarb_predict_mat)
# #Calculate the mean of the 100 models for each row (patch)
# # For past round it was done as uniform_stattest$cont_pred <- Singlemodel_Hy_arb_predict
# Hyarb_predict_df$mean_pred <- rowMeans(Hyarb_predict_mat, na.rm = FALSE)
# 
# Hyarb_predict_df$PatchID <- Hy_arb_stattest$PatchID
# 
# Hyarb_predict_df <- merge(Hy_arb_stattest, Hyarb_predict_df, by = "PatchID")
# names(Hyarb_predict_df)
# 
# ## Set discrete prediction to patches
# #Set discretization (binariation) cutpoint for each patch
# cp = cutpointr(Hyarb_predict_df, mean_pred, pres_abs, subgroup = NULL, method = minimize_metric, metric = roc01, na.rm = TRUE)
# 
# summary(cp)
# plot(cp)
# cp$optimal_cutpoint
# #Set all above threshold to 1, all below to 0 (discrete = disc)
# Hyarb_predict_df$pred_pres_abs <- ifelse(Hyarb_predict_df$mean_pred>=cp$optimal_cutpoint, 1, 0)
# 
# table(Hyarb_predict_df$pred_pres_abs)
# 
# #Add field w/ BRT Prediction score multiplied by 1000 (Mean_occurrence_index)
# #rescaled continuous prediction, to correspond to scale of HSI
# Hyarb_predict_df$Pred_1000 <- Hyarb_predict_df$mean_pred*1000
# 
# 
# #Export df to csv without the columns used only for calculation of binary prediction
# Hy_arb_stattest$mean_pred <- Hyarb_predict_df$mean_pred
# Hy_arb_stattest$Pred_1000 <-Hyarb_predict_df$Pred_1000
# Hy_arb_stattest$bin_OLpred <- Hyarb_predict_df$pred_pres_abs


### For ArcGIS importing (compatibility of fields)
Hy_arb_stattest$Value <- Hy_arb_stattest$PatchID
Hy_arb_stattest$V_s <- Hy_arb_stattest$V_Hy_arb
head(Hy_arb_stattest)

Hy_arb_stattest <- moveMe(Hy_arb_stattest, c("Value"), "after", "PatchID")
Hy_arb_stattest <- moveMe(Hy_arb_stattest, c("V_s"), "after", "V_Hy_arb")


# setwd("C:/Users/damiano/Documents/PhD/Additional_species_runs/Discrete_predictions")
# write.csv(Hy_arb_stattest, file = "Hy_arb_bin_prediction.csv")



#### Hard-save evaluation measures ####
setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
write.csv(Hyarb_AUC_cv_tab, file = "Hyarb_DefaultDispDist_AUC_cv_tab.csv")
write.csv(Hyarb_AUC_train_tab, file = "Hyarb_DefaultDispDist_AUC_train_tab.csv")
write.csv(Hyarb_var_imp_tab, file = "Hyarb_DefaultDispDist_var_imp_tab.csv")
write.csv(Hyarb_output_tab, file = "Hyarb_DefaultDispDist_BRToutput_tab.csv")
# write.csv(Hyarb_predict_df, file = "Hyarb_DefaultDispDist_predict_df.csv")

#Hard-save variable importance of each predictor
write.csv(Hyarb_HSI_tab, file = "Hyarb_DefaultDispDist_HSI_tab.csv")
write.csv(Hyarb_habAv_tab, file = "Hyarb_DefaultDispDist_habAv_tab.csv")
write.csv(Hyarb_deg_tab, file = "Hyarb_DefaultDispDist_deg_tab.csv")
write.csv(Hyarb_EgoSize_tab, file = "Hyarb_DefaultDispDist_EgoSize_tab.csv")
write.csv(Hyarb_unw_b_c_tab, file = "Hyarb_DefaultDispDist_unw_b_c_tab.csv")
write.csv(Hyarb_strength_tab, file = "Hyarb_DefaultDispDist_strength_tab.csv")
write.csv(Hyarb_Patch_Area_tab, file = "Hyarb_DefaultDispDist_Patch_Area_tab.csv")

#Hard-save model to be  able to use it in other workspaces
saveRDS(gbm_mod_Hyarb, file = "gbm_mod_Hyarb_DefaultDispDist.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# readRDS(file, refhook = NULL)
# infoRDS(file)

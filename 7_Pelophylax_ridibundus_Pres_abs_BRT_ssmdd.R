#### Script for: ##################################################################################################
#### Definition of Absences #### #### Calculation of predictors #### Boosted Regression Trees (BRT's) Modelling ###
#### For Pelophylax ridibundus #### #################################################################################
###################################################################################################################
#### Author: Damian O. Ortiz-Rodríguez, Antoine Guisan, Maarten J. van Strien #####################################
#### Article: "Sensitivity of habitat network models to changes in maximum dispersal distance" ####################
###################################################################################################################
#For application on other species/cases, run until decider plot and from the result replace 'X'Vt with the Vt value that corresponds to Vh>=1


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



#Import 'Macrotable'with all the patches and their presence or 'absence' status
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Presence_absence')

Macrotable_Pe_rid <- read.csv("Macrotable_Perid.csv")
head(Macrotable_Pe_rid)
names(Macrotable_Pe_rid)


#Drop OID  
Macrotable_Pe_rid$OID <- NULL

#Reorder columns, to be consistent with species alphabetical order
Macrotable_Pe_rid <- Macrotable_Pe_rid[,order(colnames(Macrotable_Pe_rid))] 
head(Macrotable_Pe_rid)

#########Get binary values from the counts ####
#ifelse(df>0, 1, 0 ) 
#converts it to a matrix
Pres_abs_allpatches_Pe_rid <- ifelse(Macrotable_Pe_rid >0, 1, 0 )
#recover data.frame type
Pres_abs_allpatches_Pe_rid <- as.data.frame(Pres_abs_allpatches_Pe_rid, row.names = NULL, optional = FALSE)
names(Pres_abs_allpatches_Pe_rid)
head(Pres_abs_allpatches_Pe_rid)
#Drop binarized value & count
Pres_abs_allpatches_Pe_rid$Value <- NULL
Pres_abs_allpatches_Pe_rid$Count <- NULL

#Recovering the identifier columns
Pres_abs_allpatches_Pe_rid$Value <- Macrotable_Pe_rid$Value
Pres_abs_allpatches_Pe_rid$Count <- Macrotable_Pe_rid$Count
names(Pres_abs_allpatches_Pe_rid)

#Get fields VALUE & COUNT to the beginning
Pres_abs_allpatches_Pe_rid <- moveMe(Pres_abs_allpatches_Pe_rid, c("Value", "Count"), "first")
#head(Pres_abs_allpatches_Pe_rid)
names(Pres_abs_allpatches_Pe_rid)

#### Get V_sp - Observations of Gender_species ####
#(sum(1 or 0 each year) 
Pres_abs_allpatches_Pe_rid$V_Pe_rid <- rowSums(Pres_abs_allpatches_Pe_rid[grep("Pe_rid", names(Pres_abs_allpatches_Pe_rid))])
names(Pres_abs_allpatches_Pe_rid)


#### Calculate observations of any amphibian per year ####

for(year in c("06","07","08","09","10","11","12","13","14","15")){
  Pres_abs_allpatches_Pe_rid[,paste("V_", year, sep = "")]<- rowSums(Pres_abs_allpatches_Pe_rid[grep(year, names(Pres_abs_allpatches_Pe_rid))])
}

Pres_abs_allpatches_Pe_rid$V_06 <- ifelse(Pres_abs_allpatches_Pe_rid$V_06>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_07 <- ifelse(Pres_abs_allpatches_Pe_rid$V_07>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_08 <- ifelse(Pres_abs_allpatches_Pe_rid$V_08>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_09 <- ifelse(Pres_abs_allpatches_Pe_rid$V_09>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_10 <- ifelse(Pres_abs_allpatches_Pe_rid$V_10>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_11 <- ifelse(Pres_abs_allpatches_Pe_rid$V_11>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_12 <- ifelse(Pres_abs_allpatches_Pe_rid$V_12>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_13 <- ifelse(Pres_abs_allpatches_Pe_rid$V_13>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_14 <- ifelse(Pres_abs_allpatches_Pe_rid$V_14>0, 1, 0)
Pres_abs_allpatches_Pe_rid$V_15 <- ifelse(Pres_abs_allpatches_Pe_rid$V_15>0, 1, 0)


### Get V_t - Total observations of amphibians
#By summing the vectors of spp per year
# Pres_abs_allpatches_Pe_rid$V_t <- rowSums(Pres_abs_allpatches_Pe_rid[,132:141])
Pres_abs_allpatches_Pe_rid$V_t <- rowSums(Pres_abs_allpatches_Pe_rid[,c("V_06","V_07","V_08","V_09","V_10","V_11","V_12","V_13","V_14","V_15")])


################

#Plot Vt vs Vh
boxplot(Pres_abs_allpatches_Pe_rid$V_Pe_rid ~ Pres_abs_allpatches_Pe_rid$V_t, xlab= "Vt", ylab = "Vs")
#plot bars indicating median 

write.csv(Pres_abs_allpatches_Pe_rid, file = "Pres_abs_Vs_Vt_Pe_rid.csv") 

#Total number of patches (2735)
length(Pres_abs_allpatches_Pe_rid$V_t)


################ Occurence-state (Pres_abs) decider Plots #########################

#Pres_abs_allpatches2 is the one modified to make plots, excludes all the NA 
#Pres_abs_allpatches (unchanged, with all the 0's) will be used to incorporate the graph properties
#Get only the records that are visits (Vt >= 1)
Pres_abs_allpatches_Pe_rid2 = Pres_abs_allpatches_Pe_rid[Pres_abs_allpatches_Pe_rid$V_t != 0,] 

#Bar graph with proportions, averages and text --> The main one to decide thresholds!!!

#Data for the plot #Taking mean Vs
tab = ddply(Pres_abs_allpatches_Pe_rid2, .(V_t), summarize,  mean=mean(V_Pe_rid), countH = length(which(V_Pe_rid>0)), countT = length(V_Pe_rid), lab = paste(length(which(V_Pe_rid>0)),"/",length(V_Pe_rid),sep = " "), perc = length(which(V_Pe_rid>0))/length(V_Pe_rid)*100)

## Plot to decide threshold ##
my_y_title <- "Mean Vs"

ggplot(tab, aes(factor(V_t)))+
  geom_col(aes(y=mean), fill = "transparent", colour="black")+
  scale_y_continuous(breaks = c(1:10))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(y = my_y_title, x = "Vt")

#### Count and identify absences 

#For value of mean_Vt corresponding to a Vs which looks over 1 in the plot, 
# check exact value to make surre it's 1 or over.
# If yes, then it's the Vt *threshold value* for absences
# E.g. #Vt=6 -> Vh=1.8125 #Threshold value
# (See Ortiz-Rodriguez et al 2019)
mean(Pres_abs_allpatches_Pe_rid$V_Pe_rid[Pres_abs_allpatches_Pe_rid$V_t==4])

#Number of absences at threshold value:
sum(Pres_abs_allpatches_Pe_rid$V_t>3 & Pres_abs_allpatches_Pe_rid$V_Pe_rid<1) 

which(Pres_abs_allpatches_Pe_rid$V_t >= 4 & Pres_abs_allpatches_Pe_rid$V_Pe_rid ==0)


#rename(data, d = b)
Pres_abs_allpatches_Pe_rid <- setnames(Pres_abs_allpatches_Pe_rid, "Value", "PatchID")
Pres_abs_allpatches_Pe_rid <- setnames(Pres_abs_allpatches_Pe_rid, "Count", "Patch_Area")


#### Add HSI ####
HSI = read.csv("C:/Users/damiano/Documents/PhD/Additional_species_runs/Network_setup/HSI_habPatchCode_Perid1.csv")
Pres_abs_allpatches_Pe_rid$HSI <- HSI$"MEAN"



######## Importing network ############

setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Networks')

#Import table of patches w/network properties
Perid_graph <- read_graph("Perid.graphml", format = "graphml") 

#Look at vertices and edges of a graph
E(Perid_graph)
V(Perid_graph)


#### get metrics as attributes and into a data.frame ####
Perid_graph_metrics <- data.frame(
  PatchID = V(Perid_graph)$name,
  deg=degree(Perid_graph),
  strength = strength(Perid_graph, vids = V(Perid_graph), mode = "all", loops = FALSE, weights = E(Perid_graph)$weight),
  EgoSize = ego_size(Perid_graph, order = 3, nodes = V(Perid_graph))
  
)
### EgoSize = Third-order neighborhood


#### join the tables by attribute
#inner join= join by common variable names
patches_Perid_topo_attributes <- merge(Pres_abs_allpatches_Pe_rid, Perid_graph_metrics, by = "PatchID")


#### Add Habitat Availability predictor ####
Perid_habAv <- scan("Perid_habAv.txt")
patches_Perid_topo_attributes$habAv <- Perid_habAv
head(patches_Perid_topo_attributes$habAv)


#### Measure betweenness centrality without weight included ####
#Re-import graphs, change name of 'weight' attribute
#setwd('/Networks')
Perid_graph2 <- read_graph("Perid.graphml", format = "graphml")

#rename weight attribute
#Actually copying into a new attr. and deleting original attribute
E(Perid_graph2)$Alter_wght <- E(Perid_graph2)$weight
Perid_graph2 <- remove.edge.attribute(Perid_graph2, "weight")

#Calculate explicitely unweighted b_c
#Set dataframe w/only PatchID & unw_b_c
Perid_unweighted_b_c <- data.frame(
  PatchID = V(Perid_graph2)$name,
  b_c=betweenness(Perid_graph2, weights = NULL))

#rename unw_b_c
Perid_unweighted_b_c <- rename(Perid_unweighted_b_c,  "unw_b_c" = "b_c")

#join the tables by attribute
patches_Perid_topo_attributes <- merge(x = patches_Perid_topo_attributes, y = Perid_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(patches_Perid_topo_attributes)
head(patches_Perid_topo_attributes)


#### Define occurrence_state (pres_abs) ####
### Classify the dataframe in subsets that classify all the patches in presences, absences or questionmarks/NoIdea (NA) 

#subset of questionmarks, #All the records, except the ones that comply with the following commands
patches_Perid_topo_attributes[,"pres_abs"] = NA 
#Subset of likely absences (With threshold at 5 (defined by pres_abs decider plot))
patches_Perid_topo_attributes[which(patches_Perid_topo_attributes$V_t >= 4 & patches_Perid_topo_attributes$V_Pe_rid ==0),"pres_abs"] = 0
#Subset of confirmed presences
patches_Perid_topo_attributes[which(patches_Perid_topo_attributes$V_Pe_rid > 0),"pres_abs"] = 1

#### Display the 'occurrence-state' (0/1/NA) of every patch as defined above (pres_abs) ####
patches_Perid_topo_attributes$pres_abs

#Get total amount of values for presence and for absence
table(patches_Perid_topo_attributes$pres_abs)


#### Do a subset of the columns that only has the predictors and response ####
Perid_stattest <- subset(patches_Perid_topo_attributes, select=c("PatchID", "Patch_Area", "V_Pe_rid", "V_t", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Perid_stattest)
head(Perid_stattest)

#change directory to the one with BRT things 
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/BRTs')
#write file
#t4 -> pres_abs threshold=4
write.csv(Perid_stattest, file = "Perid_stattest_t4.csv")


##### Do the 'only pres_abs (no-NoData) df's #########################
#Take out NA's, set threshold of absences (at V_t = 4 for Pe_rid), change name
Perid_4Vt_woNA = Perid_stattest
Perid_4Vt_woNA[Perid_4Vt_woNA$V_Pe_rid==0 & Perid_4Vt_woNA$V_t >= 4, 'pres_abs'] = 0
Perid_4Vt_woNA = Perid_4Vt_woNA[!is.na(Perid_4Vt_woNA$pres_abs),]

#check values of unw_b_c
mean(Perid_4Vt_woNA$unw_b_c)



#### Perform BRT's ############################################################
#### w/only unw_b_c

##Run 100 models without gridsearch, see results, how variable they are
#Hyperparameters for all models
lr = 0.001
tc = 5
bf = 0.75

n_repeats = 100

##### Perid (Full model with all the predictors) ###########################################

# Create the output table that contains all the values of each of the runs
Perid_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid
# Create vectors to record performance measures
Perid_AUC_cv_vec = vector() #Cross-validated AUC
Perid_AUC_vec = vector() #Training AUC
Perid_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Perid_HSI_vec = vector() 
Perid_EgoSize_vec = vector() 
Perid_strength_vec  = vector()
Perid_deg_vec  = vector()
Perid_habAv_vec  = vector()
Perid_unw_b_c_vec  = vector()
Perid_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Perid_predict_mat = matrix(nrow = length(Perid_stattest$PatchID), ncol = n_repeats, byrow = FALSE)


### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Perid = gbm.step(data=Perid_4Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Perid)$var, imp = summary(gbm_mod_Perid)$rel.inf)
  
  #Write var_imp results to the vectors
  Perid_HSI_vec = append(Perid_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Perid_EgoSize_vec = append(Perid_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Perid_strength_vec  = append(Perid_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Perid_deg_vec  = append(Perid_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Perid_habAv_vec  = append(Perid_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Perid_unw_b_c_vec  = append(Perid_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Perid_Patch_Area_vec  = append(Perid_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Perid_AUC_cv_vec = append(Perid_AUC_cv_vec, gbm_mod_Perid$cv.statistics$discrimination.mean)
  Perid_AUC_vec = append(Perid_AUC_vec, gbm_mod_Perid$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Perid$n.trees
  Perid_nt_vec = append(Perid_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Perid_predict_mat[,i] = predict(gbm_mod_Perid, Perid_stattest, gbm_mod_Perid$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}


Perid_output_tab[,"AUC_train"] = Perid_AUC_vec
Perid_output_tab[,"AUC_cv"] = Perid_AUC_cv_vec
Perid_output_tab[,"ntrees"] = Perid_nt_vec


#Var. importance columns
Perid_output_tab[,"HSI_imp"] = Perid_HSI_vec
Perid_output_tab[,"EgoSize_imp"] = Perid_EgoSize_vec
Perid_output_tab[,"strength_imp"] = Perid_strength_vec
Perid_output_tab[,"deg_imp"] = Perid_deg_vec
Perid_output_tab[,"habAv_imp"] = Perid_habAv_vec
Perid_output_tab[,"unw_b_c_imp"] = Perid_unw_b_c_vec
Perid_output_tab[,"Patch_Area_imp"] = Perid_Patch_Area_vec


# Make a dataframe to make plots for comparison
Perid_HSI_tab = data.frame(imp = Perid_HSI_vec)
Perid_EgoSize_tab = data.frame(imp = Perid_EgoSize_vec)
Perid_strength_tab = data.frame(imp = Perid_strength_vec)
Perid_deg_tab = data.frame(imp = Perid_deg_vec)
Perid_habAv_tab = data.frame(imp = Perid_habAv_vec)
Perid_unw_b_c_tab = data.frame(imp = Perid_unw_b_c_vec)
Perid_Patch_Area_tab = data.frame(imp = Perid_Patch_Area_vec)

Perid_HSI_tab$variable = "HSI"
Perid_EgoSize_tab$variable = "3rd. ord. neigh."
Perid_strength_tab$variable = "Strength"
Perid_deg_tab$variable = "Degree"
Perid_habAv_tab$variable = "Hab. Av."
Perid_unw_b_c_tab$variable = "B.C."
Perid_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Perid_AUC_cv_tab = data.frame(value = Perid_AUC_cv_vec)
Perid_AUC_train_tab = data.frame(value = Perid_AUC_vec)
#Make label of model for plot
Perid_AUC_cv_tab$model = "Perid"
Perid_AUC_train_tab$model = "Perid"

#combine pred. vars. into new data frame 
Perid_var_imp_tab = rbind(Perid_HSI_tab,Perid_EgoSize_tab,Perid_strength_tab,Perid_habAv_tab,Perid_deg_tab,Perid_unw_b_c_tab,Perid_Patch_Area_tab)

ggplot(Perid_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Perid_var_imp_tab$imp~Perid_var_imp_tab$variable, ylab= "Var. Importance", main = "Perid")

#Get mean var. importance of all of the vars. 
mean(Perid_HSI_vec)
mean(Perid_EgoSize_vec)
mean(Perid_strength_vec)
mean(Perid_deg_vec)
mean(Perid_habAv_vec)
mean(Perid_unw_b_c_vec)
mean(Perid_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Perid_output_tab$AUC_cv)
summary(Perid_output_tab$AUC_train)



############ NoTopo (Model with only non-topological predictors included) #################################

# Create the output table that contains all the values of each of the runs
noTopo_output_tab = data.frame(run_nr = c(1:n_repeats)) #Perid
# Create vectors to record performance measures
noTopo_AUC_cv_vec = vector() #Cross-validated AUC
noTopo_AUC_vec = vector() #Training AUC
noTopo_nt_vec = vector()# number of trees
# noTopo_predict_vec = vector()

#Create vectors with all the values of a certain predictor along all the runs
noTopo_HSI_vec = vector() 
noTopo_Patch_Area_vec  = vector()

#Loop that goes exactly for 100 iterations, to get distributions 

for(i in c(1:n_repeats)){
  
  #Perform gbm with a fixed number of trees and no cross-validation.
  gbm_mod_noTopo = gbm.step(data=Perid_4Vt_woNA, gbm.x = c('Patch_Area','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_noTopo)$var, imp = summary(gbm_mod_noTopo)$rel.inf)
  
  #Write var_imp results to the vectors
  noTopo_HSI_vec = append(noTopo_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  noTopo_Patch_Area_vec  = append(noTopo_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the CV-AUC of this model to a vector
  noTopo_AUC_cv_vec = append(noTopo_AUC_cv_vec, gbm_mod_noTopo$cv.statistics$discrimination.mean)
  noTopo_AUC_vec = append(noTopo_AUC_vec, gbm_mod_noTopo$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_noTopo$n.trees
  noTopo_nt_vec = append(noTopo_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  # noTopo_predict_vec = append(predict(gbm_mod_noTopo, Perid_stattest, gbm_mod_noTopo$n.trees, type = "response", single.tree = FALSE))
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
  
}


noTopo_output_tab[,"AUC_train"] = noTopo_AUC_vec
noTopo_output_tab[,"AUC_cv"] = noTopo_AUC_cv_vec
noTopo_output_tab[,"ntrees"] = noTopo_nt_vec
#Var. importance columns
noTopo_output_tab[,"HSI_imp"] = noTopo_HSI_vec
noTopo_output_tab[,"Patch_Area_imp"] = noTopo_Patch_Area_vec


# Make a dataframe to make plots for comparison
noTopo_HSI_tab = data.frame(imp = noTopo_HSI_vec)
noTopo_Patch_Area_tab = data.frame(imp = noTopo_Patch_Area_vec)

noTopo_HSI_tab$variable = "HSI"
noTopo_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
noTopo_AUC_cv_tab = data.frame(value = noTopo_AUC_cv_vec)
noTopo_AUC_train_tab = data.frame(value = noTopo_AUC_vec)

#Make label of model for plot
noTopo_AUC_cv_tab$model = "noTopo"
noTopo_AUC_train_tab$model = "noTopo"

#combine pred. vars. into new data frame 
noTopo_var_imp_tab = rbind(noTopo_HSI_tab,noTopo_Patch_Area_tab)

ggplot(noTopo_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(noTopo_var_imp_tab$imp~noTopo_var_imp_tab$variable, ylab= "Var. Importance", main = "noTopo_Perid")

#Get mean var. importance of all of the vars. 
mean(noTopo_HSI_vec)
mean(noTopo_Patch_Area_vec)


##Check distr. of measures of prediction accuracy 
summary(noTopo_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_cv)



##### Compare between runs ###########################

#cv AUC
#mean
AUC_cv_tab = rbind(Perid_AUC_cv_tab, noTopo_AUC_cv_tab)

#Change order of factors to display noTopo at the edge
AUC_cv_tab$model<- as.factor(AUC_cv_tab$model)
levels(AUC_cv_tab$model)
AUC_cv_tab$model<-factor(AUC_cv_tab$model, levels=c("Perid", "noTopo"))
#print out
levels(AUC_cv_tab$model)

#Plot
ggplot(AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_cv_tab$value~AUC_cv_tab$model, ylab= "Cross-validated AUC")


#training AUC
AUC_train_tab = rbind(Perid_AUC_train_tab, noTopo_AUC_train_tab)

ggplot(AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_train_tab$value~AUC_train_tab$model, ylab= "Training AUC")


summary(Perid_output_tab$AUC_cv)
summary(noTopo_output_tab$AUC_cv)

summary(Perid_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_train)


summary(Perid_output_tab$ntrees)
summary(noTopo_output_tab$ntrees)

###Get avg_#trees over the AUC>=.75 threshold
Perid_mean_acceptable_ntrees <- mean(Perid_output_tab$ntrees[Perid_output_tab$AUC_cv>=0.75])
Perid_mean_acceptable_ntrees
noTopo_mean_acceptable_ntrees <- mean(noTopo_output_tab$ntrees[noTopo_output_tab$AUC_cv>=0.75])
noTopo_mean_acceptable_ntrees
Total_mean_acceptable_ntrees <- mean(Perid_mean_acceptable_ntrees,noTopo_mean_acceptable_ntrees)
Total_mean_acceptable_ntrees



#### Get sample models for response curves (tagged with the kind of network they are) ####
#gbm.mod (No tag) = noTopo model (as it was the last one to be run)

#NoTopo
test = gbm.plot(gbm_mod_noTopo,return.grid=TRUE)
str(test)
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)
#Perid
gbm_mod_Perid = gbm.step(data=Perid_4Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','habAv','EgoSize','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

#plot
gbm.plot(gbm_mod_Perid, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,3), write.title=FALSE, cex.axis = 2,  cex.lab=2)
#NoTopo
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)


### specific stats of sample models:
summary(gbm_mod_Perid)
gbm_mod_Perid$self.statistics$discrimination

#noTopo
gbm_mod_noTopo$self.statistics$discrimination

write.csv(Perid_4Vt_woNA, file = "BinPredRun_Perid_4Vt_woNA_b_c_unw.csv")




# #################################################################################################
# #### Get discrete and rescaled predictions ####
# #################################################################################################
# #### Add continuous predictions to dataframe w/ all the predictors & occurrence-state ####
# #Do new dataframe to keep the original unchanged (only add underscore between gender & sp.)
# Pe_rid_stattest <- Perid_stattest
# 
# 
# Perid_predict_df = data.frame(Perid_predict_mat)
# #Calculate the mean of the 100 models for each row (patch)
# # For past round it was done as uniform_stattest$cont_pred <- Singlemodel_Hy_arb_predict
# Perid_predict_df$mean_pred <- rowMeans(Perid_predict_mat, na.rm = FALSE)
# 
# Perid_predict_df$PatchID <- Pe_rid_stattest$PatchID
# 
# Perid_predict_df <- merge(Pe_rid_stattest, Perid_predict_df, by = "PatchID")
# names(Perid_predict_df)
# 
# ## Set discrete prediction to patches
# #Set discretization (binariation) cutpoint for each patch
# cp = cutpointr(Perid_predict_df, mean_pred, pres_abs, subgroup = NULL, method = minimize_metric, metric = roc01, na.rm = TRUE)
# 
# summary(cp)
# plot(cp)
# cp$optimal_cutpoint
# #Set all above threshold to 1, all below to 0 (discrete = disc)
# Perid_predict_df$pred_pres_abs <- ifelse(Perid_predict_df$mean_pred>=cp$optimal_cutpoint, 1, 0)
# 
# table(Perid_predict_df$pred_pres_abs)
# 
# #Add field w/ BRT Prediction score multiplied by 1000 (Mean_occurrence_index)
# #rescaled continuous prediction, to correspond to scale of HSI
# Perid_predict_df$Pred_1000 <- Perid_predict_df$mean_pred*1000
# 
# 
# #Export df to csv without the columns used only for calculation of binary prediction
# Pe_rid_stattest$mean_pred <- Perid_predict_df$mean_pred
# Pe_rid_stattest$Pred_1000 <-Perid_predict_df$Pred_1000
# Pe_rid_stattest$bin_OLpred <- Perid_predict_df$pred_pres_abs


### For ArcGIS importing (compatibility of fields)
Pe_rid_stattest$Value <- Pe_rid_stattest$PatchID
Pe_rid_stattest$V_s <- Pe_rid_stattest$V_Pe_rid
head(Pe_rid_stattest)

Pe_rid_stattest <- moveMe(Pe_rid_stattest, c("Value"), "after", "PatchID")
Pe_rid_stattest <- moveMe(Pe_rid_stattest, c("V_s"), "after", "V_Pe_rid")

# setwd("C:/Users/damiano/Documents/PhD/Additional_species_runs/Discrete_predictions")
# write.csv(Pe_rid_stattest, file = "Pe_rid_bin_prediction.csv")


#Hard-save evaluation measures
setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
write.csv(Perid_AUC_cv_tab, file = "Perid_DefaultDispDist_AUC_cv_tab.csv")
write.csv(Perid_AUC_train_tab, file = "Perid_DefaultDispDist_AUC_train_tab.csv")
write.csv(Perid_var_imp_tab, file = "Perid_DefaultDispDist_var_imp_tab.csv")
write.csv(Perid_output_tab, file = "Perid_DefaultDispDist_BRToutput_tab.csv")
write.csv(Perid_predict_df, file = "Perid_DefaultDispDist_predict_df.csv")

#Hard-save variable importance of each predictor
write.csv(Perid_HSI_tab, file = "Perid_DefaultDispDist_HSI_tab.csv")
write.csv(Perid_habAv_tab, file = "Perid_DefaultDispDist_habAv_tab.csv")
write.csv(Perid_deg_tab, file = "Perid_DefaultDispDist_deg_tab.csv")
write.csv(Perid_EgoSize_tab, file = "Perid_DefaultDispDist_EgoSize_tab.csv")
write.csv(Perid_unw_b_c_tab, file = "Perid_DefaultDispDist_unw_b_c_tab.csv")
write.csv(Perid_strength_tab, file = "Perid_DefaultDispDist_strength_tab.csv")
write.csv(Perid_Patch_Area_tab, file = "Perid_DefaultDispDist_Patch_Area_tab.csv")

#Hard-save model to be  able to use it in other workspaces
saveRDS(gbm_mod_Perid, file = "gbm_mod_Perid_DefaultDispDist.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# readRDS(file, refhook = NULL)
# infoRDS(file)

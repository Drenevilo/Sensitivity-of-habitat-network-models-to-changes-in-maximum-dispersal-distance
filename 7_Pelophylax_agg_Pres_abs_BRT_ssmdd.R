#### Script for: ##################################################################################################
#### Definition of Absences #### #### Calculation of predictors #### Boosted Regression Trees (BRT's) Modelling ###
#### For Pelophylax agg #### #################################################################################
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

Macrotable_Pe_agg <- read.csv("Macrotable_Peagg.csv")
head(Macrotable_Pe_agg)
names(Macrotable_Pe_agg)


#Drop OID  
Macrotable_Pe_agg$OID <- NULL

#Reorder columns, to be consistent with species alphabetical order
Macrotable_Pe_agg <- Macrotable_Pe_agg[,order(colnames(Macrotable_Pe_agg))] 
head(Macrotable_Pe_agg)

#########Get binary values from the counts ####
#ifelse(df>0, 1, 0 ) 
#converts it to a matrix
Pres_abs_allpatches_Pe_agg <- ifelse(Macrotable_Pe_agg >0, 1, 0 )
#recover data.frame type
Pres_abs_allpatches_Pe_agg <- as.data.frame(Pres_abs_allpatches_Pe_agg, row.names = NULL, optional = FALSE)
names(Pres_abs_allpatches_Pe_agg)
head(Pres_abs_allpatches_Pe_agg)
#Drop binarized value & count
Pres_abs_allpatches_Pe_agg$Value <- NULL
Pres_abs_allpatches_Pe_agg$Count <- NULL

#Recovering the identifier columns
Pres_abs_allpatches_Pe_agg$Value <- Macrotable_Pe_agg$Value
Pres_abs_allpatches_Pe_agg$Count <- Macrotable_Pe_agg$Count
names(Pres_abs_allpatches_Pe_agg)

#Get fields VALUE & COUNT to the beginning
Pres_abs_allpatches_Pe_agg <- moveMe(Pres_abs_allpatches_Pe_agg, c("Value", "Count"), "first")
#head(Pres_abs_allpatches_Pe_agg)
names(Pres_abs_allpatches_Pe_agg)


#### Get V_sp - Observations of Gender_species ####
#(sum(1 or 0 each year) 
Pres_abs_allpatches_Pe_agg$V_Pe_agg <- rowSums(Pres_abs_allpatches_Pe_agg[grep("Pe_agg", names(Pres_abs_allpatches_Pe_agg))])
names(Pres_abs_allpatches_Pe_agg)


#### Calculate observations of any amphibian per year ####

for(year in c("06","07","08","09","10","11","12","13","14","15")){
  Pres_abs_allpatches_Pe_agg[,paste("V_", year, sep = "")]<- rowSums(Pres_abs_allpatches_Pe_agg[grep(year, names(Pres_abs_allpatches_Pe_agg))])
}

Pres_abs_allpatches_Pe_agg$V_06 <- ifelse(Pres_abs_allpatches_Pe_agg$V_06>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_07 <- ifelse(Pres_abs_allpatches_Pe_agg$V_07>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_08 <- ifelse(Pres_abs_allpatches_Pe_agg$V_08>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_09 <- ifelse(Pres_abs_allpatches_Pe_agg$V_09>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_10 <- ifelse(Pres_abs_allpatches_Pe_agg$V_10>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_11 <- ifelse(Pres_abs_allpatches_Pe_agg$V_11>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_12 <- ifelse(Pres_abs_allpatches_Pe_agg$V_12>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_13 <- ifelse(Pres_abs_allpatches_Pe_agg$V_13>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_14 <- ifelse(Pres_abs_allpatches_Pe_agg$V_14>0, 1, 0)
Pres_abs_allpatches_Pe_agg$V_15 <- ifelse(Pres_abs_allpatches_Pe_agg$V_15>0, 1, 0)


### Get V_t - Total observations of amphibians
#By summing the vectors of spp per year
Pres_abs_allpatches_Pe_agg$V_t <- rowSums(Pres_abs_allpatches_Pe_agg[,c("V_06","V_07","V_08","V_09","V_10","V_11","V_12","V_13","V_14","V_15")])


################

#Plot Vt vs Vh
boxplot(Pres_abs_allpatches_Pe_agg$V_Pe_agg ~ Pres_abs_allpatches_Pe_agg$V_t, xlab= "Vt", ylab = "Vs")
#plot bars indicating median 

write.csv(Pres_abs_allpatches_Pe_agg, file = "Pres_abs_Vs_Vt_Pe_agg.csv") 

#Total number of patches
length(Pres_abs_allpatches_Pe_agg$V_t)


################ Occurence-state (Pres_abs) decider Plots #########################

#Pres_abs_allpatches2 is the one modified to make plots, excludes all the NA 
#Pres_abs_allpatches (unchanged, with all the 0's) will be used to incorporate the graph properties
#Get only the records that are visits (Vt >= 1)
Pres_abs_allpatches_Pe_agg2 = Pres_abs_allpatches_Pe_agg[Pres_abs_allpatches_Pe_agg$V_t != 0,] 

#Bar graph with proportions, averages and text --> The main one to decide thresholds!!!

#Data for the plot #Taking mean Vs
tab = ddply(Pres_abs_allpatches_Pe_agg2, .(V_t), summarize,  mean=mean(V_Pe_agg), countH = length(which(V_Pe_agg>0)), countT = length(V_Pe_agg), lab = paste(length(which(V_Pe_agg>0)),"/",length(V_Pe_agg),sep = " "), perc = length(which(V_Pe_agg>0))/length(V_Pe_agg)*100)

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
mean(Pres_abs_allpatches_Pe_agg$V_Pe_agg[Pres_abs_allpatches_Pe_agg$V_t==3])

#Number of absences at threshold value:
sum(Pres_abs_allpatches_Pe_agg$V_t>2 & Pres_abs_allpatches_Pe_agg$V_Pe_agg<1) 

which(Pres_abs_allpatches_Pe_agg$V_t >= 3 & Pres_abs_allpatches_Pe_agg$V_Pe_agg ==0)


#rename(data, d = b)
Pres_abs_allpatches_Pe_agg <- setnames(Pres_abs_allpatches_Pe_agg, "Value", "PatchID")
Pres_abs_allpatches_Pe_agg <- setnames(Pres_abs_allpatches_Pe_agg, "Count", "Patch_Area")


#### Add HSI ####
HSI = read.csv("C:/Users/damiano/Documents/PhD/Additional_species_runs/Network_setup/HSI_habPatchCode_Peagg1.csv")
Pres_abs_allpatches_Pe_agg$HSI <- HSI$"MEAN"



######## Importing network ############

setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Networks')

#Import table of patches w/network properties
Peagg_graph <- read_graph("Peagg.graphml", format = "graphml") 

#Look at vertices and edges of a graph
E(Peagg_graph)
V(Peagg_graph)


#### get metrics as attributes and into a data.frame ####
Peagg_graph_metrics <- data.frame(
  PatchID = V(Peagg_graph)$name,
  deg=degree(Peagg_graph),
  strength = strength(Peagg_graph, vids = V(Peagg_graph), mode = "all", loops = FALSE, weights = E(Peagg_graph)$weight),
  EgoSize = ego_size(Peagg_graph, order = 3, nodes = V(Peagg_graph))
  
)
### EgoSize = Third-order neighborhood


#### join the tables by attribute
#inner join= join by common variable names
patches_Peagg_topo_attributes <- merge(Pres_abs_allpatches_Pe_agg, Peagg_graph_metrics, by = "PatchID")


#### Add Habitat Availability predictor ####
Peagg_habAv <- scan("Peagg_habAv.txt")
patches_Peagg_topo_attributes$habAv <- Peagg_habAv
head(patches_Peagg_topo_attributes$habAv)


#### Measure betweenness centrality without weight included ####
#Re-import graphs, change name of 'weight' attribute
#setwd('/Networks')
Peagg_graph2 <- read_graph("Peagg.graphml", format = "graphml")

#rename weight attribute
#Actually copying into a new attr. and deleting original attribute
E(Peagg_graph2)$Alter_wght <- E(Peagg_graph2)$weight
Peagg_graph2 <- remove.edge.attribute(Peagg_graph2, "weight")

#Calculate explicitely unweighted b_c
#Set dataframe w/only PatchID & unw_b_c
Peagg_unweighted_b_c <- data.frame(
  PatchID = V(Peagg_graph2)$name,
  b_c=betweenness(Peagg_graph2, weights = NULL))

#rename unw_b_c
Peagg_unweighted_b_c <- rename(Peagg_unweighted_b_c,  "unw_b_c" = "b_c")

#join the tables by attribute
patches_Peagg_topo_attributes <- merge(x = patches_Peagg_topo_attributes, y = Peagg_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(patches_Peagg_topo_attributes)
head(patches_Peagg_topo_attributes)



#### Define occurrence_state (pres_abs) ####
### Classify the dataframe in subsets that classify all the patches in presences, absences or questionmarks/NoIdea (NA) 

#subset of questionmarks, #All the records, except the ones that comply with the following commands
patches_Peagg_topo_attributes[,"pres_abs"] = NA 
#Subset of likely absences (With threshold defined by pres_abs decider plot))
patches_Peagg_topo_attributes[which(patches_Peagg_topo_attributes$V_t >= 3 & patches_Peagg_topo_attributes$V_Pe_agg ==0),"pres_abs"] = 0
#Subset of confirmed presences
patches_Peagg_topo_attributes[which(patches_Peagg_topo_attributes$V_Pe_agg > 0),"pres_abs"] = 1

#### Display the 'occurrence-state' (0/1/NA) of every patch as defined above (pres_abs) ####
patches_Peagg_topo_attributes$pres_abs

#Get total amount of values for presence and for absence
table(patches_Peagg_topo_attributes$pres_abs)


#### Do a subset of the columns that only has the predictors and response ####
Peagg_stattest <- subset(patches_Peagg_topo_attributes, select=c("PatchID", "Patch_Area", "V_Pe_agg", "V_t", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Peagg_stattest)
head(Peagg_stattest)

#change directory to the one with BRT things 
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/BRTs')
#write file
#t4 -> pres_abs threshold=4
write.csv(Peagg_stattest, file = "Peagg_stattest_t3.csv")


##### Do the 'only pres_abs (no-NoData) df's #########################
#Take out NA's, set threshold of absences (at V_t = 3 for Pe_agg), change name
Peagg_3Vt_woNA = Peagg_stattest
Peagg_3Vt_woNA[Peagg_3Vt_woNA$V_Pe_agg==0 & Peagg_3Vt_woNA$V_t >= 3, 'pres_abs'] = 0
Peagg_3Vt_woNA = Peagg_3Vt_woNA[!is.na(Peagg_3Vt_woNA$pres_abs),]

#check values of unw_b_c
mean(Peagg_3Vt_woNA$unw_b_c)



#### Perform BRT's ############################################################
#### w/only unw_b_c

##Run 100 models without gridsearch, see results, how variable they are
#Hyperparameters for all models
lr = 0.001
tc = 5
bf = 0.75

n_repeats = 100

##### Peagg (Full model with all the predictors) ###########################################

# Create the output table that contains all the values of each of the runs
Peagg_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg
# Create vectors to record performance measures
Peagg_AUC_cv_vec = vector() #Cross-validated AUC
Peagg_AUC_vec = vector() #Training AUC
Peagg_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Peagg_HSI_vec = vector() 
Peagg_EgoSize_vec = vector() 
Peagg_strength_vec  = vector()
Peagg_deg_vec  = vector()
Peagg_habAv_vec  = vector()
Peagg_unw_b_c_vec  = vector()
Peagg_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Peagg_predict_mat = matrix(nrow = length(Peagg_stattest$PatchID), ncol = n_repeats, byrow = FALSE)


### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Peagg = gbm.step(data=Peagg_3Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Peagg)$var, imp = summary(gbm_mod_Peagg)$rel.inf)
  
  #Write var_imp results to the vectors
  Peagg_HSI_vec = append(Peagg_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Peagg_EgoSize_vec = append(Peagg_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Peagg_strength_vec  = append(Peagg_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Peagg_deg_vec  = append(Peagg_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Peagg_habAv_vec  = append(Peagg_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Peagg_unw_b_c_vec  = append(Peagg_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Peagg_Patch_Area_vec  = append(Peagg_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Peagg_AUC_cv_vec = append(Peagg_AUC_cv_vec, gbm_mod_Peagg$cv.statistics$discrimination.mean)
  Peagg_AUC_vec = append(Peagg_AUC_vec, gbm_mod_Peagg$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Peagg$n.trees
  Peagg_nt_vec = append(Peagg_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Peagg_predict_mat[,i] = predict(gbm_mod_Peagg, Peagg_stattest, gbm_mod_Peagg$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}


Peagg_output_tab[,"AUC_train"] = Peagg_AUC_vec
Peagg_output_tab[,"AUC_cv"] = Peagg_AUC_cv_vec
Peagg_output_tab[,"ntrees"] = Peagg_nt_vec


#Var. importance columns
Peagg_output_tab[,"HSI_imp"] = Peagg_HSI_vec
Peagg_output_tab[,"EgoSize_imp"] = Peagg_EgoSize_vec
Peagg_output_tab[,"strength_imp"] = Peagg_strength_vec
Peagg_output_tab[,"deg_imp"] = Peagg_deg_vec
Peagg_output_tab[,"habAv_imp"] = Peagg_habAv_vec
Peagg_output_tab[,"unw_b_c_imp"] = Peagg_unw_b_c_vec
Peagg_output_tab[,"Patch_Area_imp"] = Peagg_Patch_Area_vec


# Make a dataframe to make plots for comparison
Peagg_HSI_tab = data.frame(imp = Peagg_HSI_vec)
Peagg_EgoSize_tab = data.frame(imp = Peagg_EgoSize_vec)
Peagg_strength_tab = data.frame(imp = Peagg_strength_vec)
Peagg_deg_tab = data.frame(imp = Peagg_deg_vec)
Peagg_habAv_tab = data.frame(imp = Peagg_habAv_vec)
Peagg_unw_b_c_tab = data.frame(imp = Peagg_unw_b_c_vec)
Peagg_Patch_Area_tab = data.frame(imp = Peagg_Patch_Area_vec)

Peagg_HSI_tab$variable = "HSI"
Peagg_EgoSize_tab$variable = "3rd. ord. neigh."
Peagg_strength_tab$variable = "Strength"
Peagg_deg_tab$variable = "Degree"
Peagg_habAv_tab$variable = "Hab. Av."
Peagg_unw_b_c_tab$variable = "B.C."
Peagg_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to do overlaying histograms comparing performance between models
Peagg_AUC_cv_tab = data.frame(value = Peagg_AUC_cv_vec)
Peagg_AUC_train_tab = data.frame(value = Peagg_AUC_vec)
#Make label of model for plot
Peagg_AUC_cv_tab$model = "Peagg"
Peagg_AUC_train_tab$model = "Peagg"

#combine pred. vars. into new data frame 
Peagg_var_imp_tab = rbind(Peagg_HSI_tab,Peagg_EgoSize_tab,Peagg_strength_tab,Peagg_habAv_tab,Peagg_deg_tab,Peagg_unw_b_c_tab,Peagg_Patch_Area_tab)

ggplot(Peagg_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Peagg_var_imp_tab$imp~Peagg_var_imp_tab$variable, ylab= "Var. Importance", main = "Peagg")

#Get mean var. importance of all of the vars. 
mean(Peagg_HSI_vec)
mean(Peagg_EgoSize_vec)
mean(Peagg_strength_vec)
mean(Peagg_deg_vec)
mean(Peagg_habAv_vec)
mean(Peagg_unw_b_c_vec)
mean(Peagg_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Peagg_output_tab$AUC_cv)
summary(Peagg_output_tab$AUC_train)



############ NoTopo (Model with only non-topological predictors included) #################################

# Create the output table that contains all the values of each of the runs
noTopo_output_tab = data.frame(run_nr = c(1:n_repeats)) #Peagg
# Create vectors to record performance measures
noTopo_AUC_cv_vec = vector() #Cross-validated AUC
noTopo_AUC_vec = vector() #Training AUC
noTopo_nt_vec = vector() # number of trees
# noTopo_predict_vec = vector()

#Create vectors with all the values of a certain predictor along all the runs
noTopo_HSI_vec = vector() 
noTopo_Patch_Area_vec  = vector()

#Loop that goes exactly for 100 iterations, to get distributions 

for(i in c(1:n_repeats)){

  #Perform gbm with a fixed number of trees and no cross-validation
  gbm_mod_noTopo = gbm.step(data=Peagg_3Vt_woNA, gbm.x = c('Patch_Area','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
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
  # noTopo_predict_vec = append(predict(gbm_mod_noTopo, Peagg_stattest, gbm_mod_noTopo$n.trees, type = "response", single.tree = FALSE))
  
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

#Reserve also measures in df to to compare performance between models
noTopo_AUC_cv_tab = data.frame(value = noTopo_AUC_cv_vec)
noTopo_AUC_train_tab = data.frame(value = noTopo_AUC_vec)

#Make label of model for plot
noTopo_AUC_cv_tab$model = "noTopo"
noTopo_AUC_train_tab$model = "noTopo"

#combine pred. vars. into new data frame 
noTopo_var_imp_tab = rbind(noTopo_HSI_tab,noTopo_Patch_Area_tab)

ggplot(noTopo_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(noTopo_var_imp_tab$imp~noTopo_var_imp_tab$variable, ylab= "Var. Importance", main = "noTopo_Peagg")

#Get mean var. importance of all of the vars. 
mean(noTopo_HSI_vec)
mean(noTopo_Patch_Area_vec)


##Check distr. of measures of prediction accuracy 
summary(noTopo_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_cv)



##### Compare between runs ###########################

#cv AUC
#mean, then median
AUC_cv_tab = rbind(Peagg_AUC_cv_tab, noTopo_AUC_cv_tab)

#Change order of factors to display noTopo at the edge
AUC_cv_tab$model<- as.factor(AUC_cv_tab$model)
levels(AUC_cv_tab$model)
AUC_cv_tab$model<-factor(AUC_cv_tab$model, levels=c("Peagg", "noTopo"))
#print out
levels(AUC_cv_tab$model)

#Plot
ggplot(AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_cv_tab$value~AUC_cv_tab$model, ylab= "Cross-validated AUC")


#training AUC
AUC_train_tab = rbind(Peagg_AUC_train_tab, noTopo_AUC_train_tab)

ggplot(AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_train_tab$value~AUC_train_tab$model, ylab= "Training AUC")


summary(Peagg_output_tab$AUC_cv)
summary(noTopo_output_tab$AUC_cv)

summary(Peagg_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_train)

summary(Peagg_output_tab$ntrees)
summary(noTopo_output_tab$ntrees)


###Get avg_#trees over the AUC>=.75 threshold
Peagg_mean_acceptable_ntrees <- mean(Peagg_output_tab$ntrees[Peagg_output_tab$AUC_cv>=0.75])
Peagg_mean_acceptable_ntrees
noTopo_mean_acceptable_ntrees <- mean(noTopo_output_tab$ntrees[noTopo_output_tab$AUC_cv>=0.75])
noTopo_mean_acceptable_ntrees
Total_mean_acceptable_ntrees <- mean(Peagg_mean_acceptable_ntrees,noTopo_mean_acceptable_ntrees)
Total_mean_acceptable_ntrees



#### Get sample models for response curves (tagged with the kind of network they are) ####
#gbm.mod (No tag) = noTopo model (as it was the last one to be run)

#NoTopo
test = gbm.plot(gbm_mod_noTopo,return.grid=TRUE)
str(test)
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)
#Peagg
gbm_mod_Peagg = gbm.step(data=Peagg_3Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','habAv','EgoSize','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

#plot
gbm.plot(gbm_mod_Peagg, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,3), write.title=FALSE, cex.axis = 2,  cex.lab=2)
#NoTopo
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)


### specific stats of sample models:
summary(gbm_mod_Peagg)
gbm_mod_Peagg$self.statistics$discrimination

#noTopo
gbm_mod_noTopo$self.statistics$discrimination

write.csv(Peagg_3Vt_woNA, file = "BinPredRun_Peagg_3Vt_woNA_b_c_unw.csv")



# #################################################################################################
# #### Get discrete and rescaled predictions ####
# #################################################################################################
# #### Add continuous predictions to dataframe w/ all the predictors & occurrence-state ####
# #Do new dataframe to keep the original unchanged (only add underscore between gender & sp.)
# Pe_agg_stattest <- Peagg_stattest
# 
# 
# Peagg_predict_df = data.frame(Peagg_predict_mat)
# #Calculate the mean of the 100 models for each row (patch)
# # For past round it was done as uniform_stattest$cont_pred <- Singlemodel_Hy_arb_predict
# Peagg_predict_df$mean_pred <- rowMeans(Peagg_predict_mat, na.rm = FALSE)
# 
# Peagg_predict_df$PatchID <- Pe_agg_stattest$PatchID
# 
# Peagg_predict_df <- merge(Pe_agg_stattest, Peagg_predict_df, by = "PatchID")
# names(Peagg_predict_df)
# 
# ## Set discrete prediction to patches
# #Set discretization (binariation) cutpoint for each patch
# cp = cutpointr(Peagg_predict_df, mean_pred, pres_abs, subgroup = NULL, method = minimize_metric, metric = roc01, na.rm = TRUE)
# 
# summary(cp)
# plot(cp)
# cp$optimal_cutpoint
# #Set all above threshold to 1, all below to 0 (discrete = disc)
# Peagg_predict_df$pred_pres_abs <- ifelse(Peagg_predict_df$mean_pred>=cp$optimal_cutpoint, 1, 0)
# 
# table(Peagg_predict_df$pred_pres_abs)
# 
# #Add field w/ BRT Prediction score multiplied by 1000 (Mean_occurrence_index)
# #rescaled continuous prediction, to correspond to scale of HSI
# Peagg_predict_df$Pred_1000 <- Peagg_predict_df$mean_pred*1000
# 
# 
# #Export df to csv without the columns used only for calculation of binary prediction
# Pe_agg_stattest$mean_pred <- Peagg_predict_df$mean_pred
# Pe_agg_stattest$Pred_1000 <-Peagg_predict_df$Pred_1000
# Pe_agg_stattest$bin_OLpred <- Peagg_predict_df$pred_pres_abs

### For ArcGIS importing (compatibility of fields)
Pe_agg_stattest$Value <- Pe_agg_stattest$PatchID
Pe_agg_stattest$V_s <- Pe_agg_stattest$V_Pe_agg
head(Pe_agg_stattest)

Pe_agg_stattest <- moveMe(Pe_agg_stattest, c("Value"), "after", "PatchID")
Pe_agg_stattest <- moveMe(Pe_agg_stattest, c("V_s"), "after", "V_Pe_agg")

# setwd("C:/Users/damiano/Documents/PhD/Additional_species_runs/Discrete_predictions")
# write.csv(Pe_agg_stattest, file = "Pe_agg_bin_prediction.csv")


#Hard-save evaluation measures
setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
write.csv(Peagg_AUC_cv_tab, file = "Peagg_DefaultDispDist_AUC_cv_tab.csv")
write.csv(Peagg_AUC_train_tab, file = "Peagg_DefaultDispDist_AUC_train_tab.csv")
write.csv(Peagg_var_imp_tab, file = "Peagg_DefaultDispDist_var_imp_tab.csv")
write.csv(Peagg_output_tab, file = "Peagg_DefaultDispDist_BRToutput_tab.csv")
write.csv(Peagg_predict_df, file = "Peagg_DefaultDispDist_predict_df.csv")

#Hard-save variable importance of each predictor
write.csv(Peagg_HSI_tab, file = "Peagg_DefaultDispDist_HSI_tab.csv")
write.csv(Peagg_habAv_tab, file = "Peagg_DefaultDispDist_habAv_tab.csv")
write.csv(Peagg_deg_tab, file = "Peagg_DefaultDispDist_deg_tab.csv")
write.csv(Peagg_EgoSize_tab, file = "Peagg_DefaultDispDist_EgoSize_tab.csv")
write.csv(Peagg_unw_b_c_tab, file = "Peagg_DefaultDispDist_unw_b_c_tab.csv")
write.csv(Peagg_strength_tab, file = "Peagg_DefaultDispDist_strength_tab.csv")
write.csv(Peagg_Patch_Area_tab, file = "Peagg_DefaultDispDist_Patch_Area_tab.csv")

#Hard-save model to be  able to use it in other workspaces
saveRDS(gbm_mod_Peagg, file = "gbm_mod_Peagg_DefaultDispDist.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# readRDS(file, refhook = NULL)
# infoRDS(file)
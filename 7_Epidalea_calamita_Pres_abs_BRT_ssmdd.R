#### Script for: ##################################################################################################
#### Definition of Absences #### #### Calculation of predictors #### Boosted Regression Trees (BRT's) Modelling ###
#### For Epidalea calamita #### #################################################################################
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

Macrotable_Ep_cal <- read.csv("Macrotable_Epcal.csv")
head(Macrotable_Ep_cal)
names(Macrotable_Ep_cal)


#Drop OID  
Macrotable_Ep_cal$OID <- NULL

#Reorder columns, to be consistent with species alphabetical order
Macrotable_Ep_cal <- Macrotable_Ep_cal[,order(colnames(Macrotable_Ep_cal))] 
head(Macrotable_Ep_cal)

#########Get binary values from the counts ####
#ifelse(df>0, 1, 0 ) 
#converts it to a matrix
Pres_abs_allpatches_Ep_cal <- ifelse(Macrotable_Ep_cal >0, 1, 0 )
#recover data.frame type
Pres_abs_allpatches_Ep_cal <- as.data.frame(Pres_abs_allpatches_Ep_cal, row.names = NULL, optional = FALSE)
names(Pres_abs_allpatches_Ep_cal)
head(Pres_abs_allpatches_Ep_cal)
#Drop binarized value & count
Pres_abs_allpatches_Ep_cal$Value <- NULL
Pres_abs_allpatches_Ep_cal$Count <- NULL

#Recovering the identifier columns
Pres_abs_allpatches_Ep_cal$Value <- Macrotable_Ep_cal$Value
Pres_abs_allpatches_Ep_cal$Count <- Macrotable_Ep_cal$Count
names(Pres_abs_allpatches_Ep_cal)

#Get fields VALUE & COUNT to the beginning
Pres_abs_allpatches_Ep_cal <- moveMe(Pres_abs_allpatches_Ep_cal, c("Value", "Count"), "first")
#head(Pres_abs_allpatches_Ep_cal)
names(Pres_abs_allpatches_Ep_cal)


#### Get V_sp - Observations of Gender_species ####
#(sum(1 or 0 each year) 
Pres_abs_allpatches_Ep_cal$V_Ep_cal <- rowSums(Pres_abs_allpatches_Ep_cal[grep("Ep_cal", names(Pres_abs_allpatches_Ep_cal))])
names(Pres_abs_allpatches_Ep_cal)


#### Calculate observations of any amphibian per year ####

for(year in c("06","07","08","09","10","11","12","13","14","15")){
  Pres_abs_allpatches_Ep_cal[,paste("V_", year, sep = "")]<- rowSums(Pres_abs_allpatches_Ep_cal[grep(year, names(Pres_abs_allpatches_Ep_cal))])
}

Pres_abs_allpatches_Ep_cal$V_06 <- ifelse(Pres_abs_allpatches_Ep_cal$V_06>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_07 <- ifelse(Pres_abs_allpatches_Ep_cal$V_07>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_08 <- ifelse(Pres_abs_allpatches_Ep_cal$V_08>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_09 <- ifelse(Pres_abs_allpatches_Ep_cal$V_09>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_10 <- ifelse(Pres_abs_allpatches_Ep_cal$V_10>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_11 <- ifelse(Pres_abs_allpatches_Ep_cal$V_11>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_12 <- ifelse(Pres_abs_allpatches_Ep_cal$V_12>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_13 <- ifelse(Pres_abs_allpatches_Ep_cal$V_13>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_14 <- ifelse(Pres_abs_allpatches_Ep_cal$V_14>0, 1, 0)
Pres_abs_allpatches_Ep_cal$V_15 <- ifelse(Pres_abs_allpatches_Ep_cal$V_15>0, 1, 0)


### Get V_t - Total observations of amphibians
#By summing the vectors of spp per year
Pres_abs_allpatches_Ep_cal$V_t <- rowSums(Pres_abs_allpatches_Ep_cal[,c("V_06","V_07","V_08","V_09","V_10","V_11","V_12","V_13","V_14","V_15")])


################

#Plot Vt vs Vh
boxplot(Pres_abs_allpatches_Ep_cal$V_Ep_cal ~ Pres_abs_allpatches_Ep_cal$V_t, xlab= "Vt", ylab = "Vs")
#plot bars indicating median 

write.csv(Pres_abs_allpatches_Ep_cal, file = "Pres_abs_Vs_Vt_Ep_cal.csv") 

#Total number of patches
length(Pres_abs_allpatches_Ep_cal$V_t)


################ Occurence-state (Pres_abs) decider Plots #########################

#Pres_abs_allpatches2 is the one modified to make plots, excludes all the NA 
#Pres_abs_allpatches (unchanged, with all the 0's) will be used to incorporate the graph properties
#Get only the records that are visits (Vt >= 1)
Pres_abs_allpatches_Ep_cal2 = Pres_abs_allpatches_Ep_cal[Pres_abs_allpatches_Ep_cal$V_t != 0,] 

#Bar graph with proportions, averages and text --> The main one to decide thresholds!!!

#Data for the plot #Taking mean Vs
tab = ddply(Pres_abs_allpatches_Ep_cal2, .(V_t), summarize,  mean=mean(V_Ep_cal), countH = length(which(V_Ep_cal>0)), countT = length(V_Ep_cal), lab = paste(length(which(V_Ep_cal>0)),"/",length(V_Ep_cal),sep = " "), perc = length(which(V_Ep_cal>0))/length(V_Ep_cal)*100)

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
# E.g. #Vt=6 -> Vh=1.8125 #Threshold value
# (See Ortiz-Rodriguez et al 2019)
mean(Pres_abs_allpatches_Ep_cal$V_Ep_cal[Pres_abs_allpatches_Ep_cal$V_t==4])

#Number of absences at threshold value:
sum(Pres_abs_allpatches_Ep_cal$V_t>3 & Pres_abs_allpatches_Ep_cal$V_Ep_cal<1) 

which(Pres_abs_allpatches_Ep_cal$V_t >= 4 & Pres_abs_allpatches_Ep_cal$V_Ep_cal ==0)


#### rename(data, d = b) ####
Pres_abs_allpatches_Ep_cal <- setnames(Pres_abs_allpatches_Ep_cal, "Value", "PatchID")
Pres_abs_allpatches_Ep_cal <- setnames(Pres_abs_allpatches_Ep_cal, "Count", "Patch_Area")


#### Add HSI ####
HSI = read.csv("C:/Users/damiano/Documents/PhD/Additional_species_runs/Network_setup/HSI_habPatchCode_Epcal1.csv")
Pres_abs_allpatches_Ep_cal$HSI <- HSI$"MEAN"



######## Importing network ############

setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Networks')

#Import table of patches w/network properties
Epcal_graph <- read_graph("Epcal.graphml", format = "graphml") 

#Look at vertices and edges of a graph
E(Epcal_graph)
V(Epcal_graph)


#### get metrics as attributes and into a data.frame ####
Epcal_graph_metrics <- data.frame(
  PatchID = V(Epcal_graph)$name,
  deg=degree(Epcal_graph),
  strength = strength(Epcal_graph, vids = V(Epcal_graph), mode = "all", loops = FALSE, weights = E(Epcal_graph)$weight),
  EgoSize = ego_size(Epcal_graph, order = 3, nodes = V(Epcal_graph))
  
)
### EgoSize = Third-order neighborhood


#### join the tables by attribute
#inner join= join by common variable names
patches_Epcal_topo_attributes <- merge(Pres_abs_allpatches_Ep_cal, Epcal_graph_metrics, by = "PatchID")


#### Add Habitat Availability predictor ####
Epcal_habAv <- scan("Epcal_habAv.txt")
patches_Epcal_topo_attributes$habAv <- Epcal_habAv
head(patches_Epcal_topo_attributes$habAv)


#### Measure betweenness centrality without weight included ####
#Re-import graphs, change name of 'weight' attribute
#setwd('/Networks')
Epcal_graph2 <- read_graph("Epcal.graphml", format = "graphml")

#rename weight attribute
#Actually copying into a new attr. and deleting original attribute
E(Epcal_graph2)$Alter_wght <- E(Epcal_graph2)$weight
Epcal_graph2 <- remove.edge.attribute(Epcal_graph2, "weight")

#Calculate explicitely unweighted b_c
#Set dataframe w/only PatchID & unw_b_c
Epcal_unweighted_b_c <- data.frame(
  PatchID = V(Epcal_graph2)$name,
  b_c=betweenness(Epcal_graph2, weights = NULL))

#rename unw_b_c
Epcal_unweighted_b_c <- rename(Epcal_unweighted_b_c,  "unw_b_c" = "b_c")

#join the tables by attribute
patches_Epcal_topo_attributes <- merge(x = patches_Epcal_topo_attributes, y = Epcal_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(patches_Epcal_topo_attributes)
head(patches_Epcal_topo_attributes)



#### Define occurrence_state (pres_abs) ####
### Classify the dataframe in subsets that classify all the patches in presences, absences or questionmarks/NoIdea (NA) 

#subset of questionmarks, #All the records, except the ones that comply with the following commands
patches_Epcal_topo_attributes[,"pres_abs"] = NA 
#Subset of likely absences (With threshold at 5 (defined by pres_abs decider plot))
patches_Epcal_topo_attributes[which(patches_Epcal_topo_attributes$V_t >= 4 & patches_Epcal_topo_attributes$V_Ep_cal ==0),"pres_abs"] = 0
#Subset of confirmed presences
patches_Epcal_topo_attributes[which(patches_Epcal_topo_attributes$V_Ep_cal > 0),"pres_abs"] = 1

#### Display the 'occurrence-state' (0/1/NA) of every patch as defined above (pres_abs) ####
patches_Epcal_topo_attributes$pres_abs

#Get total amount of values for presence and for absence
table(patches_Epcal_topo_attributes$pres_abs)


#### Do a subset of the columns that only has the predictors and response ####
Epcal_stattest <- subset(patches_Epcal_topo_attributes, select=c("PatchID", "Patch_Area", "V_Ep_cal", "V_t", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Epcal_stattest)
head(Epcal_stattest)

#change directory to the one with BRT things 
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/BRTs')
#write file
#t4 -> pres_abs threshold=4
write.csv(Epcal_stattest, file = "Epcal_stattest_t4.csv")


##### Do the 'only pres_abs (no-NoData) df's #########################
#Take out NA's, set threshold of absences (at V_t = 4 for Ep_cal), change name
Epcal_4Vt_woNA = Epcal_stattest
Epcal_4Vt_woNA[Epcal_4Vt_woNA$V_Ep_cal==0 & Epcal_4Vt_woNA$V_t >= 4, 'pres_abs'] = 0
Epcal_4Vt_woNA = Epcal_4Vt_woNA[!is.na(Epcal_4Vt_woNA$pres_abs),]

#check values of unw_b_c
mean(Epcal_4Vt_woNA$unw_b_c)



#### Perform BRT's ############################################################
#Hyperparameters for all models
lr = 0.001
tc = 5
bf = 0.75

n_repeats = 100

##### Epcal (Full model with all the predictors) ###########################################

# Create the output table that contains all the values of each of the runs
Epcal_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal
# Create vectors to record performance measures
Epcal_AUC_cv_vec = vector() #Cross-validated AUC
Epcal_AUC_vec = vector() #Training AUC

Epcal_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Epcal_HSI_vec = vector() 
Epcal_EgoSize_vec = vector() 
Epcal_strength_vec  = vector()
Epcal_deg_vec  = vector()
Epcal_habAv_vec  = vector()
Epcal_unw_b_c_vec  = vector()
Epcal_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Epcal_predict_mat = matrix(nrow = length(Epcal_stattest$PatchID), ncol = n_repeats, byrow = FALSE)


### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Epcal = gbm.step(data=Epcal_4Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Epcal)$var, imp = summary(gbm_mod_Epcal)$rel.inf)
  
  #Write var_imp results to the vectors
  Epcal_HSI_vec = append(Epcal_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Epcal_EgoSize_vec = append(Epcal_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Epcal_strength_vec  = append(Epcal_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Epcal_deg_vec  = append(Epcal_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Epcal_habAv_vec  = append(Epcal_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Epcal_unw_b_c_vec  = append(Epcal_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Epcal_Patch_Area_vec  = append(Epcal_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Epcal_AUC_cv_vec = append(Epcal_AUC_cv_vec, gbm_mod_Epcal$cv.statistics$discrimination.mean)
  Epcal_AUC_vec = append(Epcal_AUC_vec, gbm_mod_Epcal$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Epcal$n.trees
  Epcal_nt_vec = append(Epcal_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Epcal_predict_mat[,i] = predict(gbm_mod_Epcal, Epcal_stattest, gbm_mod_Epcal$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}


Epcal_output_tab[,"AUC_train"] = Epcal_AUC_vec
Epcal_output_tab[,"AUC_cv"] = Epcal_AUC_cv_vec
Epcal_output_tab[,"ntrees"] = Epcal_nt_vec


#Var. importance columns
Epcal_output_tab[,"HSI_imp"] = Epcal_HSI_vec
Epcal_output_tab[,"EgoSize_imp"] = Epcal_EgoSize_vec
Epcal_output_tab[,"strength_imp"] = Epcal_strength_vec
Epcal_output_tab[,"deg_imp"] = Epcal_deg_vec
Epcal_output_tab[,"habAv_imp"] = Epcal_habAv_vec
Epcal_output_tab[,"unw_b_c_imp"] = Epcal_unw_b_c_vec
Epcal_output_tab[,"Patch_Area_imp"] = Epcal_Patch_Area_vec


# Make a dataframe to make plots for comparison
Epcal_HSI_tab = data.frame(imp = Epcal_HSI_vec)
Epcal_EgoSize_tab = data.frame(imp = Epcal_EgoSize_vec)
Epcal_strength_tab = data.frame(imp = Epcal_strength_vec)
Epcal_deg_tab = data.frame(imp = Epcal_deg_vec)
Epcal_habAv_tab = data.frame(imp = Epcal_habAv_vec)
Epcal_unw_b_c_tab = data.frame(imp = Epcal_unw_b_c_vec)
Epcal_Patch_Area_tab = data.frame(imp = Epcal_Patch_Area_vec)

Epcal_HSI_tab$variable = "HSI"
Epcal_EgoSize_tab$variable = "3rd. ord. neigh."
Epcal_strength_tab$variable = "Strength"
Epcal_deg_tab$variable = "Degree"
Epcal_habAv_tab$variable = "Hab. Av."
Epcal_unw_b_c_tab$variable = "B.C."
Epcal_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to compare performance between models
Epcal_AUC_cv_tab = data.frame(value = Epcal_AUC_cv_vec)
Epcal_AUC_train_tab = data.frame(value = Epcal_AUC_vec)
#Make label of model for plot
Epcal_AUC_cv_tab$model = "Epcal"
Epcal_AUC_train_tab$model = "Epcal"

#combine pred. vars. into new data frame 
Epcal_var_imp_tab = rbind(Epcal_HSI_tab,Epcal_EgoSize_tab,Epcal_strength_tab,Epcal_habAv_tab,Epcal_deg_tab,Epcal_unw_b_c_tab,Epcal_Patch_Area_tab)

ggplot(Epcal_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Epcal_var_imp_tab$imp~Epcal_var_imp_tab$variable, ylab= "Var. Importance", main = "Epcal")

#Get mean var. importance of all of the vars. 
mean(Epcal_HSI_vec)
mean(Epcal_EgoSize_vec)
mean(Epcal_strength_vec)
mean(Epcal_deg_vec)
mean(Epcal_habAv_vec)
mean(Epcal_unw_b_c_vec)
mean(Epcal_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Epcal_output_tab$AUC_cv)
summary(Epcal_output_tab$AUC_train)



############ NoTopo (Model with only non-topological predictors included) #################################

# Create the output table that contains all the values of each of the runs
noTopo_output_tab = data.frame(run_nr = c(1:n_repeats)) #Epcal
# Create vectors to record performance measures
noTopo_AUC_cv_vec = vector() #Cross-validated AUC
noTopo_AUC_vec = vector() # Training AUC
noTopo_nt_vec = vector() # number of trees
# noTopo_predict_vec = vector()

#Create vectors with all the values of a certain predictor along all the runs
noTopo_HSI_vec = vector() 
noTopo_Patch_Area_vec  = vector()

#Loop that goes exactly for 100 iterations, to get distributions 

for(i in c(1:n_repeats)){

  #Perform gbm with a fixed number of trees and no cross-validation
  gbm_mod_noTopo = gbm.step(data=Epcal_4Vt_woNA, gbm.x = c('Patch_Area','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
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
  
  #write the continuous prediction over all the patches
  # noTopo_predict_vec = append(predict(gbm_mod_noTopo, Epcal_stattest, gbm_mod_noTopo$n.trees, type = "response", single.tree = FALSE))
  
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

#Reserve also measures in df to compare performance between models
noTopo_AUC_cv_tab = data.frame(value = noTopo_AUC_cv_vec)
noTopo_AUC_train_tab = data.frame(value = noTopo_AUC_vec)

#Make label of model for plot
noTopo_AUC_cv_tab$model = "noTopo"
noTopo_AUC_train_tab$model = "noTopo"

#combine pred. vars. into new data frame 
noTopo_var_imp_tab = rbind(noTopo_HSI_tab,noTopo_Patch_Area_tab)

ggplot(noTopo_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(noTopo_var_imp_tab$imp~noTopo_var_imp_tab$variable, ylab= "Var. Importance", main = "noTopo_Epcal")

#Get mean var. importance of all of the vars. 
mean(noTopo_HSI_vec)
mean(noTopo_Patch_Area_vec)


##Check distr. of measures of prediction accuracy 
summary(noTopo_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_cv)



##### Compare between runs ###########################

#cv AUC
#mean
AUC_cv_tab = rbind(Epcal_AUC_cv_tab, noTopo_AUC_cv_tab)

#Change order of factors to display noTopo at the edge
AUC_cv_tab$model<- as.factor(AUC_cv_tab$model)
levels(AUC_cv_tab$model)
AUC_cv_tab$model<-factor(AUC_cv_tab$model, levels=c("Epcal", "noTopo"))
#print out
levels(AUC_cv_tab$model)

#Plot
ggplot(AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_cv_tab$value~AUC_cv_tab$model, ylab= "Cross-validated AUC")


#training AUC
AUC_train_tab = rbind(Epcal_AUC_train_tab, noTopo_AUC_train_tab)

ggplot(AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_train_tab$value~AUC_train_tab$model, ylab= "Training AUC")


summary(Epcal_output_tab$AUC_cv)
summary(noTopo_output_tab$AUC_cv)

summary(Epcal_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_train)

summary(Epcal_output_tab$ntrees)
summary(noTopo_output_tab$ntrees)


###Get avg_#trees over the AUC>=.75 threshold
Epcal_mean_acceptable_ntrees <- mean(Epcal_output_tab$ntrees[Epcal_output_tab$AUC_cv>=0.75])
Epcal_mean_acceptable_ntrees
noTopo_mean_acceptable_ntrees <- mean(noTopo_output_tab$ntrees[noTopo_output_tab$AUC_cv>=0.75])
noTopo_mean_acceptable_ntrees
Total_mean_acceptable_ntrees <- mean(Epcal_mean_acceptable_ntrees,noTopo_mean_acceptable_ntrees)
Total_mean_acceptable_ntrees



#### Get sample models for response curves (tagged with the kind of network they are) ####
#gbm.mod (No tag) = noTopo model (as it was the last one to be run)

#NoTopo
test = gbm.plot(gbm_mod_noTopo,return.grid=TRUE)
str(test)
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)
#Epcal
gbm_mod_Epcal = gbm.step(data=Epcal_4Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','habAv','EgoSize','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

#plot
gbm.plot(gbm_mod_Epcal, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,3), write.title=FALSE, cex.axis = 2,  cex.lab=2)
#NoTopo
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)


### specific stats of sample models:
summary(gbm_mod_Epcal)
gbm_mod_Epcal$self.statistics$discrimination
#noTopo
gbm_mod_noTopo$self.statistics$discrimination

write.csv(Epcal_4Vt_woNA, file = "BinPredRun_Epcal_4Vt_woNA_b_c_unw.csv")



# #################################################################################################
# #### Get discrete and rescaled predictions ####
# #################################################################################################
# #### Add continuous predictions to dataframe w/ all the predictors & occurrence-state ####
# #Do new dataframe to keep the original unchanged (only add underscore between gender & sp.)
# Ep_cal_stattest <- Epcal_stattest
# 
# 
# Epcal_predict_df = data.frame(Epcal_predict_mat)
# #Calculate the mean of the 100 models for each row (patch)
# # For past round it was done as uniform_stattest$cont_pred <- Singlemodel_Hy_arb_predict
# Epcal_predict_df$mean_pred <- rowMeans(Epcal_predict_mat, na.rm = FALSE)
# 
# Epcal_predict_df$PatchID <- Ep_cal_stattest$PatchID
# 
# Epcal_predict_df <- merge(Ep_cal_stattest, Epcal_predict_df, by = "PatchID")
# names(Epcal_predict_df)
# 
# ## Set discrete prediction to patches
# #Set discretization (binariation) cutpoint for each patch
# cp = cutpointr(Epcal_predict_df, mean_pred, pres_abs, subgroup = NULL, method = minimize_metric, metric = roc01, na.rm = TRUE)
# 
# summary(cp)
# plot(cp)
# cp$optimal_cutpoint
# #Set all above threshold to 1, all below to 0 (discrete = disc)
# Epcal_predict_df$pred_pres_abs <- ifelse(Epcal_predict_df$mean_pred>=cp$optimal_cutpoint, 1, 0)
# 
# table(Epcal_predict_df$pred_pres_abs)
# 
# #Add field w/ BRT Prediction score multiplied by 1000 (Mean_occurrence_index)
# #rescaled continuous prediction, to correspond to scale of HSI
# Epcal_predict_df$Pred_1000 <- Epcal_predict_df$mean_pred*1000
# 
# 
# #Export df to csv without the columns used only for calculation of binary prediction
# Ep_cal_stattest$mean_pred <- Epcal_predict_df$mean_pred
# Ep_cal_stattest$Pred_1000 <-Epcal_predict_df$Pred_1000
# Ep_cal_stattest$bin_OLpred <- Epcal_predict_df$pred_pres_abs

### For ArcGIS importing (compatibility of fields)
Ep_cal_stattest$Value <- Ep_cal_stattest$PatchID
Ep_cal_stattest$V_s <- Ep_cal_stattest$V_Ep_cal
head(Ep_cal_stattest)

Ep_cal_stattest <- moveMe(Ep_cal_stattest, c("Value"), "after", "PatchID")
Ep_cal_stattest <- moveMe(Ep_cal_stattest, c("V_s"), "after", "V_Ep_cal")

# setwd("C:/Users/damiano/Documents/PhD/Additional_species_runs/Discrete_predictions")
# write.csv(Ep_cal_stattest, file = "Ep_cal_bin_prediction.csv")


#Hard-save evaluation measures
setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
write.csv(Epcal_AUC_cv_tab, file = "Epcal_DefaultDispDist_AUC_cv_tab.csv")
write.csv(Epcal_AUC_train_tab, file = "Epcal_DefaultDispDist_AUC_train_tab.csv")
write.csv(Epcal_var_imp_tab, file = "Epcal_DefaultDispDist_var_imp_tab.csv")
write.csv(Epcal_output_tab, file = "Epcal_DefaultDispDist_BRToutput_tab.csv")
write.csv(Epcal_predict_df, file = "Epcal_DefaultDispDist_predict_df.csv")

#Hard-save variable importance of each predictor
write.csv(Epcal_HSI_tab, file = "Epcal_DefaultDispDist_HSI_tab.csv")
write.csv(Epcal_habAv_tab, file = "Epcal_DefaultDispDist_habAv_tab.csv")
write.csv(Epcal_deg_tab, file = "Epcal_DefaultDispDist_deg_tab.csv")
write.csv(Epcal_EgoSize_tab, file = "Epcal_DefaultDispDist_EgoSize_tab.csv")
write.csv(Epcal_unw_b_c_tab, file = "Epcal_DefaultDispDist_unw_b_c_tab.csv")
write.csv(Epcal_strength_tab, file = "Epcal_DefaultDispDist_strength_tab.csv")
write.csv(Epcal_Patch_Area_tab, file = "Epcal_DefaultDispDist_Patch_Area_tab.csv")

#Hard-save model to be  able to use it in other workspaces
saveRDS(gbm_mod_Epcal, file = "gbm_mod_Epcal_DefaultDispDist.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# readRDS(file, refhook = NULL)
# infoRDS(file)
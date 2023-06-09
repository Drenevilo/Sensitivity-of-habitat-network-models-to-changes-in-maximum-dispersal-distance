#### Script for: ##################################################################################################
#### Definition of Absences #### #### Calculation of predictors #### Boosted Regression Trees (BRT's) Modelling ###
#### For Bombina variegata #### #################################################################################
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

Macrotable_Bo_var <- read.csv("Macrotable_Bovar.csv")
head(Macrotable_Bo_var)
names(Macrotable_Bo_var)

#Drop OID  
Macrotable_Bo_var$OID <- NULL

#Reorder columns, to be consistent with species alphabetical order
Macrotable_Bo_var <- Macrotable_Bo_var[,order(colnames(Macrotable_Bo_var))] 
head(Macrotable_Bo_var)


#########Get binary values from the counts ####
#ifelse(df>0, 1, 0 ) 
#converts it to a matrix
Pres_abs_allpatches_Bo_var <- ifelse(Macrotable_Bo_var >0, 1, 0 )
#recover data.frame type
Pres_abs_allpatches_Bo_var <- as.data.frame(Pres_abs_allpatches_Bo_var, row.names = NULL, optional = FALSE)
names(Pres_abs_allpatches_Bo_var)
head(Pres_abs_allpatches_Bo_var)
#Drop binarized value & count
Pres_abs_allpatches_Bo_var$Value <- NULL
Pres_abs_allpatches_Bo_var$Count <- NULL

#Recovering the identifier columns
Pres_abs_allpatches_Bo_var$Value <- Macrotable_Bo_var$Value
Pres_abs_allpatches_Bo_var$Count <- Macrotable_Bo_var$Count
names(Pres_abs_allpatches_Bo_var)

#Get fields VALUE & COUNT (patch area in hectares) to the beginning
Pres_abs_allpatches_Bo_var <- moveMe(Pres_abs_allpatches_Bo_var, c("Value", "Count"), "first")
#head(Pres_abs_allpatches_Bo_var)
names(Pres_abs_allpatches_Bo_var)


#Get V_sp - Observations of Gender_species
#(sum(1 or 0 each year)
Pres_abs_allpatches_Bo_var$V_Bo_var <- rowSums(Pres_abs_allpatches_Bo_var[grep("Bo_var", names(Pres_abs_allpatches_Bo_var))])
names(Pres_abs_allpatches_Bo_var)


#### Calculate observations of any amphibian per year ####

for(year in c("06","07","08","09","10","11","12","13","14","15")){
  Pres_abs_allpatches_Bo_var[,paste("V_", year, sep = "")]<- rowSums(Pres_abs_allpatches_Bo_var[grep(year, names(Pres_abs_allpatches_Bo_var))])
}

Pres_abs_allpatches_Bo_var$V_06 <- ifelse(Pres_abs_allpatches_Bo_var$V_06>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_07 <- ifelse(Pres_abs_allpatches_Bo_var$V_07>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_08 <- ifelse(Pres_abs_allpatches_Bo_var$V_08>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_09 <- ifelse(Pres_abs_allpatches_Bo_var$V_09>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_10 <- ifelse(Pres_abs_allpatches_Bo_var$V_10>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_11 <- ifelse(Pres_abs_allpatches_Bo_var$V_11>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_12 <- ifelse(Pres_abs_allpatches_Bo_var$V_12>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_13 <- ifelse(Pres_abs_allpatches_Bo_var$V_13>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_14 <- ifelse(Pres_abs_allpatches_Bo_var$V_14>0, 1, 0)
Pres_abs_allpatches_Bo_var$V_15 <- ifelse(Pres_abs_allpatches_Bo_var$V_15>0, 1, 0)


### Get V_t - Total observations of amphibians
#By summing the vectors of spp per year
Pres_abs_allpatches_Bo_var$V_t <- rowSums(Pres_abs_allpatches_Bo_var[,c("V_06","V_07","V_08","V_09","V_10","V_11","V_12","V_13","V_14","V_15")])


################

#Plot Vt vs Vs (Vh)
boxplot(Pres_abs_allpatches_Bo_var$V_Bo_var ~ Pres_abs_allpatches_Bo_var$V_t, xlab= "Vt", ylab = "Vs")
#plot bars indicating median 

write.csv(Pres_abs_allpatches_Bo_var, file = "Pres_abs_Vs_Vt_Bo_var.csv") 

#Total number of patches
length(Pres_abs_allpatches_Bo_var$V_t)


################ Occurence-state (Pres_abs) decider Plots #########################

#Pres_abs_allpatches2 is the one modified to make plots, excludes all the NA 
#Pres_abs_allpatches (unchanged, with all the 0's) will be used to incorporate the graph properties
#Get only the records that are visits (Vt >= 1)
Pres_abs_allpatches_Bo_var2 = Pres_abs_allpatches_Bo_var[Pres_abs_allpatches_Bo_var$V_t != 0,] 

#Bar graph with proportions, averages and text --> The main one to decide thresholds!!!

#Data for the plot #Taking mean Vs
tab = ddply(Pres_abs_allpatches_Bo_var2, .(V_t), summarize,  mean=mean(V_Bo_var), countH = length(which(V_Bo_var>0)), countT = length(V_Bo_var), lab = paste(length(which(V_Bo_var>0)),"/",length(V_Bo_var),sep = " "), perc = length(which(V_Bo_var>0))/length(V_Bo_var)*100)

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
mean(Pres_abs_allpatches_Bo_var$V_Bo_var[Pres_abs_allpatches_Bo_var$V_t==4])

#Number of absences at threshold value:
sum(Pres_abs_allpatches_Bo_var$V_t>3 & Pres_abs_allpatches_Bo_var$V_Bo_var<1) 

which(Pres_abs_allpatches_Bo_var$V_t >= 4 & Pres_abs_allpatches_Bo_var$V_Bo_var ==0)


#### rename(data, d = b) ####
Pres_abs_allpatches_Bo_var <- setnames(Pres_abs_allpatches_Bo_var, "Value", "PatchID")
Pres_abs_allpatches_Bo_var <- setnames(Pres_abs_allpatches_Bo_var, "Count", "Patch_Area")


#### Add HSI ####
HSI = read.csv("C:/Users/damiano/Documents/PhD/Additional_species_runs/Network_setup/HSI_habPatchCode_Bovar1.csv")
Pres_abs_allpatches_Bo_var$HSI <- HSI$"MEAN"



######## Importing network ############

setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/Networks')

#Import table of patches w/network properties
Bovar_graph <- read_graph("Bovar.graphml", format = "graphml") 

#Look at vertices and edges of a graph
E(Bovar_graph)
V(Bovar_graph)


#### get metrics as attributes and into a data.frame ####
Bovar_graph_metrics <- data.frame(
  PatchID = V(Bovar_graph)$name,
  deg=degree(Bovar_graph),
  strength = strength(Bovar_graph, vids = V(Bovar_graph), mode = "all", loops = FALSE, weights = E(Bovar_graph)$weight),
  EgoSize = ego_size(Bovar_graph, order = 3, nodes = V(Bovar_graph))
  
)
### EgoSize = Third-order neighborhood


#### join the tables by attribute
#inner join= join by common variable names
patches_Bovar_topo_attributes <- merge(Pres_abs_allpatches_Bo_var, Bovar_graph_metrics, by = "PatchID")


#### Add Habitat Availability predictor ####
Bovar_habAv <- scan("Bovar_habAv.txt")
patches_Bovar_topo_attributes$habAv <- Bovar_habAv
head(patches_Bovar_topo_attributes$habAv)


#### Measure betweenness centrality without weight included ####
#Re-import graphs, change name of 'weight' attribute
#setwd('/Networks')
Bovar_graph2 <- read_graph("Bovar.graphml", format = "graphml")

#rename weight attribute
#Actually copying into a new attr. and deleting original attribute
E(Bovar_graph2)$Alter_wght <- E(Bovar_graph2)$weight
Bovar_graph2 <- remove.edge.attribute(Bovar_graph2, "weight")

#Calculate explicitely unweighted b_c
#Set dataframe w/only PatchID & unw_b_c
Bovar_unweighted_b_c <- data.frame(
  PatchID = V(Bovar_graph2)$name,
  b_c=betweenness(Bovar_graph2, weights = NULL))

#rename unw_b_c
Bovar_unweighted_b_c <- rename(Bovar_unweighted_b_c,  "unw_b_c" = "b_c")

#join the tables by attribute
patches_Bovar_topo_attributes <- merge(x = patches_Bovar_topo_attributes, y = Bovar_unweighted_b_c, by = "PatchID", all.x = TRUE)
names(patches_Bovar_topo_attributes)
head(patches_Bovar_topo_attributes)



#### Define occurrence_state (pres_abs) ####
### Classify the dataframe in subsets that classify all the patches in presences, absences or questionmarks/NoIdea (NA) 

#subset of questionmarks, #All the records, except the ones that comply with the following commands
patches_Bovar_topo_attributes[,"pres_abs"] = NA 
#Subset of likely absences (With threshold at 5 (defined by pres_abs decider plot))
patches_Bovar_topo_attributes[which(patches_Bovar_topo_attributes$V_t >= 4 & patches_Bovar_topo_attributes$V_Bo_var ==0),"pres_abs"] = 0
#Subset of confirmed presences
patches_Bovar_topo_attributes[which(patches_Bovar_topo_attributes$V_Bo_var > 0),"pres_abs"] = 1

#### Display the 'occurrence-state' (0/1/NA) of every patch as defined above (pres_abs) ####
patches_Bovar_topo_attributes$pres_abs

#Get total amount of values for presence and for absence
table(patches_Bovar_topo_attributes$pres_abs)


#### Do a subset of the columns that only has the predictors and response ####
Bovar_stattest <- subset(patches_Bovar_topo_attributes, select=c("PatchID", "Patch_Area", "V_Bo_var", "V_t", "deg", "unw_b_c", "strength", "EgoSize", "HSI", "habAv", "pres_abs"))
names(Bovar_stattest)
head(Bovar_stattest)

#change directory to the one with BRT things 
setwd('C:/Users/damiano/Documents/PhD/Additional_species_runs/BRTs')
#write file
#t4 -> pres_abs threshold=4
write.csv(Bovar_stattest, file = "Bovar_stattest_t4.csv")


##### Do the 'only pres_abs (no-NoData) df's #########################
#Take out NA's, set threshold of absences (at V_t = 4 for Bo_var), change name
Bovar_4Vt_woNA = Bovar_stattest
Bovar_4Vt_woNA[Bovar_4Vt_woNA$V_Bo_var==0 & Bovar_4Vt_woNA$V_t >= 4, 'pres_abs'] = 0
Bovar_4Vt_woNA = Bovar_4Vt_woNA[!is.na(Bovar_4Vt_woNA$pres_abs),]

#check values of unw_b_c
mean(Bovar_4Vt_woNA$unw_b_c)



#### Perform BRT's ############################################################
#Hyperparameters for all models
lr = 0.001
tc = 5
bf = 0.75

n_repeats = 100

##### Bovar (Full model with all the predictors) ###########################################

# Create the output table that contains all the values of each of the runs
Bovar_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar
# Create vectors to record performance measures
Bovar_AUC_cv_vec = vector() #Cross-validated AUC
Bovar_AUC_vec = vector() #Training AUC

Bovar_nt_vec = vector() #Number of trees

#Create vectors with all the values of a certain predictor along all the runs
Bovar_HSI_vec = vector() 
Bovar_EgoSize_vec = vector() 
Bovar_strength_vec  = vector()
Bovar_role_vec  = vector()
Bovar_deg_vec  = vector()
Bovar_habAv_vec  = vector()
Bovar_unw_b_c_vec  = vector()
Bovar_Patch_Area_vec  = vector()

#Discrete Prediction of occurrence state for all patches
Bovar_predict_mat = matrix(nrow = length(Bovar_stattest$PatchID), ncol = n_repeats, byrow = FALSE)


### Loop that goes exactly for 100 iterations, to get distributions 
for(i in c(1:n_repeats)){
  
  #Perform gbm step to set number of trees, no cross-validation.
  gbm_mod_Bovar = gbm.step(data=Bovar_4Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','EgoSize','HSI', 'habAv'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
  #data frame of variable importance, to fill the vectors of the model var. importance scores
  var_imp = data.frame(var = summary(gbm_mod_Bovar)$var, imp = summary(gbm_mod_Bovar)$rel.inf)
  
  #Write var_imp results to the vectors
  Bovar_HSI_vec = append(Bovar_HSI_vec, var_imp$imp[var_imp$var=='HSI'])
  Bovar_EgoSize_vec = append(Bovar_EgoSize_vec, var_imp$imp[var_imp$var=='EgoSize'])
  Bovar_strength_vec  = append(Bovar_strength_vec, var_imp$imp[var_imp$var=='strength'])
  Bovar_deg_vec  = append(Bovar_deg_vec, var_imp$imp[var_imp$var=='deg'])
  Bovar_habAv_vec  = append(Bovar_habAv_vec, var_imp$imp[var_imp$var=='habAv'])
  Bovar_unw_b_c_vec  = append(Bovar_unw_b_c_vec, var_imp$imp[var_imp$var=='unw_b_c'])
  Bovar_Patch_Area_vec  = append(Bovar_Patch_Area_vec, var_imp$imp[var_imp$var=='Patch_Area'])
  
  #Write the AUC & CV-AUC of this model to a vector.
  Bovar_AUC_cv_vec = append(Bovar_AUC_cv_vec, gbm_mod_Bovar$cv.statistics$discrimination.mean)
  Bovar_AUC_vec = append(Bovar_AUC_vec, gbm_mod_Bovar$self.statistics$discrimination)
  
  #Write the number of trees to a vector
  nt = gbm_mod_Bovar$n.trees
  Bovar_nt_vec = append(Bovar_nt_vec, nt)
  print(nt)
  
  #write the continuous prediction over all the patches
  Bovar_predict_mat[,i] = predict(gbm_mod_Bovar, Bovar_stattest, gbm_mod_Bovar$n.trees, type = "response", single.tree = FALSE)
  
  print(paste("Finished:",i,"/",n_repeats,sep = ""))
}


Bovar_output_tab[,"AUC_train"] = Bovar_AUC_vec
Bovar_output_tab[,"AUC_cv"] = Bovar_AUC_cv_vec
Bovar_output_tab[,"ntrees"] = Bovar_nt_vec


#Var. importance columns
Bovar_output_tab[,"HSI_imp"] = Bovar_HSI_vec
Bovar_output_tab[,"EgoSize_imp"] = Bovar_EgoSize_vec
Bovar_output_tab[,"strength_imp"] = Bovar_strength_vec
Bovar_output_tab[,"deg_imp"] = Bovar_deg_vec
Bovar_output_tab[,"habAv_imp"] = Bovar_habAv_vec
Bovar_output_tab[,"unw_b_c_imp"] = Bovar_unw_b_c_vec
Bovar_output_tab[,"Patch_Area_imp"] = Bovar_Patch_Area_vec


# Make a dataframe to make plots for comparison
Bovar_HSI_tab = data.frame(imp = Bovar_HSI_vec)
Bovar_EgoSize_tab = data.frame(imp = Bovar_EgoSize_vec)
Bovar_strength_tab = data.frame(imp = Bovar_strength_vec)
Bovar_deg_tab = data.frame(imp = Bovar_deg_vec)
Bovar_habAv_tab = data.frame(imp = Bovar_habAv_vec)
Bovar_unw_b_c_tab = data.frame(imp = Bovar_unw_b_c_vec)
Bovar_Patch_Area_tab = data.frame(imp = Bovar_Patch_Area_vec)

Bovar_HSI_tab$variable = "HSI"
Bovar_EgoSize_tab$variable = "3rd. ord. neigh."
Bovar_strength_tab$variable = "Strength"
Bovar_deg_tab$variable = "Degree"
Bovar_habAv_tab$variable = "Hab. Av."
Bovar_unw_b_c_tab$variable = "B.C."
Bovar_Patch_Area_tab$variable = "Patch Area"

#Reserve also measures in df to compare performance between models
Bovar_AUC_cv_tab = data.frame(value = Bovar_AUC_cv_vec)
Bovar_AUC_train_tab = data.frame(value = Bovar_AUC_vec)
#Make label of model for plot
Bovar_AUC_cv_tab$model = "Bovar"
Bovar_AUC_train_tab$model = "Bovar"

#combine pred. vars. into new data frame 
Bovar_var_imp_tab = rbind(Bovar_HSI_tab,Bovar_EgoSize_tab,Bovar_strength_tab,Bovar_habAv_tab,Bovar_deg_tab,Bovar_unw_b_c_tab,Bovar_Patch_Area_tab)

ggplot(Bovar_var_imp_tab, aes(imp, fill = variable)) + geom_density(alpha = 0.2)
boxplot(Bovar_var_imp_tab$imp~Bovar_var_imp_tab$variable, ylab= "Var. Importance", main = "Bovar")


#Get mean var. importance of all of the vars. 
mean(Bovar_HSI_vec)
mean(Bovar_EgoSize_vec)
mean(Bovar_strength_vec)
mean(Bovar_deg_vec)
mean(Bovar_habAv_vec)
mean(Bovar_unw_b_c_vec)
mean(Bovar_Patch_Area_vec)

##Check distr. of measures of prediction accuracy 
summary(Bovar_output_tab$AUC_cv)
summary(Bovar_output_tab$AUC_train)



############ NoTopo (Model with only non-topological predictors included) #################################

# Create the output table that contains all the values of each of the runs
noTopo_output_tab = data.frame(run_nr = c(1:n_repeats)) #Bovar
# Create vectors to record performance measures
# noTopo_ce_MK_vec = vector() #classification error w MaxKappa
# noTopo_MK_th_vec = vector() #MaxKappa binarization threshold 
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
  gbm_mod_noTopo = gbm.step(data=Bovar_4Vt_woNA, gbm.x = c('Patch_Area','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE) 
  
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
  # noTopo_predict_vec = append(predict(gbm_mod_noTopo, Bovar_stattest, gbm_mod_noTopo$n.trees, type = "response", single.tree = FALSE))
  
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
boxplot(noTopo_var_imp_tab$imp~noTopo_var_imp_tab$variable, ylab= "Var. Importance", main = "noTopo_Bovar")

#Get mean var. importance of all of the vars. 
mean(noTopo_HSI_vec)
mean(noTopo_Patch_Area_vec)


##Check distr. of measures of prediction accuracy 
summary(noTopo_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_cv)



##### Compare between runs ###########################

#cv AUC
#mean
AUC_cv_tab = rbind(Bovar_AUC_cv_tab, noTopo_AUC_cv_tab)

#Change order of factors to display noTopo at the edge
AUC_cv_tab$model<- as.factor(AUC_cv_tab$model)
levels(AUC_cv_tab$model)
AUC_cv_tab$model<-factor(AUC_cv_tab$model, levels=c("Bovar", "noTopo"))
#print out
levels(AUC_cv_tab$model)

#Plot
ggplot(AUC_cv_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_cv_tab$value~AUC_cv_tab$model, ylab= "Cross-validated AUC")


#training AUC
AUC_train_tab = rbind(Bovar_AUC_train_tab, noTopo_AUC_train_tab)

ggplot(AUC_train_tab, aes(value, fill = model)) + geom_density(alpha = 0.2)
boxplot(AUC_train_tab$value~AUC_train_tab$model, ylab= "Training AUC")


summary(Bovar_output_tab$AUC_cv)
summary(noTopo_output_tab$AUC_cv)

summary(Bovar_output_tab$AUC_train)
summary(noTopo_output_tab$AUC_train)

summary(Bovar_output_tab$ntrees)
summary(noTopo_output_tab$ntrees)


###Get avg_#trees over the AUC>=.75 threshold
Bovar_mean_acceptable_ntrees <- mean(Bovar_output_tab$ntrees[Bovar_output_tab$AUC_cv>=0.75])
Bovar_mean_acceptable_ntrees
noTopo_mean_acceptable_ntrees <- mean(noTopo_output_tab$ntrees[noTopo_output_tab$AUC_cv>=0.75])
noTopo_mean_acceptable_ntrees
Total_mean_acceptable_ntrees <- mean(Bovar_mean_acceptable_ntrees,noTopo_mean_acceptable_ntrees)
Total_mean_acceptable_ntrees



#### Get sample models for response curves (tagged with the kind of network they are) ####
#gbm.mod (No tag) = noTopo model (as it was the last one to be run)

#NoTopo
test = gbm.plot(gbm_mod_noTopo,return.grid=TRUE)
str(test)
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)
#Bovar
gbm_mod_Bovar = gbm.step(data=Bovar_4Vt_woNA, gbm.x = c('Patch_Area','deg','unw_b_c','strength','habAv','EgoSize','HSI'), gbm.y = 'pres_abs', family = "bernoulli", learning.rate = lr, tree.complexity = tc, bag.fraction = bf, silent = TRUE)

#plot
gbm.plot(gbm_mod_Bovar, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,3), write.title=FALSE, cex.axis = 2,  cex.lab=2)
#NoTopo
gbm.plot(gbm_mod_noTopo, smooth = FALSE, common.scale = FALSE, plot.layout=c(1,2), write.title=FALSE, cex.axis = 1.5,  cex.lab=1.5)

### specific stats of sample models:
summary(gbm_mod_Bovar)
gbm_mod_Bovar$self.statistics$discrimination

#noTopo
gbm_mod_noTopo$self.statistics$discrimination

write.csv(Bovar_4Vt_woNA, file = "BinPredRun_Bovar_4Vt_woNA_b_c_unw.csv")



# #################################################################################################
# #### Get discrete and rescaled predictions ####
# #################################################################################################
# #### Add continuous predictions to dataframe w/ all the predictors & occurrence-state ####
# #Do new dataframe to keep the original unchanged (only add underscore between gender & sp.)
# Bo_var_stattest <- Bovar_stattest
# 
# 
# Bovar_predict_df = data.frame(Bovar_predict_mat)
# #Calculate the mean of the 100 models for each row (patch)
# # For past round it was done as uniform_stattest$cont_pred <- Singlemodel_Hy_arb_predict
# Bovar_predict_df$mean_pred <- rowMeans(Bovar_predict_mat, na.rm = FALSE)
# 
# Bovar_predict_df$PatchID <- Bo_var_stattest$PatchID
# 
# Bovar_predict_df <- merge(Bo_var_stattest, Bovar_predict_df, by = "PatchID")
# names(Bovar_predict_df)
# 
# ## Set discrete prediction to patches
# #Set discretization (binariation) cutpoint for each patch
# cp = cutpointr(Bovar_predict_df, mean_pred, pres_abs, subgroup = NULL, method = minimize_metric, metric = roc01, na.rm = TRUE)
# 
# summary(cp)
# plot(cp)
# cp$optimal_cutpoint
# #Set all above threshold to 1, all below to 0 (discrete = disc)
# Bovar_predict_df$pred_pres_abs <- ifelse(Bovar_predict_df$mean_pred>=cp$optimal_cutpoint, 1, 0)
# 
# table(Bovar_predict_df$pred_pres_abs)
# 
# #Add field w/ BRT Prediction score multiplied by 1000 (Mean_occurrence_index)
# #rescaled continuous prediction, to correspond to scale of HSI
# Bovar_predict_df$Pred_1000 <- Bovar_predict_df$mean_pred*1000
# 
# 
# #Export df to csv without the columns used only for calculation of binary prediction
# Bo_var_stattest$mean_pred <- Bovar_predict_df$mean_pred
# Bo_var_stattest$Pred_1000 <-Bovar_predict_df$Pred_1000
# Bo_var_stattest$bin_OLpred <- Bovar_predict_df$pred_pres_abs


### For ArcGIS importing (compatibility of fields)
Bo_var_stattest$Value <- Bo_var_stattest$PatchID
Bo_var_stattest$V_s <- Bo_var_stattest$V_Bo_var
head(Bo_var_stattest)

Bo_var_stattest <- moveMe(Bo_var_stattest, c("Value"), "after", "PatchID")
Bo_var_stattest <- moveMe(Bo_var_stattest, c("V_s"), "after", "V_Bo_var")

# setwd("C:/Users/damiano/Documents/PhD/Additional_species_runs/Discrete_predictions")
# write.csv(Bo_var_stattest, file = "Bo_var_bin_prediction.csv")


#Hard-save evaluation measures
setwd("C:/Users/damiano/Documents/PhD/Sensitivity_DispDist/BRTs")
write.csv(Bovar_AUC_cv_tab, file = "Bovar_DefaultDispDist_AUC_cv_tab.csv")
write.csv(Bovar_AUC_train_tab, file = "Bovar_DefaultDispDist_AUC_train_tab.csv")
write.csv(Bovar_var_imp_tab, file = "Bovar_DefaultDispDist_var_imp_tab.csv")
write.csv(Bovar_output_tab, file = "Bovar_DefaultDispDist_BRToutput_tab.csv")
write.csv(Bovar_predict_df, file = "Bovar_DefaultDispDist_predict_df.csv")

#Hard-save variable importance of each predictor
write.csv(Bovar_HSI_tab, file = "Bovar_DefaultDispDist_HSI_tab.csv")
write.csv(Bovar_habAv_tab, file = "Bovar_DefaultDispDist_habAv_tab.csv")
write.csv(Bovar_deg_tab, file = "Bovar_DefaultDispDist_deg_tab.csv")
write.csv(Bovar_EgoSize_tab, file = "Bovar_DefaultDispDist_EgoSize_tab.csv")
write.csv(Bovar_unw_b_c_tab, file = "Bovar_DefaultDispDist_unw_b_c_tab.csv")
write.csv(Bovar_strength_tab, file = "Bovar_DefaultDispDist_strength_tab.csv")
write.csv(Bovar_Patch_Area_tab, file = "Bovar_DefaultDispDist_Patch_Area_tab.csv")

#Hard-save model to be  able to use it in other workspaces
saveRDS(gbm_mod_Bovar, file = "gbm_mod_Bovar_DefaultDispDist.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# readRDS(file, refhook = NULL)
# infoRDS(file)

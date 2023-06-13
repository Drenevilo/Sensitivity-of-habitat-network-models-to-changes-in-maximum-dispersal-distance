# Sensitivity-of-habitat-network-models-to-changes-in-maximum-dispersal-distance
Code used for developing the research presented in the article "Sensitivity of habitat network models to changes in maximum dispersal distance", by Damian O. Ortiz‐Rodríguez, Antoine Guisan and Maarten J. van Strien

Processes ordered by sequence, coded by number at the start of the file name:

1. Definition of study area and species presence records, pre-processing of predictor variables for HSM, defintion of environmental mask
2. Habitat Suitability Models for all species and a generic one
3.  Mask application, habitat patch definition, calculation of mean HSI per habitat patch
4.  Cost surface definition
5.  Edge definition (Generation of networks) & calculation of Habitat availability predictor with Species-specific maximum dispersal distance
6.  Edge definition (Generation of networks) & calculation of Habitat availability predictor, changing the d0 parameter to implement different maximum dispersal distances
7.  Definition of Absences, Calculation of topological predictors, Boosted Regression Trees (BRT's) Modelling - For species-specific maximum dispersal distance
8.  Definition of Absences, Calculation of topological predictors, Boosted Regression Trees (BRT's) Modelling - For all the alternative max. disp. dists. +  Comparison of the performance among all the network models with different maximum dispersal distances
8a. Table used in 8 (Line 8453), which was generated with the same data, but outside of R   

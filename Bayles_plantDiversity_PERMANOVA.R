# This analyzes differences in plant community composition directly using PERMANOVA 
# using Morisita's dissimilarity index
# (instead of Shannon diversity index) for changes with different treatment combinations
# (+/- AMF, +/- biochar, location in each plot). This returns the results reported in the 
# manuscript.

if(!(require("vegan"))){
    install.packages("vegan")
} 

library(vegan)

Bayles_plantDiversity_2013 <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_aggregatedPlantDiversityData_2013.tsv", header = TRUE, sep = "\t")

Bayles_plantDiversity_2014 <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_aggregatedPlantDiversityData_2014.tsv", header = TRUE, sep = "\t")

Bayles_plantDiversity_2016 <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_aggregatedPlantDiversityData_2016.tsv", header = TRUE, sep = "\t")

# Now load the plot information (AMF/Biochar/Blocks)
Bayles_plotTreatments <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_plotTreatments.csv", header = TRUE, sep = ",")

plotTreatments_2013 <- data.frame("Plot" = as.factor(rep(seq(1,36,1), each = 2)), "plotLoc" = as.factor(rep(c("C","R"),18)), "AMFTreatment" = as.factor(rep(Bayles_plotTreatments$AMF, each = 2)),
                                  "BiocharTreatment" = as.factor(rep(Bayles_plotTreatments$Biochar, each = 2)), "Block" = as.factor(rep(Bayles_plotTreatments$Block, each = 2)))

plotTreatments_2014 <- data.frame("Plot" = as.factor(rep(seq(1,36,1), each = 2)), "plotLoc" = as.factor(rep(c("C","R"),18)), "AMFTreatment" = as.factor(rep(Bayles_plotTreatments$AMF, each = 2)),
                                  "BiocharTreatment" = as.factor(rep(Bayles_plotTreatments$Biochar, each = 2)), "Block" = as.factor(rep(Bayles_plotTreatments$Block, each = 2)))

plotTreatments_2016 <- data.frame("Plot" = as.factor(Bayles_plantDiversity_2016$plotNum), "plotLoc" = as.factor(rep("C",30)), "AMFTreatment" = as.factor(Bayles_plotTreatments$AMF[1:30]),
                                  "BiocharTreatment" = as.factor(Bayles_plotTreatments$Biochar[1:30]), "Block" = as.factor(Bayles_plotTreatments$Block[1:30]))

# ==========================

# Across all 3 years
uniquePlantNames <- sort(unique(c(colnames(Bayles_plantDiversity_2013[-1]), colnames(Bayles_plantDiversity_2014[-1]), colnames(Bayles_plantDiversity_2016[-1]))))


# Native (Nat) designations for each of the unique plant names. Any with ambiguous species
# are assigned to non-native. If 'both' in USDA plants, giving it 'Nat' designation.

nativeHolder <- c("Nat","No","No","Nat", "Nat", "Nat", "Nat", "Nat", "Nat", "Nat", "Nat", "No",
                  "No", "Nat", "No", "No", "No", "No", "No", "Nat", "No", "No", "Nat",
                  "Nat", "Nat", "No", "Nat", "No", "Nat", "Nat", "No", "Nat", "Nat", "No", "No",
                  "No", "No", "Nat", "Nat", "No", "Nat", "Nat", "Nat", "Nat", "No", "No", "Nat",
                  "Nat", "No", "No", "Nat", "No", "No", "Nat", "No", "Nat", "Nat", "Nat", "Nat",
                  "No", "No", "Nat", "No", "No", "No", "Nat", "Nat", "No", "No", "No", "Nat",
                  "No", "No", "Nat", "Nat", "No", "Nat", "Nat", "No", "Nat", "Nat", "No", "No",
                  "No", "No", "No", "No", "Nat", "Nat", "No", "No", "No", "No", "No", "No",
                  "Nat", "No", "Nat", "No", "No", "No")

plantNameStatus <- cbind(uniquePlantNames, nativeHolder)

# Prep data for diversity calc. Will have sep calc for native and overall (all plants).
# The diversity calc with vegan can't have a plot column, so remove that here (the rows are in order by plot from 1-36 anyway)

BaylesPlantDiv_aggregated_2013_overall <- Bayles_plantDiversity_2013[,-1]
BaylesPlantDiv_aggregated_2014_overall <- Bayles_plantDiversity_2014[,-1]
BaylesPlantDiv_aggregated_2016_overall <- Bayles_plantDiversity_2016[,-1]

# For each of the plants represented for each of the years, this returns the native/non status for each by matching
# it by name with the plant name/status in the plantNameStatus
plantNameStatus_2013 <- sapply(X = colnames(BaylesPlantDiv_aggregated_2013_overall), FUN = function(x){plantNameStatus[which(plantNameStatus[,1] == x),2]})
BaylesPlantDiv_aggregated_2013_nativeONLY <- BaylesPlantDiv_aggregated_2013_overall[,which(plantNameStatus_2013 == "Nat")]

plantNameStatus_2014 <- sapply(X = colnames(BaylesPlantDiv_aggregated_2014_overall), FUN = function(x){plantNameStatus[which(plantNameStatus[,1] == x),2]})
BaylesPlantDiv_aggregated_2014_nativeONLY <- BaylesPlantDiv_aggregated_2014_overall[,which(plantNameStatus_2014 == "Nat")]

plantNameStatus_2016 <- sapply(X = colnames(BaylesPlantDiv_aggregated_2016_overall), FUN = function(x){plantNameStatus[which(plantNameStatus[,1] == x),2]})
BaylesPlantDiv_aggregated_2016_nativeONLY <- BaylesPlantDiv_aggregated_2016_overall[,which(plantNameStatus_2016 == "Nat")]

# 2013
# ======================
# Overall. Full plot Plot loc p = 0.07778
set.seed(729)
bayles_plant_permanova_2013_all <- adonis2(formula = BaylesPlantDiv_aggregated_2013_overall ~ plotTreatments_2013$plotLoc * plotTreatments_2013$AMFTreatment * plotTreatments_2013$BiocharTreatment, strata = plotTreatments_2013$Block,
                                            permutations = 5000, method = "morisita", binary = FALSE)
print(bayles_plant_permanova_2013_all)
# =====================

# 2014
# ======================
# Overall. Nothing sign. Full plot Plot loc p = 0.3555
set.seed(729)
bayles_plant_permanova_2014_all <- adonis2(formula = BaylesPlantDiv_aggregated_2014_overall ~ plotTreatments_2014$plotLoc * plotTreatments_2014$AMFTreatment * plotTreatments_2014$BiocharTreatment, strata = plotTreatments_2014$Block,
                                           permutations = 5000, method = "morisita", binary = FALSE)
print(bayles_plant_permanova_2014_all)

# Native ONLY; plot location is sign. p = 0.03579
set.seed(729)
bayles_plant_permanova_2014_native <- adonis2(formula = BaylesPlantDiv_aggregated_2014_nativeONLY ~ plotTreatments_2014$plotLoc * plotTreatments_2014$AMFTreatment * plotTreatments_2014$BiocharTreatment, strata = plotTreatments_2014$Block,
                                              permutations = 5000, method = "morisita", binary = FALSE)
print(bayles_plant_permanova_2014_native)
# =====================
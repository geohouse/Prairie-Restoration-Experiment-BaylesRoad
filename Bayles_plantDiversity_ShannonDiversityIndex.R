# This uses linear mixed effects models to test for differences in the Shannon diversity
# index of plant community composition between different treatment combinations. This 
# produces the results cited in the manuscript and plots the graph for Fig. S2.

# overall richness: 2013 - 44, 2014 - 67, 2016 - 66
# native richness: 2013 - 27, 2014 - 32, 2016 - 31

if(!(require("vegan"))){
    install.packages("vegan")
} 

if(!(require("lme4"))){
    install.packages("lme4")
} 

if(!(require("lmerTest"))){
    install.packages("lmerTest")
} 

if(!(require("ggplot2"))){
    install.packages("ggplot2")
} 

library(vegan)
library(lme4)
library(lmerTest)
library(ggplot2)

# SAS - like contrasts (dropping last factor)
options(contrasts=c(factor = "contr.SAS", ordered = "contr.poly"))

# Import the plant counts data
Bayles_plantDiversity_2013 <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_aggregatedPlantDiversityData_2013.tsv", header = TRUE, sep = "\t")

Bayles_plantDiversity_2014 <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_aggregatedPlantDiversityData_2014.tsv", header = TRUE, sep = "\t")

Bayles_plantDiversity_2016 <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_aggregatedPlantDiversityData_2016.tsv", header = TRUE, sep = "\t")

# Now load the plot information (AMF/Biochar/Blocks)
Bayles_plotTreatments <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_plotTreatments.csv", header = TRUE, sep = ",")

# Across all 3 years
uniquePlantNames <- sort(unique(c(colnames(Bayles_plantDiversity_2013[-1]), colnames(Bayles_plantDiversity_2014[-1]), colnames(Bayles_plantDiversity_2016[-1]))))

# Native (Nat) designations for each of the unique plant names. Any with ambiguous species
# are assigned to non-native. If 'both' in USDA plants, giving it 'Nat' designation.

nativeHolder <- c("Nat","No","No","Nat", "Nat", "Nat", "Nat", "Nat", "Nat", "Nat", "Nat", "No",
                  "No", "Nat", "No", "No", "No", "No", "No", "Nat", "No", "No", "Nat",
                  "Nat", "Nat", "No", "Nat", "No", "Nat", "Nat", "No", "Nat", "Nat", "Nat", "No",
                  "No", "No", "Nat", "Nat", "No", "Nat", "Nat", "Nat", "Nat", "No", "No", "Nat",
                  "Nat", "No", "No", "Nat", "No", "No", "Nat", "No", "Nat", "Nat", "Nat", "Nat",
                  "No", "No", "Nat", "No", "No", "No", "Nat", "Nat", "No", "No", "No", "Nat",
                  "No", "No", "Nat", "Nat", "No", "Nat", "Nat", "No", "Nat", "Nat", "No", "No",
                  "No", "No", "No", "No", "Nat", "Nat", "No", "No", "No", "No", "No", "No",
                  "Nat", "No", "Nat", "No", "No", "No")

plantNameStatus <- cbind(uniquePlantNames, nativeHolder)

# Prep data for diversity calc. Will have sep calc for native and overall (all plants).
# The diversity calc with vegan can't have a plot column, so remove that here (the rows are in order by plot from 1-36 anyway)

BaylesPlantDiv_aggregated_2013_forDiversity_overall <- Bayles_plantDiversity_2013[,-1]
BaylesPlantDiv_aggregated_2014_forDiversity_overall <- Bayles_plantDiversity_2014[,-1]
BaylesPlantDiv_aggregated_2016_forDiversity_overall <- Bayles_plantDiversity_2016[,-1]

# For each of the plants represented for each of the years, this returns the native/non status for each by matching
# it by name with the plant name/status in the plantNameStatus
plantNameStatus_2013 <- sapply(X = colnames(BaylesPlantDiv_aggregated_2013_forDiversity_overall), FUN = function(x){plantNameStatus[which(plantNameStatus[,1] == x),2]})
BaylesPlantDiv_aggregated_2013_forDiversity_nativeONLY <- BaylesPlantDiv_aggregated_2013_forDiversity_overall[,which(plantNameStatus_2013 == "Nat")]

plantNameStatus_2014 <- sapply(X = colnames(BaylesPlantDiv_aggregated_2014_forDiversity_overall), FUN = function(x){plantNameStatus[which(plantNameStatus[,1] == x),2]})
BaylesPlantDiv_aggregated_2014_forDiversity_nativeONLY <- BaylesPlantDiv_aggregated_2014_forDiversity_overall[,which(plantNameStatus_2014 == "Nat")]

plantNameStatus_2016 <- sapply(X = colnames(BaylesPlantDiv_aggregated_2016_forDiversity_overall), FUN = function(x){plantNameStatus[which(plantNameStatus[,1] == x),2]})
BaylesPlantDiv_aggregated_2016_forDiversity_nativeONLY <- BaylesPlantDiv_aggregated_2016_forDiversity_overall[,which(plantNameStatus_2016 == "Nat")]

# Now run the diversity calculations

plantDiversity_2013_overall <- diversity(BaylesPlantDiv_aggregated_2013_forDiversity_overall, index = "shannon", MARGIN = 1)
plantDiversity_2013_nativeONLY <- diversity(BaylesPlantDiv_aggregated_2013_forDiversity_nativeONLY, index = "shannon", MARGIN = 1)

plantDiversity_2014_overall <- diversity(BaylesPlantDiv_aggregated_2014_forDiversity_overall, index = "shannon", MARGIN = 1)
plantDiversity_2014_nativeONLY <- diversity(BaylesPlantDiv_aggregated_2014_forDiversity_nativeONLY, index = "shannon", MARGIN = 1)

plantDiversity_2016_overall <- diversity(BaylesPlantDiv_aggregated_2016_forDiversity_overall, index = "shannon", MARGIN = 1)
plantDiversity_2016_nativeONLY <- diversity(BaylesPlantDiv_aggregated_2016_forDiversity_nativeONLY, index = "shannon", MARGIN = 1)

# Need to put together the experimental information with the diversity values. This is the same for 2013 and 2014, but diff. for 2016.

fullPlantDiversity_2013_2014 <- data.frame("Plot" = rep(rep(seq(1,36,1), each = 2),2), "plotLoc" = rep(rep(c("C","R"),18),2), "overallDiversity" = c(plantDiversity_2013_overall,plantDiversity_2014_overall),
                                      "nativeDiversityONLY" = c(plantDiversity_2013_nativeONLY,plantDiversity_2014_nativeONLY), "AMFTreatment" = rep(rep(Bayles_plotTreatments$AMF, each = 2),2),
                                      "BiocharTreatment" = rep(rep(Bayles_plotTreatments$Biochar, each = 2),2), "Block" = rep(rep(Bayles_plotTreatments$Block, each = 2),2), "Year" = c(rep(2013,72),rep(2014,72)))

fullPlantDiversity_2013 <- data.frame("Plot" = rep(seq(1,36,1), each = 2), "plotLoc" = rep(c("C","R"),18), "overallDiversity" = plantDiversity_2013_overall,
                                           "nativeDiversityONLY" = plantDiversity_2013_nativeONLY, "AMFTreatment" = rep(Bayles_plotTreatments$AMF, each = 2),
                                           "BiocharTreatment" = rep(Bayles_plotTreatments$Biochar, each = 2), "Block" = rep(Bayles_plotTreatments$Block, each = 2), "Year" = rep(2013,72))

fullPlantDiversity_2014 <- data.frame("Plot" = rep(seq(1,36,1), each = 2), "plotLoc" = rep(c("C","R"),18), "overallDiversity" = plantDiversity_2014_overall,
                                           "nativeDiversityONLY" = plantDiversity_2014_nativeONLY, "AMFTreatment" = rep(Bayles_plotTreatments$AMF, each = 2),
                                           "BiocharTreatment" = rep(Bayles_plotTreatments$Biochar, each = 2), "Block" = rep(Bayles_plotTreatments$Block, each = 2), "Year" = rep(2014,72))



fullPlantDiversity_2016 <- data.frame("Plot" = Bayles_plantDiversity_2016$plotNum, "plotLoc" = rep("C",30), "overallDiversity" = plantDiversity_2016_overall,
                                      "nativeDiversityONLY" = plantDiversity_2016_nativeONLY, "AMFTreatment" = Bayles_plotTreatments$AMF[1:30],
                                      "BiocharTreatment" = Bayles_plotTreatments$Biochar[1:30], "Block" = Bayles_plotTreatments$Block[1:30], "Year" = rep("2016",30))

# 2013, 2014
# ================

# Overall diversity
# -----------------
# Plot location x AMF inoculation interaction
full_2013_2014_overall <- lmer(overallDiversity ~ AMFTreatment * BiocharTreatment * plotLoc + (1 | Year/Block/Plot), data = fullPlantDiversity_2013_2014)
print(lmerTest::anova(full_2013_2014_overall, type = 3, ddf="Kenward-Roger"))

#=========================
# Create the plot for this significant interaction (Figure S2 in the manuscript)
noAMF_ring <- data.frame(meanHPrime = mean(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "N" & fullPlantDiversity_2013_2014$plotLoc == "R"], na.rm = TRUE),
                         seHPrime = sd(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "N" & fullPlantDiversity_2013_2014$plotLoc == "R"], na.rm = TRUE)/sqrt(sum(!is.na(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "N" & fullPlantDiversity_2013_2014$plotLoc == "R"]))))

noAMF_center <- data.frame(meanHPrime = mean(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "N" & fullPlantDiversity_2013_2014$plotLoc == "C"], na.rm = TRUE),
                           seHPrime = sd(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "N" & fullPlantDiversity_2013_2014$plotLoc == "C"], na.rm = TRUE)/sqrt(sum(!is.na(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "N" & fullPlantDiversity_2013_2014$plotLoc == "C"]))))

yesAMF_ring <- data.frame(meanHPrime = mean(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "Y" & fullPlantDiversity_2013_2014$plotLoc == "R"], na.rm = TRUE),
                          seHPrime = sd(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "Y" & fullPlantDiversity_2013_2014$plotLoc == "R"], na.rm = TRUE)/sqrt(sum(!is.na(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "Y" & fullPlantDiversity_2013_2014$plotLoc == "R"]))))

yesAMF_center <- data.frame(meanHPrime = mean(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "Y" & fullPlantDiversity_2013_2014$plotLoc == "C"], na.rm = TRUE),
                            seHPrime = sd(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "Y" & fullPlantDiversity_2013_2014$plotLoc == "C"], na.rm = TRUE)/sqrt(sum(!is.na(fullPlantDiversity_2013_2014$overallDiversity[fullPlantDiversity_2013_2014$AMFTreatment == "Y" & fullPlantDiversity_2013_2014$plotLoc == "C"]))))

HPrime_YesAMF <- rbind(yesAMF_center, yesAMF_ring)
rownames(srp_YesAMF) <- c("Middle" ,"Margin")

HPrime_NoAMF <- rbind(noAMF_center, noAMF_ring)
rownames(srp_NoAMF) <- c("Middle" ,"Margin")

interactPlot <- ggplot() + 
    geom_line(aes(x = c(0.95,1.95), y = HPrime_NoAMF$meanHPrime), color = "#DE77AE", size = 2) + 
    geom_line(aes(x = c(1.05,2.05), y = HPrime_YesAMF$meanHPrime), color = "#7FBC41", size = 2) + 
    geom_errorbar(aes(x = c(0.95,1.95),
                      ymin = HPrime_NoAMF$meanHPrime - HPrime_NoAMF$seHPrime, ymax = HPrime_NoAMF$meanHPrime + HPrime_NoAMF$seHPrime),
                  width = 0.05) +
    geom_errorbar(aes(x = c(1.05,2.05),
                      ymin = HPrime_YesAMF$meanHPrime - HPrime_YesAMF$seHPrime, ymax = HPrime_YesAMF$meanHPrime + HPrime_YesAMF$seHPrime),
                  width = 0.05) + 
    
    geom_point(aes(x = c(0.95,1.95), y = HPrime_NoAMF$meanHPrime), color = "#DE77AE", size = 5) +
    geom_point(aes(x = c(1.05,2.05), y = HPrime_YesAMF$meanHPrime), color = "#7FBC41", size = 5) + 
    theme_classic() + 
    xlab("Location in plot") + ylab("Shannon diversity index (H')") + 
    scale_fill_manual(breaks = "+ AM fungi", "- AM fungi") + 
    scale_x_discrete(limits = c(1,2), labels = c("Middle", "Margin")) + 
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))# + ylim(0,1)

print(interactPlot)
#=========================

# Native plant diversity higher in middle compared to margins of plots regardless of AM fungal inoc.
# --------------------
full_2013_2014_native <- lmer(nativeDiversityONLY ~ AMFTreatment * BiocharTreatment * plotLoc + (1 | Year/Block/Plot), data = fullPlantDiversity_2013_2014)
print(lmerTest::anova(full_2013_2014_native, type = 3, ddf="Kenward-Roger"))
# --------------------

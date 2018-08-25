# Using linear mixed effect model to test for differences in nurse plant growth over all years at the same time.
# Growth is relative to each plant's height at planting and is log10-transformed.
# This produces the output in Table S3 and the LMER results noted in Fig. 2 (asterisks after
# plant species names)

if(!(require("lme4"))){
    install.packages("lme4")
} 

if(!(require("lmerTest"))){
    install.packages("lmerTest")
} 

library(lme4)
library(lmerTest)

# SAS - like contrasts (dropping last factor)
options(contrasts=c(factor = "contr.SAS", ordered = "contr.poly"))

nurseGrowth <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_NursePlantGrowthAllYears.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Add a column for the block number for each plot based on the plot number. This is working correctly
nurseGrowth$blockNum <- floor((nurseGrowth$plotNum - 1)/6) + 1

nurseGrowth_2013 <- nurseGrowth[nurseGrowth$yearSampled == 2013,]

nurseGrowth_2014 <- nurseGrowth[nurseGrowth$yearSampled == 2014,]

nurseGrowth_2015 <- nurseGrowth[nurseGrowth$yearSampled == 2015,]

# Only And ger growth in 2016
nurseGrowth_2016 <- nurseGrowth[nurseGrowth$yearSampled == 2016,]

nurseGrowth_2013_noNA <- nurseGrowth_2013[which(!is.na(nurseGrowth_2013$heightNormToPlantingHeight)),]
# Only keep the measurements for measurement dates after the initial (planting) date, when all relative heights
# are 1.
nurseGrowth_2013_noNA_afterPlantingMeasures <- nurseGrowth_2013_noNA[nurseGrowth_2013_noNA$dateSampled != 41407,]

nurseGrowth_2014_noNA <- nurseGrowth_2014[which(!is.na(nurseGrowth_2014$heightNormToPlantingHeight)),]

nurseGrowth_2015_noNA <- nurseGrowth_2015[which(!is.na(nurseGrowth_2015$heightNormToPlantingHeight)),]

nurseGrowth_2016_noNA <- nurseGrowth_2016[which(!is.na(nurseGrowth_2016$heightNormToPlantingHeight)),]

# All years combined. This is 2013-2014 for most nurse plants due to low survival after that,
# and 2013 - 2015 for And ger and Les cap.
# ===================

dataToUse_allYears <- rbind(nurseGrowth_2013_noNA_afterPlantingMeasures, nurseGrowth_2014_noNA, nurseGrowth_2015_noNA)

dataToUse_allYears$plotNum <- as.factor(dataToUse_allYears$plotNum)
dataToUse_allYears$yearSampled <- as.factor(dataToUse_allYears$yearSampled)

# Only for use with And ger
dataToUse_allYears_through2016 <- rbind(nurseGrowth_2013_noNA_afterPlantingMeasures, nurseGrowth_2014_noNA, nurseGrowth_2015_noNA, nurseGrowth_2016_noNA)

dataToUse_allYears_through2016$plotNum <- as.factor(dataToUse_allYears_through2016$plotNum)
dataToUse_allYears_through2016$yearSampled <- as.factor(dataToUse_allYears_through2016$yearSampled)

for(nurseSpecies in unique(dataToUse_allYears$speciesName)){
    
    print(paste0(" @@@ This is for: ", nurseSpecies))
    
    # Remove the June 29, 2013 measures (day 41454) for Liatris, Desmodium, Lespedeza, Allium, and Echinacea because purposefully
    # did not look for them to measure during that measurement. This was the only time this was the case.
    if(nurseSpecies %in% c("Lia spi", "Les cap")){
        nurseSpeciesData <- dataToUse_allYears[which(dataToUse_allYears$speciesName == nurseSpecies & dataToUse_allYears$dateSampled != 41454),]
    } else if(nurseSpecies %in% c("Ech pur", "All cer", "Pet pur", "Cor tri", "Des ill", "Sch sco")){
        # Only use the 2013,2014 measures for Ech pur, All cer, Dal pur, Des ill, Sch sco (2015 test for these is omitted below because the 
        # model fit fails due to lack of treatment representation.). Also only use 2013 and 2014 for Cor tri because couldn't confidently ascribe nurse plants in 2015.
        nurseSpeciesData <- dataToUse_allYears[which(dataToUse_allYears$speciesName == nurseSpecies & dataToUse_allYears$dateSampled != 41454 & dataToUse_allYears$yearSampled != 2015),]
    } else if (nurseSpecies == "And ger"){
        # Need to drop the 2016 measurements for And ger, otherwise model doesn't converge.
        nurseSpeciesData <- dataToUse_allYears_through2016[which(dataToUse_allYears_through2016$speciesName == nurseSpecies & dataToUse_allYears_through2016$yearSampled != 2016),]
    }    else{
        nurseSpeciesData <- dataToUse_allYears[dataToUse_allYears$speciesName == nurseSpecies,]
    }
    
    # Do not run model for Lia spi (very poor survival)
    if(nurseSpecies != "Lia spi"){
        #The random effect of year is set to help account for variability in the slope of the lines, not the intercepts
        # because from Fig 2, the intercepts are fairly similar for the plants across the years, but the slopes can vary
        # by quite a bit. 
        # log10 transformation makes huge difference in removing funnel shapes in residual plots (all models still converge)
        # Added interaction possibility with biochar to model, and then need to set REML = TRUE to keep all models converging.

        allYears <- lmer(log10(heightNormToPlantingHeight) ~  as.factor(AMFInoculation) * as.factor(biocharLevel) + (1+ yearSampled|blockNum/plotNum), data = nurseSpeciesData, na.action = na.omit, REML = TRUE)
        
        print(lmerTest::anova(allYears, type = 3, ddf="Kenward-Roger"))
    }
}
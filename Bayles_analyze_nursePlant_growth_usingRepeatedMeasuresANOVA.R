# Using repeated measures ANOVA to test for differences in nurse plant growth over 
# each of the first two growing seasons. 
# Growth is relative to each plant's height at planting and is log10-transformed.
# This produces the output in Table S4 and the repeated measures ANOVA results noted in Fig. 2
# (asterisks after each growing season label)

# All results for Liatris spicata (Lia spi) are unreliable for testing AMF inoculation effects
# due to low survival of non-inoculated plants (and marginally better for inoculated plants),
# and this is reflected in its labeling in Table S4 and Fig. 2. Survival across biochar treatments
# was generally better and allows for more representative comparisons (although no
# biochar treatment was a significant predictor.)

# This is a repeated measures anova (repeated measures of height within each year; the years analyzed separately)
# using the most basic error structure (the plot number) to control for the repeated measures

# SAS - like contrasts (dropping last factor)
options(contrasts=c(factor = "contr.SAS", ordered = "contr.poly"))

nurseGrowth <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_NursePlantGrowthAllYears.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

nurseGrowth_2013 <- nurseGrowth[nurseGrowth$yearSampled == 2013,]

nurseGrowth_2014 <- nurseGrowth[nurseGrowth$yearSampled == 2014,]

nurseGrowth_2015 <- nurseGrowth[nurseGrowth$yearSampled == 2015,]

# Only And ger growth in 2016
nurseGrowth_2016 <- nurseGrowth[nurseGrowth$yearSampled == 2016,]
# andgerGrowth_2016 <- andgerGrowth[andgerGrowth$yearSampled == 2016,]

nurseGrowth_2013_noNA <- nurseGrowth_2013[which(!is.na(nurseGrowth_2013$heightNormToPlantingHeight)),]
nurseGrowth_2013_noNA_afterPlantingMeasures <- nurseGrowth_2013_noNA[nurseGrowth_2013_noNA$dateSampled != 41407,]

nurseGrowth_2014_noNA <- nurseGrowth_2014[which(!is.na(nurseGrowth_2014$heightNormToPlantingHeight)),]

nurseGrowth_2015_noNA <- nurseGrowth_2015[which(!is.na(nurseGrowth_2015$heightNormToPlantingHeight)),]

nurseGrowth_2016_noNA <- nurseGrowth_2016[which(!is.na(nurseGrowth_2016$heightNormToPlantingHeight)),]


dataToUse_allYears <- rbind(nurseGrowth_2013_noNA_afterPlantingMeasures, nurseGrowth_2014_noNA, nurseGrowth_2015_noNA)

dataToUse_allYears$plotNum <- as.factor(dataToUse_allYears$plotNum)
dataToUse_allYears$yearSampled <- as.factor(dataToUse_allYears$yearSampled)

dataToUse_allYears$AMFInoculation <- as.factor(dataToUse_allYears$AMFInoculation)
dataToUse_allYears$biocharLevel <- as.factor(dataToUse_allYears$biocharLevel)

# Need to add column that is days since Jan 1st of the sampling year
# 2013
serialDatesForConversion_2013 <- dataToUse_allYears$dateSampled[dataToUse_allYears$yearSampled == 2013]
Jan_1_2013_date <- as.Date("2013-01-01")
datesForConversion_2013 <- as.Date(serialDatesForConversion_2013, origin = "1900-01-01")
# Gives a difftime class object
convertedDatesSince_Jan_1_2013 <- datesForConversion_2013 - Jan_1_2013_date
convertedDatesSince_Jan_1_2013 <- as.numeric(convertedDatesSince_Jan_1_2013)

#2014
serialDatesForConversion_2014 <- dataToUse_allYears$dateSampled[dataToUse_allYears$yearSampled == 2014]
Jan_1_2014_date <- as.Date("2014-01-01")
datesForConversion_2014 <- as.Date(serialDatesForConversion_2014, origin = "1900-01-01")
# Gives a difftime class object
convertedDatesSince_Jan_1_2014 <- datesForConversion_2014 - Jan_1_2014_date
convertedDatesSince_Jan_1_2014 <- as.numeric(convertedDatesSince_Jan_1_2014)

#2015
serialDatesForConversion_2015 <- dataToUse_allYears$dateSampled[dataToUse_allYears$yearSampled == 2015]
Jan_1_2015_date <- as.Date("2015-01-01")
datesForConversion_2015 <- as.Date(serialDatesForConversion_2015, origin = "1900-01-01")
# Gives a difftime class object
convertedDatesSince_Jan_1_2015 <- datesForConversion_2015 - Jan_1_2015_date
convertedDatesSince_Jan_1_2015 <- as.numeric(convertedDatesSince_Jan_1_2015)

dataToUse_allYears$numDaysSinceJan1 <- as.factor(c(convertedDatesSince_Jan_1_2013,convertedDatesSince_Jan_1_2014, convertedDatesSince_Jan_1_2015))

dataToUse_allYears_through2016$numDaysSinceJan1 <- as.factor(c(convertedDatesSince_Jan_1_2013,convertedDatesSince_Jan_1_2014, convertedDatesSince_Jan_1_2015, convertedDatesSince_Jan_1_2016))


# 2013 only
# ===================

print("###### 2013 ONLY ###########")

dataToUse_2013 <- nurseGrowth_2013_noNA_afterPlantingMeasures

dataToUse_2013$plotNum <- as.factor(dataToUse_2013$plotNum)
dataToUse_2013$yearSampled <- as.factor(dataToUse_2013$yearSampled)

for(nurseSpecies in unique(dataToUse_2013$speciesName)){
    
    print(paste0(" @@@ This is for: ", nurseSpecies))
    
    
    # Remove the June 29, 2013 measures (day 41454) for Liatris, Desmodium, Lespedeza, Allium, and Echinacea because purposefully
    # did not look for them to measure during that measurement. This was the only time this was the case.
    if(nurseSpecies %in% c("Lia spi", "Des ill", "Les cap", "Ech pur", "All cer")){
        nurseSpeciesData <- dataToUse_2013[which(dataToUse_2013$speciesName == nurseSpecies & dataToUse_2013$dateSampled != 41454),]
    } else{
        
        nurseSpeciesData <- dataToUse_2013[dataToUse_2013$speciesName == nurseSpecies,]
    }
    
    allYears_anova <- aov(log10(heightNormToPlantingHeight) ~  AMFInoculation * biocharLevel + Error(plotNum), data = nurseSpeciesData)
    
    print("****2013*****")
    print(summary(allYears_anova, type = 3))
    #print(model.tables(allYears_anova,"means"),digits=3)
    
    
}

# 2014 only NOTE - LIMITED POWER FOR SOME SP. B/C HIGH MORTALITY REFLECTED WITH 
# LOW RESID DF VALUES. Many with enough power are still sign. for AMF, Cor is for biochar too (grows better with no biochar; but no interact)
# ===================

print("###### 2014 ONLY ###########")

dataToUse_2014 <- nurseGrowth_2014_noNA

dataToUse_2014$plotNum <- as.factor(dataToUse_2014$plotNum)
dataToUse_2014$yearSampled <- as.factor(dataToUse_2014$yearSampled)

for(nurseSpecies in unique(dataToUse_2014$speciesName)){
    
    print(paste0(" @@@ This is for: ", nurseSpecies))
    
    nurseSpeciesData <- dataToUse_2014[dataToUse_2014$speciesName == nurseSpecies,]
    
    allYears_anova <- aov(log10(heightNormToPlantingHeight) ~  AMFInoculation * biocharLevel + Error(plotNum), data = nurseSpeciesData)
    
    print("****2014*****")
    print(summary(allYears_anova, type = 3))
    #print(model.tables(allYears_anova,"means"),digits=3)
    
    
}

# 2015 only; No sign. interacts (other than little blue based on very few remaining plants) so using main effects; NOTE - LIMITED POWER FOR SOME SP. B/C HIGH MORTALITY REFLECTED WITH 
# LOW RESID DF VALUES.  AMF still sign for And ger.
# ===================

print("###### 2015 ONLY ###########")

dataToUse_2015 <- nurseGrowth_2015_noNA

dataToUse_2015$plotNum <- as.factor(dataToUse_2015$plotNum)
dataToUse_2015$yearSampled <- as.factor(dataToUse_2015$yearSampled)

for(nurseSpecies in unique(dataToUse_2015$speciesName)){
    
    
    # Need to limit the comparisons due to high mortality, otherwise tests can't run.
    # 061618 The Sch sco test runs, and returns a sign. biochar x AMF interaction, but there's only 1 plant without
    # AMF and one plant without biochar at the second measurement of 2015, so I'm not trusting these results - the figure
    # with the Sch sco survival conveys this message very well, and is the main result for Sch sco - that biochar x AMF affected its
    # survival
    if(nurseSpecies != "All cer" & nurseSpecies != "Ech pur" & nurseSpecies != "Pet pur" & nurseSpecies != "Cor tri" & nurseSpecies != "Lia spi"){
        
        print(paste0(" @@@ This is for: ", nurseSpecies))
        
        nurseSpeciesData <- dataToUse_2015[dataToUse_2015$speciesName == nurseSpecies,]
        
        allYears_anova <- aov(log10(heightNormToPlantingHeight) ~  AMFInoculation * biocharLevel + Error(plotNum), data = nurseSpeciesData)
        
        print("****2015*****")
        print(summary(allYears_anova, type = 3))
        #print(model.tables(allYears_anova,"means"),digits=3)
        
    }
    
}


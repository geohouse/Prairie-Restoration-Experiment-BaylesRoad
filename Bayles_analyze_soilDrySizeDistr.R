# This analyzes the soil size distribution data to give the results reported in Table S5 and
# the graph that is Fig. 5A in the manuscript. The soil fractions were measured using 
# 8mm, 4mm, 2mm, 1mm, 0.5mm, and 0.25mm 
# meshes. For these analyses, the soil fractions > 1mm and those <= 1mm were combined 
# to give two groups. This 1mm cutoff between groups corresponds to the size class used
# for the water stable aggregates analysis. The mass of these size fractions is first
# converted to a proportion of the total sample mass, and then logit-transformed 
# before analysis. This tests using the >1mm size group, but because there are only two
# groups, results should be identical to those if using the <= 1mm size group instead.

# Install required packages if necessary.
if(!(require("lme4"))){
    install.packages("lme4")
} 

if(!(require("lmerTest"))){
    install.packages("lmerTest")
} 

if(!(require("ggplot2"))){
    install.packages("ggplot2")
} 

library(lme4)
library(lmerTest)
library(ggplot2)

# SAS - like contrasts (dropping last factor) - this only affects the contrast tested for the biochar.
options(contrasts=c(factor = "contr.SAS", ordered = "contr.poly"))

soilDistrData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/BaylesBiochar_drySizeDist_2013_2014.csv", header = TRUE, sep = ",")

# These are the two size groups used
mass_gt_1mm <- soilDistrData$mass_8mm + soilDistrData$mass_4mm + soilDistrData$mass_2mm
mass_lt_1mm <- soilDistrData$mass_1mm + soilDistrData$mass_0.5mm + soilDistrData$mass_0.25mm + soilDistrData$mass_fine

# Convert mass to proportion
prop_gt_1mm <- mass_gt_1mm/soilDistrData$totalMass
prop_lt_1mm <- mass_lt_1mm/soilDistrData$totalMass

# Logit-transform the proportion
logit_gt_1mm <- log10(prop_gt_1mm/(1-prop_gt_1mm))
logit_lt_1mm <- log10(prop_lt_1mm/(1-prop_lt_1mm))

# Enter these transformed values to the data frame used for the analysis.
soilDistrData$mass_gt_1mm <- mass_gt_1mm
soilDistrData$prop_gt_1mm <- prop_gt_1mm
soilDistrData$logit_gt_1mm <- logit_gt_1mm

soilDistrData$mass_lt_1mm <- mass_lt_1mm
soilDistrData$prop_lt_1mm <- prop_lt_1mm
soilDistrData$logit_lt_1mm <- logit_lt_1mm

soilDistrData$sampleNum <- seq(1,nrow(soilDistrData),1)


soilDistrData$block <- as.factor(soilDistrData$block)
soilDistrData$plot <- as.factor(soilDistrData$plot)

#-----------------
# For both years. 
print(">= 1mm logit for both growing seasons")
soilDrySizeTest_gt1_logit <- lmer(logit_gt_1mm ~ AMFInoc * biochar * plotLoc + (1 | year/block/plot), data = soilDistrData)
print(lmerTest::anova(soilDrySizeTest_gt1_logit, type = 3, ddf="Kenward-Roger"))
#-----------------

# Test just 2013 and 2014 separately

soilDistrData_2013 <- soilDistrData[soilDistrData$year == 2013,]
soilDistrData_2014 <- soilDistrData[soilDistrData$year == 2014,]

# ------------------
# 2013
print("> 1mm logit 2013 only")
soilDrySizeTest_gt1_logit_2013 <- lmer(logit_gt_1mm ~ AMFInoc * biochar * plotLoc + (1 | block/plot), data = soilDistrData_2013)
print(lmerTest::anova(soilDrySizeTest_gt1_logit_2013, type = 3, ddf="Kenward-Roger"))
# ------------------

# ------------------
print("> 1mm logit 2014 only")
soilDrySizeTest_gt1_logit_2014 <- lmer(logit_gt_1mm ~ AMFInoc * biochar * plotLoc + (1 | block/plot), data = soilDistrData_2014)
print(lmerTest::anova(soilDrySizeTest_gt1_logit_2014, type = 3, ddf="Kenward-Roger"))
# ------------------

# This is the plotter to make the plot in Figure 5A (plot location x biochar interaction
# for the first growing season)

noBiochar_ring <- data.frame(meanGt1mm = mean(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "N" & soilDistrData_2013$plotLoc == "Ring"], na.rm = TRUE),
                             seGt1mm = sd(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "N" & soilDistrData_2013$plotLoc == "Ring"], na.rm = TRUE)/sqrt(sum(!is.na(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "N" & soilDistrData_2013$plotLoc == "Ring"]))))

noBiochar_center <-  data.frame(meanGt1mm = mean(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "N" & soilDistrData_2013$plotLoc == "Center"], na.rm = TRUE),
                                seGt1mm = sd(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "N" & soilDistrData_2013$plotLoc == "Center"], na.rm = TRUE)/sqrt(sum(!is.na(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "N" & soilDistrData_2013$plotLoc == "Center"]))))

lowBiochar_ring <- data.frame(meanGt1mm = mean(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "L" & soilDistrData_2013$plotLoc == "Ring"], na.rm = TRUE),
                              seGt1mm = sd(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "L" & soilDistrData_2013$plotLoc == "Ring"], na.rm = TRUE)/sqrt(sum(!is.na(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "L" & soilDistrData_2013$plotLoc == "Ring"]))))

lowBiochar_center <-  data.frame(meanGt1mm = mean(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "L" & soilDistrData_2013$plotLoc == "Center"], na.rm = TRUE),
                                 seGt1mm = sd(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "L" & soilDistrData_2013$plotLoc == "Center"], na.rm = TRUE)/sqrt(sum(!is.na(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "L" & soilDistrData_2013$plotLoc == "Center"]))))

highBiochar_ring <- data.frame(meanGt1mm = mean(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "H" & soilDistrData_2013$plotLoc == "Ring"], na.rm = TRUE),
                               seGt1mm = sd(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "H" & soilDistrData_2013$plotLoc == "Ring"], na.rm = TRUE)/sqrt(sum(!is.na(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "H" & soilDistrData_2013$plotLoc == "Ring"]))))

highBiochar_center <-  data.frame(meanGt1mm = mean(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "H" & soilDistrData_2013$plotLoc == "Center"], na.rm = TRUE),
                                  seGt1mm = sd(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "H" & soilDistrData_2013$plotLoc == "Center"], na.rm = TRUE)/sqrt(sum(!is.na(soilDistrData_2013$logit_gt_1mm[soilDistrData_2013$biochar == "H" & soilDistrData_2013$plotLoc == "Center"]))))

logit_noBiochar <- rbind(noBiochar_center, noBiochar_ring)
rownames(logit_noBiochar) <- c("Middle", "Margin")

logit_lowBiochar <- rbind(lowBiochar_center, lowBiochar_ring)
rownames(logit_lowBiochar) <- c("Middle", "Margin")

logit_highBiochar <- rbind(highBiochar_center, highBiochar_ring)
rownames(logit_highBiochar) <- c("Middle", "Margin")

interactPlot <- ggplot() + 
    geom_line(aes(x = c(0.95,1.95), y = logit_noBiochar$meanGt1mm), color = "#D0D0D0", size = 2) +
    geom_line(aes(x = c(1,2), y = logit_lowBiochar$meanGt1mm), color = "#999999", size = 2) + 
    geom_line(aes(x = c(1.05,2.05), y = logit_highBiochar$meanGt1mm), color = "#000000", size = 2) + 
    geom_errorbar(aes(x = c(0.95,1.95),
                      ymin = logit_noBiochar$meanGt1mm - logit_noBiochar$seGt1mm, ymax = logit_noBiochar$meanGt1mm + logit_noBiochar$seGt1mm),
                  width = 0.05) + 
    geom_errorbar(aes(x = c(1,2),
                      ymin = logit_lowBiochar$meanGt1mm - logit_lowBiochar$seGt1mm, ymax = logit_lowBiochar$meanGt1mm + logit_lowBiochar$seGt1mm),
                  width = 0.05) +
    geom_errorbar(aes(x = c(1.05,2.05),
                      ymin = logit_highBiochar$meanGt1mm - logit_highBiochar$seGt1mm, ymax = logit_highBiochar$meanGt1mm + logit_highBiochar$seGt1mm),
                  width = 0.05) +
    
    geom_point(aes(x = c(0.95,1.95), y = logit_noBiochar$meanGt1mm), color = "#D0D0D0", size = 5) +
    
    
    geom_point(aes(x = c(1,2), y = logit_lowBiochar$meanGt1mm), color = "#999999", size = 5) +
    
    
    geom_point(aes(x = c(1.05,2.05), y = logit_highBiochar$meanGt1mm), color = "#000000", size = 5) +
    
    theme_classic() +
    xlab("Location in plot") + ylab("Fraction of aggregates (logit)") +
    scale_fill_manual(breaks = "No biochar", "Low biochar", "High biochar") +
    scale_x_discrete(limits = c(1,2), labels = c("Middle", "Margin")) +
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))# + ylim(0,1)


print(interactPlot)



# This analyzes the soil size distribution data to give the results reported in Table S5 and
# the graph that is Fig. 5A in the manuscript. The soil fractions were measured using 
# 8mm, 4mm, 2mm, 1mm, 0.5mm, and 0.25mm 
# meshes. For these analyses, the soil fractions > 1mm and those <= 1mm were combined 
# to give two groups. This 1mm cutoff between groups corresponds to the size class used
# for the water stable aggregates analysis. The mass of these size fractions is first
# converted to a proportion of the total sample mass, and then logit-transformed 
# before analysis.

# Using the logit transformed proportions with a > 1mm group and a <= 1mm group to correspond to the 
# WSA size fractions.

# Problem with changing ddf between 2013 and 2014 datasets fixed 112717.
# This is all caused by instability of the Satterthwate approx of df in the lmerTest::anova (I don't
# know if the Satterthwate approx is inherently unstable or whether this is a bug (more likely))
# With shuffling of data, the denominator df can change wildly (and it shouldn't).
# when use ddf="Kenward-Roger" in the anova call, the ddf becomes stable across shuffling, and 
# is consistent between years. This doesn't change p values much (another reason it seems to be a 
# bug), but can have a dramatic change on ddf. 112717.

# Install required packages if necessary.
if(!(require("lme4"))){
    install.packages("lme4")
} 

if(!(require("lmerTest"))){
    install.packages("lmerTest")
} 

library(lme4)
library(lmerTest)

# Default R contrasts (dropping first factor) 
options(contrasts = c(factor = "contr.treatment", ordered = "contr.treatment"))

# SAS - like contrasts (dropping last factor) - this only affects the contrast tested for the biochar.
options(contrasts=c(factor = "contr.SAS", ordered = "contr.poly"))

soilDistrData <- read.table("~/Box Sync/R_code/Bayles_biochar/BaylesBiochar_drySizeDist_2013_2014.csv", header = TRUE, sep = ",")

mass_gt_1mm <- soilDistrData$mass_8mm + soilDistrData$mass_4mm + soilDistrData$mass_2mm

mass_lt_1mm <- soilDistrData$mass_1mm + soilDistrData$mass_0.5mm + soilDistrData$mass_0.25mm + soilDistrData$mass_fine

prop_gt_1mm <- mass_gt_1mm/soilDistrData$totalMass
prop_lt_1mm <- mass_lt_1mm/soilDistrData$totalMass

logit_gt_1mm <- log10(prop_gt_1mm/(1-prop_gt_1mm))
logit_lt_1mm <- log10(prop_lt_1mm/(1-prop_lt_1mm))

soilDistrData$mass_gt_1mm <- mass_gt_1mm
soilDistrData$prop_gt_1mm <- prop_gt_1mm
soilDistrData$logit_gt_1mm <- logit_gt_1mm

soilDistrData$mass_lt_1mm <- mass_lt_1mm
soilDistrData$prop_lt_1mm <- prop_lt_1mm
soilDistrData$logit_lt_1mm <- logit_lt_1mm

soilDistrData$prop_fine <- soilDistrData$mass_fine/soilDistrData$totalMass
soilDistrData$prop_025mm <- soilDistrData$mass_0.25mm/soilDistrData$totalMass
soilDistrData$prop_05mm <- soilDistrData$mass_0.5mm/soilDistrData$totalMass
soilDistrData$prop_1mm <- soilDistrData$mass_1mm/soilDistrData$totalMass
soilDistrData$prop_2mm <- soilDistrData$mass_2mm/soilDistrData$totalMass
soilDistrData$prop_4mm <- soilDistrData$mass_4mm/soilDistrData$totalMass
soilDistrData$prop_8mm <- soilDistrData$mass_8mm/soilDistrData$totalMass

soilDistrData$logit_fine <- log10(soilDistrData$prop_fine/(1-soilDistrData$prop_fine))
soilDistrData$logit_025mm <- log10(soilDistrData$prop_025mm/(1-soilDistrData$prop_025mm))
soilDistrData$logit_05mm <- log10(soilDistrData$prop_05mm/(1-soilDistrData$prop_05mm))
soilDistrData$logit_1mm <- log10(soilDistrData$prop_1mm/(1-soilDistrData$prop_1mm))
soilDistrData$logit_2mm <- log10(soilDistrData$prop_2mm/(1-soilDistrData$prop_2mm))
soilDistrData$logit_4mm <- log10(soilDistrData$prop_4mm/(1-soilDistrData$prop_4mm))
soilDistrData$logit_8mm <- log10(soilDistrData$prop_8mm/(1-soilDistrData$prop_8mm))

soilDistrData$sampleNum <- seq(1,nrow(soilDistrData),1)

# soilDistrData$year <- as.factor(soilDistrData$year)
soilDistrData$block <- as.factor(soilDistrData$block)
soilDistrData$plot <- as.factor(soilDistrData$plot)

soilDrySizeTest_gt1 <- lmer(mass_gt_1mm ~ AMFInoc + biochar + plotLoc + (1 | year/block), data = soilDistrData)
print(">= 1mm mass")
print(lmerTest::summary(soilDrySizeTest_gt1))

soilDrySizeTest_gt1_prop <- lmer(prop_gt_1mm ~ AMFInoc + biochar + plotLoc + (1 | year/block), data = soilDistrData)
print(">= 1mm proportion")
print(lmerTest::summary(soilDrySizeTest_gt1_prop))

# Sign. greater > 1mm size fractions with H biochar compared to L biochar (p = 0.0466), but not compared to N biochar (0.0859); 
# No diff between L and N biochar (p = 0.7812). No sign. interaction.
soilDrySizeTest_gt1_logit <- lmer(logit_gt_1mm ~ AMFInoc + biochar + plotLoc + (1 | year/block), data = soilDistrData)
print(">= 1mm logit")
print(lmerTest::summary(soilDrySizeTest_gt1_logit))


soilDrySizeTest_lt1 <- lmer(mass_lt_1mm ~ AMFInoc + biochar + plotLoc + (1 | year/block), data = soilDistrData)
print("< 1mm mass")
print(lmerTest::summary(soilDrySizeTest_lt1))
plot(soilDrySizeTest_lt1)

soilDrySizeTest_lt1_prop <- lmer(prop_lt_1mm ~ AMFInoc + biochar + plotLoc + (1 | year/block), data = soilDistrData)
print("< 1mm proportion")
print(lmerTest::summary(soilDrySizeTest_lt1_prop))
plot(soilDrySizeTest_lt1_prop)

# Sign. less < 1mm size fractions with H biochar compared to L biochar (p = 0.0466), but not compared to N biochar (0.0859); 
# No diff between L and N biochar (p = 0.7812). No sign. interaction.
# soilDrySizeTest_lt1_logit <- lmer(logit_lt_1mm ~ AMFInoc + biochar + plotLoc + (1 | year/block), data = soilDistrData)
# print("< 1mm logit")
# print(lmerTest::summary(soilDrySizeTest_lt1_logit))
# print(lmerTest::anova(soilDrySizeTest_lt1_logit))
# plot(soilDrySizeTest_lt1_logit)

#-----------------
# For both years. 060518 Nothing significant; this is reported in Table S3 in the manuscript.
# Confirmed 081918
soilDrySizeTest_gt1_logit <- lmer(logit_gt_1mm ~ AMFInoc * biochar * plotLoc + (1 | year/block/plot), data = soilDistrData)
print(">= 1mm logit")
# Modified - makes no difference for this combined test.
print(lmerTest::anova(soilDrySizeTest_gt1_logit, type = 3, ddf="Kenward-Roger"))
#-----------------

print(lmerTest::summary(soilDrySizeTest_gt1_logit))
# Old.
# print(lmerTest::anova(soilDrySizeTest_gt1_logit, type = 3))
# Modified - makes no difference for this combined test.
print(lmerTest::anova(soilDrySizeTest_gt1_logit, type = 3, ddf="Kenward-Roger"))
plot(soilDrySizeTest_gt1_logit)


interaction.plot(as.factor(soilDistrData$AMFInoc),as.factor(soilDistrData$plotLoc),
                 response = soilDistrData$logit_gt_1mm, xlab = "AMF inoculation", ylab = "Dry size distr", main = "> 1mm")

interaction.plot(as.factor(soilDistrData$biochar),as.factor(soilDistrData$plotLoc),
                 response = soilDistrData$logit_gt_1mm, xlab = "Biochar", ylab = "Dry size distr", main = "> 1mm")


plot(logit_lt_1mm ~ biochar, data = soilDistrData)
plot(logit_gt_1mm ~ biochar, data = soilDistrData)


# Test just 2013 and 2014 separately

soilDistrData_2013 <- soilDistrData[soilDistrData$year == 2013,]
soilDistrData_2014 <- soilDistrData[soilDistrData$year == 2014,]

# ------------------
# 2013
# 060518 This is the model reported in the manuscript.
# Confirmed 081918
soilDrySizeTest_gt1_logit_2013 <- lmer(logit_gt_1mm ~ AMFInoc * biochar * plotLoc + (1 | block/plot), data = soilDistrData_2013)
print("> 1mm logit 2013 only")
# Modified. Didn't change results
print(lmerTest::anova(soilDrySizeTest_gt1_logit_2013, type = 3, ddf="Kenward-Roger"))
# ------------------


print(lmerTest::summary(soilDrySizeTest_gt1_logit_2013))
# Old.
# print(lmerTest::anova(soilDrySizeTest_gt1_logit_2013, type = 3))
# Modified. Didn't change results
print(lmerTest::anova(soilDrySizeTest_gt1_logit_2013, type = 3, ddf="Kenward-Roger"))

plot(soilDrySizeTest_gt1_logit_2013)

interaction.plot(as.factor(soilDistrData_2013$AMFInoc),as.factor(soilDistrData_2013$plotLoc),
                                   response = soilDistrData_2013$logit_gt_1mm, xlab = "AMF inoculation", ylab = "Dry size distr", main = "> 1mm")

interaction.plot(as.factor(soilDistrData_2013$biochar),as.factor(soilDistrData_2013$plotLoc),
                 response = soilDistrData_2013$logit_gt_1mm, xlab = "Biochar", ylab = "Dry size distr", main = "> 1mm")


# Nothing here.
# soilDrySizeTest_lt1_logit_2013 <- lmer(logit_lt_1mm ~ AMFInoc * biochar + plotLoc + (1 | block), data = soilDistrData_2013)
# print("< 1mm logit 2013 only")
# print(lmerTest::summary(soilDrySizeTest_lt1_logit_2013))
# plot(soilDrySizeTest_lt1_logit_2013)


# OLD - Sign. interaction between biochar and AMF. For H biochar, greater fraction with - AMF, less with + AMF.
# For L Biochar, less fraction with -AMF, more with + AMF. This interaction only appears with the summary(), not quite with anova()

# Sign. greater > 1mm fraction with H biochar compared to N biochar (p = 0.00648) and L biochar (p = 0.01588) (No diff between L and N biochar)

# *************
# For some reason, the 2014 data with plot as a nested random factor doesn't take it into account when
# calculating the df or the P value (I have no idea why - I've tried everything I can think of and 
# somehow it's just the data subsetting step for whatever reason doesn't change the dfs on the second
# half of the data when plot it included as a random var.).

# One way to *ugly* get around this, computing the regularized beta function in Wolfram Alpha like this:
# betaregularized(55/(55+1*0.08), 55/2,1/2)
# Where 55 here is the denom df,
# 1 is the numer df,
# 0.08 is the F value. Replace these as necess to calc the corresponding p value.

# ***************
# But I'm not doing this here. I'll report the DF that R gives, even if it seems wrong. (The p-val)

# ------------------
# 060518 This is the model reported in the manuscript
# Confirmed 081918.
soilDrySizeTest_gt1_logit_2014 <- lmer(logit_gt_1mm ~ AMFInoc * biochar * plotLoc + (1 | block/plot), data = soilDistrData_2014)
print("> 1mm logit 2014 only")
# Modified - this keeps the ddf the same between 2013 and 2014.
print(lmerTest::anova(soilDrySizeTest_gt1_logit_2014, type = 3, ddf="Kenward-Roger"))
# ------------------

print(lmerTest::summary(soilDrySizeTest_gt1_logit_2014))
# Old - this is critical NOT to use
# print(lmerTest::anova(soilDrySizeTest_gt1_logit_2014, type = 3))
# Modified - this keeps the ddf the same between 2013 and 2014.
print(lmerTest::anova(soilDrySizeTest_gt1_logit_2014, type = 3, ddf="Kenward-Roger"))

print(lmerTest::lsmeansLT(soilDrySizeTest_gt1_logit_2014))
plot(soilDrySizeTest_gt1_logit_2014)

boxplot(soilDistrData_2014$logit_gt_1mm ~ soilDistrData_2014$biochar)



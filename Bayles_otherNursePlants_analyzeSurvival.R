# Analyze Bayles nurse plant loss data for species other than Andropogon gerardii 
# (Cor tri, Les cap, Des ill, Dal pur, Sch sco, Ech pur) 
# in glm mixed model (blocks as random effects), with a binomial
# (logit) link to fit a logistic regression. Survival coded as 0 (survived) or 1 (died) as of the last
# measure in 2014 (survival very poor for these after that, and couldn't tell nurse Cor tri confidently
# from seeded plants). Also testing survival after 2013 (using 365 days as cutoff between 2013 and 2014)

# This also generates the interaction plot shown as Fig 4 in the manuscript.

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

# SAS - like contrasts (dropping last factor)
options(contrasts=c(factor = "contr.SAS", ordered = "contr.poly"))

dalpur_LossData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_Dalpur_lossDates_through2014_ForR.csv", header = TRUE, sep = ",")
desill_LossData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_Desill_lossDates_through2014_ForR.csv", header = TRUE, sep = ",")
echpur_LossData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_Echpur_lossDates_through2014_ForR.csv", header = TRUE, sep = ",")
cortri_LossData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_Cortri_lossDates_through2014_ForR.csv", header = TRUE, sep = ",")
schsco_LossData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_Schsco_lossDates_through2014_ForR.csv", header = TRUE, sep = ",")
lescap_LossData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_Lescap_lossDates_through2014_ForR.csv", header = TRUE, sep = ",")
allcer_LossData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_Allcer_lossDates_through2014_ForR.csv", header = TRUE, sep = ",")
liaspi_LossData <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_Liaspi_lossDates_through2014_ForR.csv", header = TRUE, sep = ",")

dalpur_LossData$mortalityAtEnd <- ifelse(!is.na(dalpur_LossData$lossDate_daysSince_05_13_2013),1,0)
desill_LossData$mortalityAtEnd <- ifelse(!is.na(desill_LossData$lossDate_daysSince_05_13_2013),1,0)
echpur_LossData$mortalityAtEnd <- ifelse(!is.na(echpur_LossData$lossDate_daysSince_05_13_2013),1,0)
cortri_LossData$mortalityAtEnd <- ifelse(!is.na(cortri_LossData$lossDate_daysSince_05_13_2013),1,0)
schsco_LossData$mortalityAtEnd <- ifelse(!is.na(schsco_LossData$lossDate_daysSince_05_13_2013),1,0)
lescap_LossData$mortalityAtEnd <- ifelse(!is.na(lescap_LossData$lossDate_daysSince_05_13_2013),1,0)
allcer_LossData$mortalityAtEnd <- ifelse(!is.na(allcer_LossData$lossDate_daysSince_05_13_2013),1,0)
liaspi_LossData$mortalityAtEnd <- ifelse(!is.na(liaspi_LossData$lossDate_daysSince_05_13_2013),1,0)

dalpur_LossData$mortalityAtEnd_2013 <- ifelse(!is.na(dalpur_LossData$lossDate_daysSince_05_13_2013) & dalpur_LossData$lossDate_daysSince_05_13_2013 < 365,1,0)
desill_LossData$mortalityAtEnd_2013 <- ifelse(!is.na(desill_LossData$lossDate_daysSince_05_13_2013) & desill_LossData$lossDate_daysSince_05_13_2013 < 365,1,0)
echpur_LossData$mortalityAtEnd_2013 <- ifelse(!is.na(echpur_LossData$lossDate_daysSince_05_13_2013) & echpur_LossData$lossDate_daysSince_05_13_2013 < 365,1,0)
cortri_LossData$mortalityAtEnd_2013 <- ifelse(!is.na(cortri_LossData$lossDate_daysSince_05_13_2013) & cortri_LossData$lossDate_daysSince_05_13_2013 < 365,1,0)
schsco_LossData$mortalityAtEnd_2013 <- ifelse(!is.na(schsco_LossData$lossDate_daysSince_05_13_2013) & schsco_LossData$lossDate_daysSince_05_13_2013 < 365,1,0)
lescap_LossData$mortalityAtEnd_2013 <- ifelse(!is.na(lescap_LossData$lossDate_daysSince_05_13_2013) & lescap_LossData$lossDate_daysSince_05_13_2013 < 365,1,0)
allcer_LossData$mortalityAtEnd_2013 <- ifelse(!is.na(allcer_LossData$lossDate_daysSince_05_13_2013) & allcer_LossData$lossDate_daysSince_05_13_2013 < 365,1,0)
liaspi_LossData$mortalityAtEnd_2013 <- ifelse(!is.na(liaspi_LossData$lossDate_daysSince_05_13_2013) & liaspi_LossData$lossDate_daysSince_05_13_2013 < 365,1,0)

# Change biochar from 3 level to binary (+/- any biochar in the soil) for these tests so the output
# is easily comparable to the ANOVA tables from the plant growth and the soil structure

dalpur_LossData$binaryBiochar <- ifelse(dalpur_LossData$biocharLevel == 'H' |  dalpur_LossData$biocharLevel == 'L', 'Y', 'N')
desill_LossData$binaryBiochar <- ifelse(desill_LossData$biocharLevel == 'H' |  desill_LossData$biocharLevel == 'L', 'Y', 'N')
echpur_LossData$binaryBiochar <- ifelse(echpur_LossData$biocharLevel == 'H' |  echpur_LossData$biocharLevel == 'L', 'Y', 'N')
cortri_LossData$binaryBiochar <- ifelse(cortri_LossData$biocharLevel == 'H' |  cortri_LossData$biocharLevel == 'L', 'Y', 'N')
schsco_LossData$binaryBiochar <- ifelse(schsco_LossData$biocharLevel == 'H' |  schsco_LossData$biocharLevel == 'L', 'Y', 'N')
lescap_LossData$binaryBiochar <- ifelse(lescap_LossData$biocharLevel == 'H' |  lescap_LossData$biocharLevel == 'L', 'Y', 'N')
allcer_LossData$binaryBiochar <- ifelse(allcer_LossData$biocharLevel == 'H' |  allcer_LossData$biocharLevel == 'L', 'Y', 'N')
liaspi_LossData$binaryBiochar <- ifelse(liaspi_LossData$biocharLevel == 'H' |  liaspi_LossData$biocharLevel == 'L', 'Y', 'N')

# Survival through 2014 - Sign interaction for Schsco
# ======================
print("Through 2014")

print("Cor tri")
# No sign. interaction. Testing main effects.
cortri_fullSurvival <- glmer(mortalityAtEnd ~ as.factor(AMFInoculation) + as.factor(binaryBiochar) + (1 | block), 
                             data = cortri_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(cortri_fullSurvival))

print("Des ill")
# No sign. interaction. Testing main effects.
desill_fullSurvival <- glmer(mortalityAtEnd ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = desill_LossData, family = "binomial", control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(desill_fullSurvival))

print("Dal pur")
dalpur_fullSurvival <- glmer(mortalityAtEnd ~ as.factor(AMFInoculation) + as.factor(binaryBiochar) + (1 | block), 
                             data = dalpur_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(dalpur_fullSurvival))

print("Ech pur")
# No sign. interaction. Testing main effects.
echpur_fullSurvival <- glmer(mortalityAtEnd ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = echpur_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(echpur_fullSurvival))

print("Les cap")
# No sign. interaction. Testing main effects.
lescap_fullSurvival <- glmer(mortalityAtEnd ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = lescap_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(lescap_fullSurvival))

print("Sch sco")

schsco_fullSurvival <- glmer(mortalityAtEnd ~ as.factor(AMFInoculation) * as.factor(biocharLevel) + (1 | block), 
                             data = schsco_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(schsco_fullSurvival))

# ====================
# Create a plot of the interaction between biochar treatment and AMF 
# inoculation affecting Sch sco survival (Fig 4 in manuscript.)

noBio_NoAMF <- data.frame(propSurviving = sum(schsco_LossData$biocharLevel == "N" & 
                                                  schsco_LossData$mortalityAtEnd == 0 & 
                                                  schsco_LossData$AMFInoculation == "N")/6,
                          se = sd(schsco_LossData$biocharLevel == "N" & 
                                      schsco_LossData$mortalityAtEnd == 0 & 
                                      schsco_LossData$AMFInoculation == "N")/sqrt(6))

LBio_NoAMF <- data.frame(propSurviving = sum(schsco_LossData$biocharLevel == "L" & 
                                                 schsco_LossData$mortalityAtEnd == 0 & 
                                                 schsco_LossData$AMFInoculation == "N")/6,
                         se = sd(schsco_LossData$biocharLevel == "L" & 
                                     schsco_LossData$mortalityAtEnd == 0 & 
                                     schsco_LossData$AMFInoculation == "N")/sqrt(6))

HBio_NoAMF <- data.frame(propSurviving = sum(schsco_LossData$biocharLevel == "H" & 
                                                 schsco_LossData$mortalityAtEnd == 0 & 
                                                 schsco_LossData$AMFInoculation == "N")/6,
                         se = sd(schsco_LossData$biocharLevel == "H" & 
                                     schsco_LossData$mortalityAtEnd == 0 & 
                                     schsco_LossData$AMFInoculation == "N")/sqrt(6))

noBio_YesAMF <- data.frame(propSurviving = sum(schsco_LossData$biocharLevel == "N" & 
                                                   schsco_LossData$mortalityAtEnd == 0 & 
                                                   schsco_LossData$AMFInoculation == "Y")/6,
                           se = sd(schsco_LossData$biocharLevel == "N" & 
                                       schsco_LossData$mortalityAtEnd == 0 & 
                                       schsco_LossData$AMFInoculation == "Y")/sqrt(6))

LBio_YesAMF <- data.frame(propSurviving = sum(schsco_LossData$biocharLevel == "L" & 
                                                  schsco_LossData$mortalityAtEnd == 0 & 
                                                  schsco_LossData$AMFInoculation == "Y")/6,
                          se = sd(schsco_LossData$biocharLevel == "L" & 
                                      schsco_LossData$mortalityAtEnd == 0 & 
                                      schsco_LossData$AMFInoculation == "Y")/sqrt(6))

HBio_YesAMF <- data.frame(propSurviving = sum(schsco_LossData$biocharLevel == "H" & 
                                                  schsco_LossData$mortalityAtEnd == 0 & 
                                                  schsco_LossData$AMFInoculation == "Y")/6,
                          se = sd(schsco_LossData$biocharLevel == "H" & 
                                      schsco_LossData$mortalityAtEnd == 0 & 
                                      schsco_LossData$AMFInoculation == "Y")/sqrt(6))

survival_YesAMF <- rbind(noBio_YesAMF, LBio_YesAMF, HBio_YesAMF)
rownames(survival_YesAMF) <- c("No_biochar", "Low_biochar", "High_biochar")

survival_NoAMF <- rbind(noBio_NoAMF, LBio_NoAMF, HBio_NoAMF)
rownames(survival_NoAMF) <- c("No_biochar", "Low_biochar", "High_biochar")

interactPlot <- ggplot() + geom_errorbar(aes(x = c(1,2,3),
                                             ymin = survival_NoAMF$propSurviving - survival_NoAMF$se, ymax = survival_NoAMF$propSurviving + survival_NoAMF$se),
                                         width = 0.1) + 
    geom_line(aes(x = c(1,2,3), y = survival_NoAMF$propSurviving), color = "#DE77AE", size = 2) + 
    geom_point(aes(x = c(1,2,3), y = survival_NoAMF$propSurviving), color = "#DE77AE", size = 5) +
    geom_errorbar(aes(x = c(1,2,3),
                      ymin = survival_YesAMF$propSurviving - survival_YesAMF$se, ymax = survival_YesAMF$propSurviving + survival_YesAMF$se),
                  width = 0.1) + 
    geom_line(aes(x = c(1,2,3), y = survival_YesAMF$propSurviving), color = "#7FBC41", size = 2) + 
    geom_point(aes(x = c(1,2,3), y = survival_YesAMF$propSurviving), color = "#7FBC41", size = 5) + 
    theme_classic() + 
    xlab("Biochar treatment") + ylab("Proportion surviving") + 
    scale_fill_manual(breaks = "+ AM fungi", "- AM fungi") + 
    scale_x_discrete(limits = c(1,2,3), labels = c("No", "Low", "High")) + 
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) + ylim(0,1)


print(interactPlot)

#========================

print("All cer")
# No sign. interaction. Testing main effects.
allcer_fullSurvival <- glmer(mortalityAtEnd ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = allcer_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(allcer_fullSurvival))

print("Lia spi")
# No sign. interaction. Testing main effects.
liaspi_fullSurvival <- glmer(mortalityAtEnd ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = liaspi_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(liaspi_fullSurvival))

# ======================


# Survival through 2013; Dal pur mortality is higher with no biochar in the soil - p = 0.05
# ======================
print("Through 2013 ONLY")

print("Cor tri")
# No sign. interaction. Testing main effects.
cortri_2013Survival <- glmer(mortalityAtEnd_2013 ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = cortri_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(cortri_2013Survival))

print("Des ill")
# No sign. interaction. Testing main effects.
desill_2013Survival <- glmer(mortalityAtEnd_2013 ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = desill_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(desill_2013Survival))

print("Dal pur")
# Can't test the interaction of AMF/biochar - model fails to converge. 
dalpur_2013Survival <- glmer(mortalityAtEnd_2013 ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = dalpur_LossData, family = binomial("logit"), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
print(summary(dalpur_2013Survival))

numDiedHBiochar <- sum(dalpur_LossData$mortalityAtEnd_2013[dalpur_LossData$biocharLevel == "H"])
numDiedNBiochar <- sum(dalpur_LossData$mortalityAtEnd_2013[dalpur_LossData$biocharLevel == "N"])
numDiedLBiochar <- sum(dalpur_LossData$mortalityAtEnd_2013[dalpur_LossData$biocharLevel == "L"])

print("Ech pur")
# No sign. interaction. Testing main effects.
echpur_2013Survival <- glmer(mortalityAtEnd_2013 ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = echpur_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(echpur_2013Survival))

print("Les cap")
# No sign. interaction. Testing main effects.
lescap_2013Survival <- glmer(mortalityAtEnd_2013 ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                             data = lescap_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(lescap_2013Survival))


print("Sch sco") 
schsco_2013Survival <- glmer(mortalityAtEnd_2013 ~ as.factor(AMFInoculation) * as.factor(biocharLevel) + (1 | block), 
                             data = schsco_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                               optCtrl=list(maxfun=100000)))
print(summary(schsco_2013Survival))


print("All cer") 
# No sign. interaction. Testing main effects.
    allcer_2013Survival <- glmer(mortalityAtEnd_2013 ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                                 data = allcer_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                                   optCtrl=list(maxfun=100000)))
print(summary(allcer_2013Survival))

print("Lia spi")
# No sign. interaction. Testing main effects.
    liaspi_2013Survival <- glmer(mortalityAtEnd_2013 ~ as.factor(AMFInoculation) + as.factor(biocharLevel) + (1 | block), 
                                 data = liaspi_LossData, family = binomial(link = "logit"), control=glmerControl(optimizer="bobyqa",
                                                                                                   optCtrl=list(maxfun=100000)))
print(summary(liaspi_2013Survival))
# ======================


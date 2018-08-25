# Analyze Bayles big bluestem loss data in glm mixed model (blocks as random effects), with a binomial
# (logit) link to fit a logistic regression. Survival coded as 0 (survived) or 1 (died) as of the last
# measure in fall 2016 (082116)

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

andgerLossData <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_Andger_lossDatesForR.csv", header = TRUE, sep = ",")

andgerLossData$block <- as.factor(rep(seq(1,6,1), each = 6))
andgerLossData$mortalityAtEnd <- ifelse(!is.na(andgerLossData$lossDate_daysSince_05_13_2013),1,0)

# Combined AMF and Biochar
# Testing interaction doesn't converge. This is testing main effects
andger_fullSurvival <- glmer(mortalityAtEnd ~ AMFInoculation + biocharLevel + (1 | block), 
                             data = andgerLossData, family = binomial(link = "logit"), 
                             control=glmerControl(optimizer="bobyqa",
                                optCtrl=list(maxfun=100000)))

summary(andger_fullSurvival)

numDiedYAMF <- sum(andgerLossData$mortalityAtEnd[andgerLossData$AMFInoculation == "Y"])
numDiedNAMF <- sum(andgerLossData$mortalityAtEnd[andgerLossData$AMFInoculation == "N"])

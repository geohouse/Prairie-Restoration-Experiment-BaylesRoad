# This analyzes the water stable aggregates (WSA) data to give the results 
# reported in Table S5. WSA data measured from blocks 2-6 for 2013 (first growing season)
# and from blocks 1-3 for 2014 (second growing season).

# Install required packages if necessary.
if(!(require("lme4"))){
    install.packages("lme4")
} 

if(!(require("lmerTest"))){
    install.packages("lmerTest")
} 

library(lme4)
library(lmerTest)

# SAS - like contrasts (dropping last factor) - this only affects the contrast tested for the biochar.
options(contrasts=c(factor = "contr.SAS", ordered = "contr.poly"))

WSAData_2013 <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/BaylesBiochar_2013_WSAMeasures.csv", header = TRUE, sep = ",")
WSAData_2014 <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/BaylesBiochar_2014_WSAMeasures.csv", header = TRUE, sep = ",")

WSA_combined <- rbind(WSAData_2013, WSAData_2014)

# Logit transform the data
logit_WSA <- log10(WSA_combined$WSA/(1-WSA_combined$WSA))

WSA_combined$logit_WSA <- logit_WSA

WSA_forTest_2013 <- WSA_combined[WSA_combined$year == "2013",]

WSA_forTest_2014 <- WSA_combined[WSA_combined$year == "2014",]

#-----------------
# Both years - nothing sign.
print("Both years.")
WSA_test <- lmer(logit_WSA ~ AMFInoc * biochar * plotLoc + (1 | year/block/plot), data = WSA_combined)
print(lmerTest::anova(WSA_test, type = 3, ddf="Kenward-Roger"))
#-----------------

#-----------------
# 2013 only - nothing sign.
print("2013")
WSA_test_2013 <- lmer(logit_WSA ~ AMFInoc * biochar * plotLoc + (1 | block/plot), data = WSA_forTest_2013)
# Modified
print(lmerTest::anova(WSA_test_2013, type = 3, ddf="Kenward-Roger"))
#-----------------

#-----------------------
# 2014 only - nothing sign.
print("2014")
WSA_test_2014 <- lmer(logit_WSA ~ AMFInoc * biochar * plotLoc +  (1 | block/plot), data = WSA_forTest_2014)
print(lmerTest::anova(WSA_test_2014, type = 3, ddf="Kenward-Roger"))
#----------------------
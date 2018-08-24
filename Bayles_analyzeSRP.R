# This analyzes the soluble reactive phosphorus (SRP) flux data from near the soil surface
# measured at the end of the first growing season, and generates an interaction plot
# (Fig. 5B in the manuscript)

library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)

# SAS - like contrasts (dropping last factor) - this only affects the contrast tested for the biochar.
options(contrasts=c(factor = "contr.SAS", ordered = "contr.poly"))

SRPData <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_SRP_Jan2014_forR.csv", header = TRUE, sep = ",")

SRPData$AMFInoc <- as.factor(SRPData$AMFInoc)
SRPData$biocharAmt <- as.factor(SRPData$biocharAmt)
SRPData$plotLocation <- as.factor(SRPData$plotLocation)
SRPData$block <- as.factor(SRPData$block)
SRPData$plotNum <- as.factor(SRPData$plotNum)

SRP_noNA <- SRPData[which(!is.na(SRPData$SRP_flux_ugPerSqCm_over29Days)),]

SRP_mixModel <- lmer(sqrt(SRP_flux_ugPerSqCm_over29Days) ~ AMFInoc * biocharAmt * plotLocation + (1 | block/plotLocation), data = SRP_noNA)

# Printing analysis results.
print(lmerTest::anova(SRP_mixModel, type = 3, ddf="Kenward-Roger"))

# Plotter (Fig 5B)
#--------

noAMF_ring <- data.frame(meanSRPFlux = mean(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "N" & SRPData$plotLocation == "R"], na.rm = TRUE),
                         seSRPFlux = sd(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "N" & SRPData$plotLocation == "R"], na.rm = TRUE)/sqrt(sum(!is.na(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "N" & SRPData$plotLocation == "R"]))))

noAMF_center <- data.frame(meanSRPFlux = mean(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "N" & SRPData$plotLocation == "C"], na.rm = TRUE),
                           seSRPFlux = sd(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "N" & SRPData$plotLocation == "C"], na.rm = TRUE)/sqrt(sum(!is.na(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "N" & SRPData$plotLocation == "C"]))))

yesAMF_ring <- data.frame(meanSRPFlux = mean(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "Y" & SRPData$plotLocation == "R"], na.rm = TRUE),
                          seSRPFlux = sd(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "Y" & SRPData$plotLocation == "R"], na.rm = TRUE)/sqrt(sum(!is.na(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "Y" & SRPData$plotLocation == "R"]))))

yesAMF_center <- data.frame(meanSRPFlux = mean(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "Y" & SRPData$plotLocation == "C"], na.rm = TRUE),
                            seSRPFlux = sd(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "Y" & SRPData$plotLocation == "C"], na.rm = TRUE)/sqrt(sum(!is.na(SRPData$SRP_flux_ugPerSqCm_over29Days[SRPData$AMFInoc == "Y" & SRPData$plotLocation == "C"]))))


srp_YesAMF <- rbind(yesAMF_center, yesAMF_ring)
rownames(srp_YesAMF) <- c("Middle" ,"Margin")

srp_NoAMF <- rbind(noAMF_center, noAMF_ring)
rownames(srp_NoAMF) <- c("Middle" ,"Margin")


interactPlot <- ggplot() + 
    geom_line(aes(x = c(0.95,1.95), y = srp_NoAMF$meanSRPFlux), color = "#DE77AE", size = 2) + 
    geom_line(aes(x = c(1.05,2.05), y = srp_YesAMF$meanSRPFlux), color = "#7FBC41", size = 2) + 
    geom_errorbar(aes(x = c(0.95,1.95),
                      ymin = srp_NoAMF$meanSRPFlux - srp_NoAMF$seSRPFlux, ymax = srp_NoAMF$meanSRPFlux + srp_NoAMF$seSRPFlux),
                  width = 0.05) +
    geom_errorbar(aes(x = c(1.05,2.05),
                      ymin = srp_YesAMF$meanSRPFlux - srp_YesAMF$seSRPFlux, ymax = srp_YesAMF$meanSRPFlux + srp_YesAMF$seSRPFlux),
                  width = 0.05) + 
    
    geom_point(aes(x = c(0.95,1.95), y = srp_NoAMF$meanSRPFlux), color = "#DE77AE", size = 5) +
    geom_point(aes(x = c(1.05,2.05), y = srp_YesAMF$meanSRPFlux), color = "#7FBC41", size = 5) + 
    theme_classic() + 
    xlab("Location in plot") + ylab(expression(paste("SRP flux ", expression(mu), "g cm2"))) + 
    scale_fill_manual(breaks = "+ AM fungi", "- AM fungi") + 
    scale_x_discrete(limits = c(1,2), labels = c("Middle", "Margin")) + 
    theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))# + ylim(0,1)

print(interactPlot)


# This is a plotter for the Bayles nurse plant root colonization (arbuscules)
# This plot is Fig S1 in the manuscript.

rootScores <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_nurseRootColonizationScores_forR.csv", header = TRUE, sep = ",")

# order by functional group then by species name
rootScores <- rootScores[order(rootScores$FunctionalGroup, rootScores$Species),]

rootScores$speciesPlusTreatment <- as.factor(paste(rootScores$Species, rootScores$sterileOrInoc, sep = "_"))

fillColors <- rep(c("#7FBC41","#DE77AE"),9)

# Saving at 6.5W and 4H seems to look alright.

boxplot(rootScores$arbuscColonization ~ rootScores$speciesPlusTreatment, boxfill = fillColors,
        xaxt = "n", cex.axis = 1.2, ylab = "Fraction of arbuscular colonization in roots ", cex.lab = 1)


axLims <- par("usr")

axis(side = 1, las = 1, at = seq(1.5,17.5,2), labels = unique(rootScores$Species), cex.axis = 0.5)

abline(v=seq(2.5,16.5,2))
# This is the combination boxplot script to plot the nurse plant growth
# by AMF inoculation treatment shown in Fig. 2.

# Asterisk after plant name is sign. effect of AMF across the years plotted for it (Linear mixed effects model).
# Asterisk after the year label is sign. AMF effect for growth that year (repeated measures ANOVA).

# Saves pretty well as 10x 10

if(!(require("dplyr"))){
    install.packages("dplyr")
} 

if(!(require("lubridate"))){
    install.packages("lubridate")
} 

if(!(require("ggplot2"))){
    install.packages("ggplot2")
} 

if(!(require("tidyr"))){
    install.packages("tidyr")
} 

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

nurseGrowth <- read.table("~/Prairie-Restoration-Experiment-BaylesRoad-master/Bayles_NursePlantGrowthAllYears.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

serialMay1Date_2013 <- 41395
serialMay1Date_2014 <- 41760
serialMay1Date_2015 <- 42125
serialMay1Date_2016 <- 42491
serialMay1Date_2017 <- 42856

daySpacer <- 30

num2013Days <- yday("2013-11-01") - yday("2013-05-01")
num2014Days <- yday("2014-09-01") - yday("2014-05-01")
num2015Days <- yday("2015-09-01") - yday("2015-05-01")
num2016Days <- yday("2016-09-01") - yday("2016-05-01")

setDatesForPlotting <- function(currRow){
    sampleDate <- as.numeric(currRow[5])
    sampleYear <- currRow[6]
    
    if(sampleYear == 2013){
        return(sampleDate - serialMay1Date_2013) 
    } else if(sampleYear == 2014){
        return(sampleDate - serialMay1Date_2014 + num2013Days + daySpacer)
    } else if(sampleYear == 2015){
        return(sampleDate - serialMay1Date_2015 + num2013Days + daySpacer + num2014Days + daySpacer)
    } else if(sampleYear == 2016){
        return(sampleDate - serialMay1Date_2016 + num2013Days + daySpacer + num2014Days + daySpacer + num2015Days + daySpacer)
    }
}

recodedDatesForGraph <- apply(X = nurseGrowth, MARGIN = 1, FUN = setDatesForPlotting)

nurseGrowth$recodedDateForPlot <- recodedDatesForGraph

# split into nurse plant data sets
# Remove the June 29, 2013 measures (day 41454) for Liatris, Desmodium, Lespedeza, Allium, and Echinacea because purposefully
# did not look for them to measure during that measurement. This was the only time this was the case.

andgerGrowth <- nurseGrowth[nurseGrowth$speciesName == "And ger",]
# only 2013 and 2014 for Ech pur
echpurGrowth <- nurseGrowth[nurseGrowth$speciesName == "Ech pur" & nurseGrowth$yearSampled %in% c(2013, 2014) & nurseGrowth$dateSampled != 41454,]
# only 2013 and 2014 for All cer
allcerGrowth <- nurseGrowth[nurseGrowth$speciesName == "All cer" & nurseGrowth$yearSampled %in% c(2013, 2014) & nurseGrowth$dateSampled != 41454,]

lescapGrowth <- nurseGrowth[nurseGrowth$speciesName == "Les cap" & nurseGrowth$dateSampled != 41454,]
# only 2013 and 2014 for Dal pur
dalpurGrowth <- nurseGrowth[nurseGrowth$speciesName == "Pet pur" & nurseGrowth$yearSampled %in% c(2013, 2014),]

desillGrowth <- nurseGrowth[nurseGrowth$speciesName == "Des ill" & nurseGrowth$dateSampled != 41454,]
schscoGrowth <- nurseGrowth[nurseGrowth$speciesName == "Sch sco",]
# only 2013 and 2014 for Cor tri
cortriGrowth <- nurseGrowth[nurseGrowth$speciesName == "Cor tri" & nurseGrowth$yearSampled %in% c(2013, 2014),]
# Also measures missing for days 41483, 41507, and 41557 for Liatris
liaspiGrowth <- nurseGrowth[nurseGrowth$speciesName == "Lia spi" & nurseGrowth$yearSampled %in% c(2013, 2014) & !(nurseGrowth$dateSampled %in% c(41454,41483,41507,41557)),]

par(mfrow = c(3,3))

boxWidth <- 8

axisHolder_2013_2014_2015_2016 <- c(yday("2013-05-01") - yday("2013-05-01"), 
                                    yday("2013-06-01") - yday("2013-05-01"), 
                                    yday("2013-07-01") - yday("2013-05-01"), 
                                    yday("2013-08-01") - yday("2013-05-01"), 
                                    yday("2013-09-01") - yday("2013-05-01"), 
                                    yday("2013-10-01") - yday("2013-05-01"), 
                                    yday("2013-11-01") - yday("2013-05-01"),
                                    yday("2014-05-01") - yday("2014-05-01") + num2013Days + daySpacer,
                                    yday("2014-06-01") - yday("2014-05-01") + num2013Days + daySpacer,
                                    yday("2014-07-01") - yday("2014-05-01") + num2013Days + daySpacer,
                                    yday("2014-08-01") - yday("2014-05-01") + num2013Days + daySpacer,
                                    yday("2014-09-01") - yday("2014-05-01") + num2013Days + daySpacer,
                                    yday("2015-05-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                                    yday("2015-06-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                                    yday("2015-07-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                                    yday("2015-08-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                                    yday("2015-09-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                                    yday("2016-05-01") - yday("2016-05-01") + num2013Days + daySpacer + num2014Days + daySpacer + num2015Days + daySpacer,
                                    yday("2016-06-01") - yday("2016-05-01") + num2013Days + daySpacer + num2014Days + daySpacer + num2015Days + daySpacer,
                                    yday("2016-07-01") - yday("2016-05-01") + num2013Days + daySpacer + num2014Days + daySpacer + num2015Days + daySpacer,
                                    yday("2016-08-01") - yday("2016-05-01") + num2013Days + daySpacer + num2014Days + daySpacer + num2015Days + daySpacer,
                                    yday("2016-09-01") - yday("2016-05-01") + num2013Days + daySpacer + num2014Days + daySpacer + num2015Days + daySpacer)


axisHolder_2013_2014 <- c(yday("2013-05-01") - yday("2013-05-01"), 
                          yday("2013-06-01") - yday("2013-05-01"), 
                          yday("2013-07-01") - yday("2013-05-01"), 
                          yday("2013-08-01") - yday("2013-05-01"), 
                          yday("2013-09-01") - yday("2013-05-01"), 
                          yday("2013-10-01") - yday("2013-05-01"), 
                          yday("2013-11-01") - yday("2013-05-01"),
                          yday("2014-05-01") - yday("2014-05-01") + num2013Days + daySpacer,
                          yday("2014-06-01") - yday("2014-05-01") + num2013Days + daySpacer,
                          yday("2014-07-01") - yday("2014-05-01") + num2013Days + daySpacer,
                          yday("2014-08-01") - yday("2014-05-01") + num2013Days + daySpacer,
                          yday("2014-09-01") - yday("2014-05-01") + num2013Days + daySpacer)

axisHolder_2013_2014_2015 <- c(yday("2013-05-01") - yday("2013-05-01"), 
                               yday("2013-06-01") - yday("2013-05-01"), 
                               yday("2013-07-01") - yday("2013-05-01"), 
                               yday("2013-08-01") - yday("2013-05-01"), 
                               yday("2013-09-01") - yday("2013-05-01"), 
                               yday("2013-10-01") - yday("2013-05-01"), 
                               yday("2013-11-01") - yday("2013-05-01"),
                               yday("2014-05-01") - yday("2014-05-01") + num2013Days + daySpacer,
                               yday("2014-06-01") - yday("2014-05-01") + num2013Days + daySpacer,
                               yday("2014-07-01") - yday("2014-05-01") + num2013Days + daySpacer,
                               yday("2014-08-01") - yday("2014-05-01") + num2013Days + daySpacer,
                               yday("2014-09-01") - yday("2014-05-01") + num2013Days + daySpacer,
                               yday("2015-05-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                               yday("2015-06-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                               yday("2015-07-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                               yday("2015-08-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer,
                               yday("2015-09-01") - yday("2015-05-01") + num2013Days + daySpacer + num2014Days + daySpacer)



# passed to cex.main for setting the size of the titles
titleSize <- 2
# Used to set the size of the year labels in each panel.
yearLabelSize <- 1.5
# Used to set the size of the axes labels (months and relative heights)
axisLabelSize <- 1.3
# Used to scale the y axis title
yAxisTitleSize <- 1.5

# And ger
# ===============
par(mar=c(4,3.5,2,1.5))
boxplot(andgerGrowth$heightNormToPlantingHeight ~ andgerGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(andgerGrowth$recodedDateForPlot), main = expression(bold(italic("Andropogon gerardii*"))), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, ylim = c(0,22.2))

axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013*", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014*", pos = 4, cex = yearLabelSize)
text(354,yAxTextPlacement, "2015*", pos = 4, cex = yearLabelSize)
text(505,yAxTextPlacement, "2016^", pos = 4, cex = yearLabelSize)

boxplot(andgerGrowth$heightNormToPlantingHeight[andgerGrowth$AMFInoculation == "Y"] ~ andgerGrowth$dateSampled[andgerGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(andgerGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n")
boxplot(andgerGrowth$heightNormToPlantingHeight[andgerGrowth$AMFInoculation == "N"] ~ andgerGrowth$dateSampled[andgerGrowth$AMFInoculation == "N"],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(andgerGrowth$recodedDateForPlot) - 5,
        xaxt = "n", yaxt = "n")

axis(side = 1,at = axisHolder_2013_2014_2015_2016, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov","May","Jun","Jul","Aug","Sep","May","Jun","Jul","Aug","Sep",
                                                                       "May","Jun","Jul","Aug","Sep"), cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(v = num2013Days + daySpacer + num2014Days + (daySpacer/2), col = "black", lty = 1)
abline(v = num2013Days + daySpacer + num2014Days + daySpacer + num2015Days + (daySpacer/2), col = "black", lty = 1)
abline(h = 1, col = "red", lty = 2)
# ========================

# Sch sco
# ===============
boxplot(schscoGrowth$heightNormToPlantingHeight ~ schscoGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(schscoGrowth$recodedDateForPlot), main = expression(italic("Schizachyrium scoparium*")), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, xlim = c(-6.62,490),
        ylim = c(0,8))

axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013*", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014*", pos = 4, cex = yearLabelSize)
text(354,yAxTextPlacement, "2015^", pos = 4, cex = yearLabelSize)

boxplot(schscoGrowth$heightNormToPlantingHeight[schscoGrowth$AMFInoculation == "Y"] ~ schscoGrowth$dateSampled[schscoGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(schscoGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n", xlim = c(-6.62,490))

boxplot(schscoGrowth$heightNormToPlantingHeight[schscoGrowth$AMFInoculation == "N"] ~ schscoGrowth$dateSampled[schscoGrowth$AMFInoculation == "N"],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(schscoGrowth$recodedDateForPlot) - 5,
        xaxt = "n", yaxt = "n", xlim = c(-6.62,490))

axis(side = 1,at = axisHolder_2013_2014_2015, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov",
                                                                  "May","Jun","Jul","Aug","Sep",
                                                                  "May","Jun","Jul","Aug","Sep"), cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(v = num2013Days + daySpacer + num2014Days + (daySpacer/2), col = "black", lty = 1)

abline(h = 1, col = "red", lty = 2)
# ========================

# Dal pur
# ===============
boxplot(dalpurGrowth$heightNormToPlantingHeight ~ dalpurGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(dalpurGrowth$recodedDateForPlot) , main = expression(italic("Dalea purpurea")), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, ylim = c(0,8))
axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013*", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014*", pos = 4, cex = yearLabelSize)

boxplot(dalpurGrowth$heightNormToPlantingHeight[dalpurGrowth$AMFInoculation == "Y"] ~ dalpurGrowth$dateSampled[dalpurGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(dalpurGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n", ylim = c(0.1788166,5.6829345))
# Very poor survival in 2015; plotting just 2013 and 2014.
boxplot(dalpurGrowth$heightNormToPlantingHeight[dalpurGrowth$AMFInoculation == "N"] ~ dalpurGrowth$dateSampled[dalpurGrowth$AMFInoculation == "N"],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(dalpurGrowth$recodedDateForPlot) - 5,
        xaxt = "n", yaxt = "n", ylim = c(0.1788166,5.6829345))

axis(side = 1,at = axisHolder_2013_2014, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov",
                                                             "May","Jun","Jul","Aug","Sep"), cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(v = num2013Days + daySpacer + num2014Days + (daySpacer/2), col = "black", lty = 1)

abline(h = 1, col = "red", lty = 2)
# ========================

# Des ill
# ===============
boxplot(desillGrowth$heightNormToPlantingHeight ~ desillGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(desillGrowth$recodedDateForPlot), main = expression(italic("Desmodium illinoense*")), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, ylim = c(-1.53, 50.2), xlim = c(-6.62,490))

axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013*", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014*", pos = 4, cex = yearLabelSize)
text(354,yAxTextPlacement, "2015^", pos = 4, cex = yearLabelSize)

boxplot(desillGrowth$heightNormToPlantingHeight[desillGrowth$AMFInoculation == "Y"] ~ desillGrowth$dateSampled[desillGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(desillGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n", ylim = c(-1.53, 50.2), xlim = c(-6.62,490))

boxplot(desillGrowth$heightNormToPlantingHeight[desillGrowth$AMFInoculation == "N"] ~ desillGrowth$dateSampled[desillGrowth$AMFInoculation == "N"],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(desillGrowth$recodedDateForPlot) - 5,
        xaxt = "n", yaxt = "n", ylim = c(0, 50.2), xlim = c(-6.62,490))

axis(side = 1,at = axisHolder_2013_2014_2015, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov",
                                                                  "May","Jun","Jul","Aug","Sep",
                                                                  "May","Jun","Jul","Aug","Sep"), cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(v = num2013Days + daySpacer + num2014Days + (daySpacer/2), col = "black", lty = 1)

abline(h = 1, col = "red", lty = 2)
# ========================

# Les cap
# ===============
boxplot(lescapGrowth$heightNormToPlantingHeight ~ lescapGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(lescapGrowth$recodedDateForPlot), main = expression(italic("Lespedeza capitata*")), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, xlim = c(-6.62,490),
        ylim = c(0,22.2))

axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013*", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014", pos = 4, cex = yearLabelSize)
text(354,yAxTextPlacement, "2015", pos = 4, cex = yearLabelSize)

boxplot(lescapGrowth$heightNormToPlantingHeight[lescapGrowth$AMFInoculation == "Y"] ~ lescapGrowth$dateSampled[lescapGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(lescapGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n", xlim = c(-6.62,490))

boxplot(lescapGrowth$heightNormToPlantingHeight[lescapGrowth$AMFInoculation == "N"] ~ lescapGrowth$dateSampled[lescapGrowth$AMFInoculation == "N"],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(lescapGrowth$recodedDateForPlot) - 5,
        xaxt = "n", yaxt = "n", xlim = c(-6.62,490))

axis(side = 1,at = axisHolder_2013_2014_2015, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov",
                                                                  "May","Jun","Jul","Aug","Sep",
                                                                  "May","Jun","Jul","Aug","Sep"), cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(v = num2013Days + daySpacer + num2014Days + (daySpacer/2), col = "black", lty = 1)

abline(h = 1, col = "red", lty = 2)
# ========================

# All cer
# ===============
boxplot(allcerGrowth$heightNormToPlantingHeight ~ allcerGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(allcerGrowth$recodedDateForPlot) , main = expression(italic("Allium cernuum")), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, ylim = c(0,8))

axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014", pos = 4, cex = yearLabelSize)

boxplot(allcerGrowth$heightNormToPlantingHeight[allcerGrowth$AMFInoculation == "Y"] ~ allcerGrowth$dateSampled[allcerGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(allcerGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n", ylim = c(0.1788166,5.6829345))
# No non-AMF plants measured in the first measurement of 2014 so remove that
boxplot(allcerGrowth$heightNormToPlantingHeight[allcerGrowth$AMFInoculation == "N" & allcerGrowth$dateSampled != 41785] ~ allcerGrowth$dateSampled[allcerGrowth$AMFInoculation == "N"& allcerGrowth$dateSampled != 41785],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(allcerGrowth$recodedDateForPlot)[c(1:5,7:9)] - 5,
        xaxt = "n", yaxt = "n", ylim = c(0.1788166,5.6829345))

axis(side = 1,at = axisHolder_2013_2014, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov",
                                                             "May","Jun","Jul","Aug","Sep"), cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)

abline(h = 1, col = "red", lty = 2)
# ========================

# Cor tri
# ===============
boxplot(cortriGrowth$heightNormToPlantingHeight ~ cortriGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(cortriGrowth$recodedDateForPlot), main = expression(italic("Coreopsis tripteris*")), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, ylim = c(0,22.2))

axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013*", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014*", pos = 4, cex = yearLabelSize)

boxplot(cortriGrowth$heightNormToPlantingHeight[cortriGrowth$AMFInoculation == "Y"] ~ cortriGrowth$dateSampled[cortriGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(cortriGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n")
boxplot(cortriGrowth$heightNormToPlantingHeight[cortriGrowth$AMFInoculation == "N"] ~ cortriGrowth$dateSampled[cortriGrowth$AMFInoculation == "N"],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(cortriGrowth$recodedDateForPlot) - 5,
        xaxt = "n", yaxt = "n")

axis(side = 1,at = axisHolder_2013_2014, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov","May","Jun","Jul","Aug","Sep"),
     cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(h = 1, col = "red", lty = 2)
# ========================

# Ech pur
# ===============

boxplot(echpurGrowth$heightNormToPlantingHeight ~ echpurGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(echpurGrowth$recodedDateForPlot) , main = expression(italic("Echinacea purpurea")), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, ylim = c(0,22.2))
axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014", pos = 4, cex = yearLabelSize)
boxplot(echpurGrowth$heightNormToPlantingHeight[echpurGrowth$AMFInoculation == "Y"] ~ echpurGrowth$dateSampled[echpurGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(echpurGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n")
boxplot(echpurGrowth$heightNormToPlantingHeight[echpurGrowth$AMFInoculation == "N"] ~ echpurGrowth$dateSampled[echpurGrowth$AMFInoculation == "N"],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(echpurGrowth$recodedDateForPlot) - 5,
        xaxt = "n", yaxt = "n")

axis(side = 1,at = axisHolder_2013_2014, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov","May","Jun","Jul","Aug","Sep"),
     cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(h = 1, col = "red", lty = 2)
# ========================

# Lia spi
# ===============
boxplot(liaspiGrowth$heightNormToPlantingHeight ~ liaspiGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = unique(liaspiGrowth$recodedDateForPlot), main = expression(italic("Liatris spicata")), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize, ylim = c(0,8))

axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013^", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014^", pos = 4, cex = yearLabelSize)

boxplot(liaspiGrowth$heightNormToPlantingHeight[liaspiGrowth$AMFInoculation == "Y"] ~ liaspiGrowth$dateSampled[liaspiGrowth$AMFInoculation == "Y"],
        boxfill = "#7FBC41", add = TRUE, boxwex = boxWidth, at = unique(liaspiGrowth$recodedDateForPlot) + 5,
        xaxt = "n", yaxt = "n", ylim = c(0.1788166,5.6829345))
boxplot(liaspiGrowth$heightNormToPlantingHeight[liaspiGrowth$AMFInoculation == "N"] ~ liaspiGrowth$dateSampled[liaspiGrowth$AMFInoculation == "N"],
        boxfill = "#DE77AE", add = TRUE, boxwex = boxWidth, at = unique(liaspiGrowth$recodedDateForPlot)[c(1,2,4)] - 5,
        xaxt = "n", yaxt = "n", ylim = c(0.1788166,5.6829345))

axis(side = 1,at = axisHolder_2013_2014, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov","May","Jun","Jul","Aug","Sep"),
     cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(h = 1, col = "red", lty = 2)
# ========================

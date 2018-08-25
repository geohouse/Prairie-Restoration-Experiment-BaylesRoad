# This plots the growth of Coreopsis tripteris nurse-plants with different biochar
# treatments during the first two growing seasons to illustrate the significant biochar 
# effect on its growth. Growth is relative to planting height for each plant. 

# This is the plot in Figure 3 of the manuscript.

# Asterisk after the year label is sign. biochar effect for growth that year.

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

nurseGrowth <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_NursePlantGrowthAllYears.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

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

# only 2013 and 2014 for Cor tri
cortriGrowth <- nurseGrowth[nurseGrowth$speciesName == "Cor tri" & nurseGrowth$yearSampled %in% c(2013, 2014),]

boxWidth <- 6
boxOffset <- 7


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


# passed to cex.main for setting the size of the titles
titleSize <- 2
# Used to set the size of the year labels in each panel.
yearLabelSize <- 1.5
# Used to set the size of the axes labels (months and relative heights)
axisLabelSize <- 1.3
# Used to scale the y axis title
yAxisTitleSize <- 1.5


# OK to save as 6W x 4H.
# Cor tri
# ===============
boxplot(cortriGrowth$heightNormToPlantingHeight ~ cortriGrowth$dateSampled, boxfill = rgb(1,1,1, alpha = 1),
        border = rgb(1,1,1, alpha = 1), xaxt = "n", at = c(0,unique(cortriGrowth$recodedDateForPlot)[-1]), cex.main = titleSize,
        cex.axis = axisLabelSize, ylab = "Relative Growth", cex.lab = yAxisTitleSize)#, ylim = c(0,22.2))

axLims <- par("usr")
yAxSpan <- axLims[4] - axLims[3]
yAxTextPlacement <- axLims[4] - 0.05*yAxSpan
text(1,yAxTextPlacement, "2013", pos = 4, cex = yearLabelSize)
text(199,yAxTextPlacement, "2014*", pos = 4, cex = yearLabelSize)

boxplot(cortriGrowth$heightNormToPlantingHeight[cortriGrowth$biocharLevel == "N"] ~ cortriGrowth$dateSampled[cortriGrowth$biocharLevel == "N"],
        boxfill = "#FFFFFF", add = TRUE, boxwex = boxWidth, at = c(0,unique(cortriGrowth$recodedDateForPlot[cortriGrowth$biocharLevel == "N" & !is.na(cortriGrowth$heightNormToPlantingHeight)])[-1]) - boxOffset,
        xaxt = "n", yaxt = "n")
boxplot(cortriGrowth$heightNormToPlantingHeight[cortriGrowth$biocharLevel == "L"] ~ cortriGrowth$dateSampled[cortriGrowth$biocharLevel == "L"],
        boxfill = "#AAAAAA", add = TRUE, boxwex = boxWidth, at = c(0,unique(cortriGrowth$recodedDateForPlot[cortriGrowth$biocharLevel == "L" & !is.na(cortriGrowth$heightNormToPlantingHeight)])[-1]),
        xaxt = "n", yaxt = "n")
boxplot(cortriGrowth$heightNormToPlantingHeight[cortriGrowth$biocharLevel == "H"] ~ cortriGrowth$dateSampled[cortriGrowth$biocharLevel == "H"],
        boxfill = "#555555", add = TRUE, boxwex = boxWidth, at = c(0,unique(cortriGrowth$recodedDateForPlot[cortriGrowth$biocharLevel == "H" & !is.na(cortriGrowth$heightNormToPlantingHeight)])[-1]) + boxOffset,
        xaxt = "n", yaxt = "n")
axis(side = 1,at = axisHolder_2013_2014, las = 2, labels = c("May","Jun","Jul","Aug","Sep","Oct","Nov","May","Jun","Jul","Aug","Sep"),
     cex.axis = axisLabelSize)

abline(v = num2013Days + (daySpacer/2), col = "black", lty = 1)
abline(h = 1, col = "red", lty = 2)
# ========================

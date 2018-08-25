# This plots rarefaction curves for the rate of OTU accumulation per increase in number of sequences
# obtained for each of the biological (different DNA extractions) and technical (different PCR amplifications of the 
# same DNA extraction) replicates for the initial AMF inoculum used.

# This makes the plot in Fig S4 of the manuscript.

if(!(require("vegan"))){
    install.packages("vegan")
} 

library(vegan)

rawOTUTable_Bayles <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_AbundOTU_97PercClustered_OTUTable_allEntries.tsv",
                                 header = TRUE, sep = "\t")

rawOTUFrame_Bayles <- data.frame(rawOTUTable_Bayles[, 7:length(rawOTUTable_Bayles[1,])])

# These are IDd from the phylo
AMFOTUNumbers_Bayles <- scan("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_ListOfAMFOTUs_fromPhylogeny.txt")

# This works correctly to subset only the AMF OTU columns by moving them by name (eg X1,X12) 
OTUframeBayles_onlyAMFOTUs <- rawOTUFrame_Bayles[,paste("X",AMFOTUNumbers_Bayles,sep="")]

# 94 is Inoc rep 1 from 2013
# 95 is Inoc rep 2 from 2013
# 177 is Inoc rep 1 tech rep 1 from 2014
# 178 is Inoc rep 1 tech rep 2 from 2014
# 179 is Inoc rep 1 tech rep 3 from 2014
# 180 is Inoc rep 2 tech rep 1 from 2014
# 181 is Inoc rep 2 tech rep 2 from 2014
# 182 is Inoc rep 2 tech rep 3 from 2014

OTUframe_Bayles_onlyInoculumSamples <- OTUframeBayles_onlyAMFOTUs[c(94:95,177:182),]

minInocSeqNum <- min(rowSums(OTUframe_Bayles_onlyInoculumSamples))

# Still fairly steep OTU accumulation for most of these (not a partic. strong diff between seq runs, which is good.)
rarefactionCurve <- rarecurve(OTUframe_Bayles_onlyInoculumSamples, step = 20, sample = minInocSeqNum)

# Get the x limits
xMax <- sapply(rarefactionCurve, function(x){max(attr(x, "Subsample"))})
yMax <- sapply(rarefactionCurve, function(x){max(x)})

plot(c(1,xMax), c(1,yMax), type = "n", xlab = "Number of sequences", ylab = "Number of OTUs", 
     cex.lab = 1.3, cex.axis = 1.3)

rarefyColors <- c("#9970ab","#5aae61","#9970ab","#9970ab","#9970ab","#5aae61","#5aae61","#5aae61")

for(index in seq_along(rarefactionCurve)){
    lines(x = attr(rarefactionCurve[[index]], "Subsample"), rarefactionCurve[[index]], col = rarefyColors[index], lwd = 3)
}

legend("bottomright", lty = 1, lwd = 3, c("Replicate 1", "Replicate 2"), col = rarefyColors[1:2])

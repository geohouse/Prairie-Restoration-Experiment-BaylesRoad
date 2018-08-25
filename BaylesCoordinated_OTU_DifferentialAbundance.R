# This analyzes differential abundance of AM fungal OTUs with different experimental treatments
# (+/- AM fungi, +/- biochar, and the middle versus the margin of each plot). This gives 
# the figures that were combined to make the Fig. S3.

# Install required packages from Bioconductor if necessary.
if(!(require("phyloseq"))){
    source("https://bioconductor.org/biocLite.R")
    biocLite("phyloseq")
} 

if(!(require("DESeq2"))){
    source("https://bioconductor.org/biocLite.R")
    biocLite("DESeq2")
} 

library('phyloseq')
library('DESeq2')

# Input data is from three sources:
#   1) OTU table from which the differential abundances are calculated
#   2) Metadata about each sample (AM fungal treatment, biochar treatment, where in 
#   each plot the sample was collected [middle vs. margin])
#   3) Taxonomic attribution of each OTU

OTUTable <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_AbundOTU_97PercClustered_OTUTable_onlyAMFOTUs_onlyFieldSamples.tsv", header = T, sep = "\t")

OTUTable_metadata <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_AbundOTU_97PercClustered_OTUTable_onlyAMFOTUs_onlyFieldSamples_metadata.tsv", header = T, sep = "\t")

# These are the OTU taxonomic attributions.
AMF_OTU_attr <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_AMF_OTU_TaxonomicAttributions.csv", sep = ",", header = TRUE)

samplingYear_factor <- as.factor(OTUTable_metadata$Year)

AMFTreatment_factor <- as.factor(OTUTable_metadata$AM_fungal_inoc.)

biocharTreatment_factor <- as.factor(OTUTable_metadata$Biochar)

plotLocation_factor <- as.factor(OTUTable_metadata$Plot_Location)
plotLocation_factor_2013 <- plotLocation_factor[which(samplingYear_factor == "2013")]
plotLocation_factor_2014 <- plotLocation_factor[which(samplingYear_factor == "2014")]

block_factor <- as.factor(c(rep(seq(1,6,1), each = 12), rep(seq(1,6,1), each = 12)))

# Combined 2013 + 2014
# Added year to model to help control for changes between years.
# =========================

OTUTable_metadata$Year <- as.factor(OTUTable_metadata$Year)

OTUTable_metadata$block <- as.factor(rep(rep(seq(1,6,1), each = 12),2))

# Need to add a pseudocount to all entries in the OTU table
OTUTable <- OTUTable + 1

# Convert OTU table to phyloseq OTU table
phyloseq_OTUTable_2013_2014 <- otu_table(OTUTable, taxa_are_rows = FALSE)

# This seems to work, but View() fails on it....
phyloseq_sampleData_2013_2014 <- sample_data(OTUTable_metadata)

phyloseq_combined_2013_2014 <- phyloseq(phyloseq_OTUTable_2013_2014, phyloseq_sampleData_2013_2014)

# Now set up models for only the middle of the plots
OTUTable_center <- OTUTable[which(plotLocation_factor == "C"),]
OTUTable_metadata_center <- OTUTable_metadata[which(plotLocation_factor == "C"),]

phyloseq_OTUTable_2013_2014_center <- otu_table(OTUTable_center, taxa_are_rows = FALSE)
phyloseq_sampleData_2013_2014_center <- sample_data(OTUTable_metadata_center)

phyloseq_combined_2013_2014_center <- phyloseq(phyloseq_OTUTable_2013_2014_center, phyloseq_sampleData_2013_2014_center)

# This is a function to plot the alternating light/dark bars indicating OTUs attributed
# to different genera/families.
plotBars <- function(){
    lineThickness <- 3
    abLineColor <- "gray"
    abLineThickness <- 0.75
    segmentBarWidth <- 1.2
    # Keep at 0.3 and save as 16"W by 4" H
    barOffset <- 0.3
    
    contrastBarColor1 <- "black"
    contrastBarColor2 <- "gray"
    
    # For testing
    #yLocation <- 0
    
    # For real
    yLocation = -2.98
    
    # Updated attributions using the Mor elo outgroup 080217.
    
    # For Div (5)
    segments(x0 = barOffset, 
             y0 = yLocation, 
             x1 = barOffset + 5 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 5 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Aca. (12)
    segments(x0 = barOffset + 5 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 17 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor2, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 17 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Scu. (5)
    segments(x0 = barOffset + 17 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 22 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 22 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Cet. (2)
    segments(x0 = barOffset + 22 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 24 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor2, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 24 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Rac. (2)
    segments(x0 = barOffset + 24 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 26 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 26 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Den. (1)
    segments(x0 = barOffset + 26 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 27 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor2, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 27 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Gig. (2)
    segments(x0 = barOffset + 27 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 29 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 29 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Claroideoglomeraceae. (21)
    segments(x0 = barOffset + 29 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 50 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor2, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 50 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Cla. (7)
    segments(x0 = barOffset + 50 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 57 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 57 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Glomeraceae (47)
    segments(x0 = barOffset + 57 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 104 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor2, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 104 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Rhi. (65)
    segments(x0 = barOffset + 104 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 169 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 169 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Glomus (3)
    segments(x0 = barOffset + 169 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 172 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor2, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 172 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Fun. (6)
    segments(x0 = barOffset + 172 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 178 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 178 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Sep. (2)
    segments(x0 = barOffset + 178 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 180 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor2, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 180 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Archaeosporaceae (11)
    segments(x0 = barOffset + 180 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 191 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 191 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Archaeospora (1)
    segments(x0 = barOffset + 191 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 192 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor2, lwd = lineThickness, lend = "butt")
    abline(v = barOffset + 192 * segmentBarWidth, col = abLineColor, lty = 2, lwd = abLineThickness)
    
    # For Par. (8)
    segments(x0 = barOffset + 192 * segmentBarWidth, 
             y0 = yLocation, 
             x1 = barOffset + 200 * segmentBarWidth, 
             y1 = yLocation, 
             col = contrastBarColor1, lwd = lineThickness, lend = "butt")
}

# Fig S2C
# AMF full plot 
# ---------------
# Including block first controls for block effects, but only returns the results of the last var in the formula (either AMF or Biochar)
OTU_sampleDataForDESeq_2013_2014 <- phyloseq_to_deseq2(phyloseq_combined_2013_2014, ~ Year + block + AM_fungal_inoc.)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# in west; - val = higher in east sites), along with p-values for the log2 fold change.
DESeqOutput_2013_2014_AMF <- DESeq(OTU_sampleDataForDESeq_2013_2014, fitType = 'local')

# The first level to contrast is the numerator in the log2 ratio so values >0 are skewed towards that level, 
# while values < 0 are skewed towards the second (denominator) level given to contrast.
DESeqResults_2013_2014_AMF <- results(DESeqOutput_2013_2014_AMF, contrast = c("AM_fungal_inoc.", "Y", "N"))

# Re-order the results to match the order of the taxonomic attribution labels.
labelMatch_2013_2014_AMF <- match(paste0("X",AMF_OTU_attr$OTUNum), rownames(DESeqResults_2013_2014_AMF))

DESeqResultsLabelSorted_2013_2014_AMF <- DESeqResults_2013_2014_AMF[labelMatch_2013_2014_AMF,]

DESeqResultsPValSorted_2013_2014_AMF <- DESeqResults_2013_2014_AMF[order(DESeqResults_2013_2014_AMF$padj),]
DESeqResultsEffectSizeSorted_2013_2014_AMF <- DESeqResults_2013_2014_AMF[order(abs(DESeqResults_2013_2014_AMF$log2FoldChange), decreasing = TRUE),]

labelColor <- ifelse(DESeqResultsLabelSorted_2013_2014_AMF$padj < 0.05, 'orange', 'white')
barplot(DESeqResultsLabelSorted_2013_2014_AMF$log2FoldChange, col = labelColor,
        ylab = 'Pos (+ AMF)/ Neg (- AMF)', main = '2013_2014 AMF full plot',
        ylim = c(-3,3))

# Add the alternating light/dark bars along the bottom of the plot
plotBars()

# --------------------

# Fig S2B
# Biochar full plot 
# ---------------
# Including block first controls for block effects, but only returns the results of the last var in the formula (either AMF or Biochar)
OTU_sampleDataForDESeq_2013_2014_biochar <- phyloseq_to_deseq2(phyloseq_combined_2013_2014, ~ Year + block + Biochar)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# in west; - val = higher in east sites), along with p-values for the log2 fold change.
DESeqOutput_2013_2014_biochar <- DESeq(OTU_sampleDataForDESeq_2013_2014_biochar, fitType = 'local')

# The first level to contrast is the numerator in the log2 ratio so values >0 are skewed towards that level, 
# while values < 0 are skewed towards the second (denominator) level given to contrast.
DESeqResults_2013_2014_biochar <- results(DESeqOutput_2013_2014_biochar, contrast = c("Biochar", "H", "N"))

labelMatch_2013_2014_biochar <- match(paste0("X",AMF_OTU_attr$OTUNum), rownames(DESeqResults_2013_2014_biochar))

DESeqResultsLabelSorted_2013_2014_biochar <- DESeqResults_2013_2014_biochar[labelMatch_2013_2014_biochar,]

DESeqResultsPValSorted_2013_2014_biochar <- DESeqResults_2013_2014_biochar[order(DESeqResults_2013_2014_biochar$padj),]
DESeqResultsEffectSizeSorted_2013_2014_biochar <- DESeqResults_2013_2014_biochar[order(abs(DESeqResults_2013_2014_biochar$log2FoldChange), decreasing = TRUE),]

labelColor <- ifelse(DESeqResultsLabelSorted_2013_2014_biochar$padj < 0.05, 'orange', 'white')
barplot(DESeqResultsLabelSorted_2013_2014_biochar$log2FoldChange, col = labelColor,
        ylab = 'Pos (High)/ Neg (No)', main = '2013_2014 Biochar full plot',
        ylim = c(-3,3))
# Add the alternating light/dark bars along the bottom of the plot
plotBars()


# --------------------

# Fig S2A
# plot location 
# ---------------
# Including block first controls for block effects, but only returns the results of the last var in the formula (either AMF or Biochar)
OTU_sampleDataForDESeq_2013_2014_plotLoc <- phyloseq_to_deseq2(phyloseq_combined_2013_2014, ~ Year + block + Plot_Location)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# in west; - val = higher in east sites), along with p-values for the log2 fold change.
DESeqOutput_2013_2014_plotLoc <- DESeq(OTU_sampleDataForDESeq_2013_2014_plotLoc, fitType = 'local')

# The first level to contrast is the numerator in the log2 ratio so values >0 are skewed towards that level, 
# while values < 0 are skewed towards the second (denominator) level given to contrast.
DESeqResults_2013_2014_plotLoc <- results(DESeqOutput_2013_2014_plotLoc, contrast = c("Plot_Location", "C", "R"))

labelMatch_2013_2014_plotLoc <- match(paste0("X",AMF_OTU_attr$OTUNum), rownames(DESeqResults_2013_2014_plotLoc))

DESeqResultsLabelSorted_2013_2014_plotLoc <- DESeqResults_2013_2014_plotLoc[labelMatch_2013_2014_plotLoc,]

DESeqResultsPValSorted_2013_2014_plotLoc <- DESeqResults_2013_2014_plotLoc[order(DESeqResults_2013_2014_plotLoc$padj),]
DESeqResultsEffectSizeSorted_2013_2014_plotLoc <- DESeqResults_2013_2014_plotLoc[order(abs(DESeqResults_2013_2014_plotLoc$log2FoldChange), decreasing = TRUE),]

labelColor <- ifelse(DESeqResultsLabelSorted_2013_2014_plotLoc$padj < 0.05, 'orange', 'white')
barplot(DESeqResultsLabelSorted_2013_2014_plotLoc$log2FoldChange, col = labelColor,
        ylab = 'Pos (Center)/ Neg (Border)', main = '2013_2014 location in plot',
        ylim = c(-3,3))
# Add the alternating light/dark bars along the bottom of the plot
plotBars()

# --------------------


# Fig S2D
# AMF center of plots - no OTUs sign.
# ---------------
# Including block first controls for block effects, but only returns the results of the last var in the formula (either AMF or Biochar)
OTU_sampleDataForDESeq_2013_2014_center <- phyloseq_to_deseq2(phyloseq_combined_2013_2014_center, ~ Year + block + AM_fungal_inoc.)

# This tests the log2 fold change of the number of seqs in remnant compared to disturbed sites for each OTU (+ val = higher
# in west; - val = higher in east sites), along with p-values for the log2 fold change.
DESeqOutput_2013_2014_AMF_center <- DESeq(OTU_sampleDataForDESeq_2013_2014_center, fitType = 'local')

# The first level to contrast is the numerator in the log2 ratio so values >0 are skewed towards that level, 
# while values < 0 are skewed towards the second (denominator) level given to contrast.
DESeqResults_2013_2014_AMF_center <- results(DESeqOutput_2013_2014_AMF_center, contrast = c("AM_fungal_inoc.", "Y", "N"))
DESeqResultsPValSorted_2013_2014_AMF_center <- DESeqResults_2013_2014_AMF_center[order(DESeqResults_2013_2014_AMF_center$padj),]

labelMatch_2013_2014_AMF_center <- match(paste0("X",AMF_OTU_attr$OTUNum), rownames(DESeqResults_2013_2014_AMF_center))

DESeqResultsLabelSorted_2013_2014_AMF_center <- DESeqResults_2013_2014_AMF_center[labelMatch_2013_2014_AMF_center,]
DESeqResultsEffectSizeSorted_2013_2014_AMF_center <- DESeqResults_2013_2014_AMF_center[order(abs(DESeqResults_2013_2014_AMF_center$log2FoldChange), decreasing = TRUE),]

labelColor <- ifelse(DESeqResultsLabelSorted_2013_2014_AMF_center$padj < 0.05, 'orange', 'white')
barplot(DESeqResultsLabelSorted_2013_2014_AMF_center$log2FoldChange, col = labelColor,
        ylab = 'Pos (+ AMF)/ Neg (-AMF)', main = '2013_2014 AMF center of plots',
        ylim = c(-3,3))

# Add the alternating light/dark bars along the bottom of the plot
plotBars()




# --------------------


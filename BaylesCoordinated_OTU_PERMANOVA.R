# This script tests for direct differences in the AM fungal community (as represented
# by OTUs) using PERMANOVA. This runs the test reported in the manuscript, where AM fungal
# inoculation in the middle of plots during the second growing season is nearly a 
# significant predictor of AM fungal community composition during the
# second growing season. All other tests run indicated that the treatments used were not
# significant predictors of the AM fungal community composition.

if(!(require("vegan"))){
    install.packages("vegan")
} 

library("vegan")

# This requires the OTU table and the metadata for the OTU table
OTUTable <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_AbundOTU_97PercClustered_OTUTable_onlyAMFOTUs_onlyFieldSamples.tsv", header = T, sep = "\t")

OTUTable_metadata <- read.table("~/Box Sync/R_code/Bayles_biochar/Bayles_forGitHub/Bayles_AbundOTU_97PercClustered_OTUTable_onlyAMFOTUs_onlyFieldSamples_metadata.tsv", header = T, sep = "\t")

samplingYear_factor <- as.factor(OTUTable_metadata$Year)

OTUTable_2014 <- OTUTable[which(samplingYear_factor == "2014"), ]

log10_numSeqs_2014 <- log10(rowSums(OTUTable[which(samplingYear_factor == "2014"), ]))

# Only 2014 data
AMFTreatment_factor_2014 <- as.factor(OTUTable_metadata$AM_fungal_inoc.[which(samplingYear_factor == "2014")])
biocharTreatment_factor_2014 <- as.factor(OTUTable_metadata$Biochar[which(samplingYear_factor == "2014")])
plotLocation_factor_2014 <- as.factor(OTUTable_metadata$Plot_Location[which(samplingYear_factor == "2014")])

# Only 2014 Middle of plots
# ------------------------
OTUTable_2014_center <- OTUTable_2014[which(plotLocation_factor_2014 == "C"),]

AMFTreatment_factor_2014_center <- AMFTreatment_factor_2014[which(plotLocation_factor_2014 == "C")]
biocharTreatment_factor_2014_center <- biocharTreatment_factor_2014[which(plotLocation_factor_2014 == "C")]
block_factor_2014_center <- as.factor(c(rep(seq(1,6,1), each = 6)))
log10_numSeqs_2014_center <- log10_numSeqs_2014[which(plotLocation_factor_2014 == "C")]


# Use log10(seqNumber) as a covariate and stratify by block.
# Set the seed to ensure identical permutation results each time this is run.
set.seed(434)
bayles_permanova_2014_main_center <- adonis2(formula = OTUTable_2014_center ~  AMFTreatment_factor_2014_center + 
                                                 biocharTreatment_factor_2014_center + log10_numSeqs_2014_center, strata = block_factor_2014_center,
                                             permutations = 5000, method = "morisita", binary = FALSE, by = "margin")

print(bayles_permanova_2014_main_center)
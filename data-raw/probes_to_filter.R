# M. price 2013 list of cross-hybridizing probes -------------------------

# Load annotation
UBC_anno <- read.table('Z:/Victor/Documents/UBC Annotation Combined.txt', fill =
                         TRUE, header = TRUE)

# 641 XY cross hybridizing probes
UBC_annoXY <- as.character(UBC_anno[grep("YES",UBC_anno$XY_Hits), 'IlmnID'])
# 2049 Autosomal probes
UBC_annoAuto <- as.character(UBC_anno[grep("YES",UBC_anno$Autosomal_Hits),
                                      'IlmnID'])

CH_mprice <- union(UBC_annoAuto, UBC_annoXY) # AND/OR

# Overlapping EPIC and 450k probes ------------------------------------
library(minfiData)
library(minfiDataEPIC)

data(RGsetEx)
cpg_450k <- rownames(getBeta(RGsetEx))
snp_450k <- rownames(getSnpBeta(RGsetEx))

data(RGsetEPIC)
cpg_EPIC <- rownames(getBeta(RGsetEPIC))
snp_EPIC <- rownames(getSnpBeta(RGsetEPIC))

overlap_450k_EPIC <- intersect(c(cpg_450k, snp_450k), c(cpg_EPIC, snp_EPIC))

# Invariable probes --------------------------------------------------------
library(RCurl) 
base <- paste("https://raw.githubusercontent.com/redgar598/",
        "Tissue_Invariable_450K_CpGs/master/", sep = "")
csv <- c("Invariant_Buccal_CpGs.csv", 
         "Invariant_Blood_CpGs.csv", 
         "Invariant_Placenta_CpGs.csv")
invar_Bu <- read.csv(text = getURL(paste(base, csv[1], sep = "")))$CpG
invar_Bl <- read.csv(text = getURL(paste(base, csv[2], sep = "")))$CpG
invar_Pl <- read.csv(text = getURL(paste(base, csv[3], sep = "")))$CpG



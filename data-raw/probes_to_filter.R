# M. price 2013 list of cross-hybridizing probes -------------------------
# Load annotation
library(dplyr)
library(GEOquery)
base <-'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL16304&format=file'
end <- '&file=GPL16304%5FGene%5Ffeatures%5FPlatformTable%2Etxt%2Egz'

mprice <- as_tibble(getGEO('GPL16304')@dataTable@table)
mprice_XY <- mprice %>% filter(XY_Hits == 'XY_YES') %>% pull(ID)
mprice_AU <- mprice %>% filter(Autosomal_Hits == 'A_YES') %>% pull(ID)

mprice_CH <- union(mprice_XY, mprice_AU)

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

library(devtools)
use_data(mprice_CH, overlap_450k_EPIC, invar_Bu, invar_Bl, invar_Pl, 
         overwrite = T)
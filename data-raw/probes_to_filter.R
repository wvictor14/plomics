# M. price 2013 list of cross-hybridizing probes -------------------------

# Load annotation
UBC_anno <- read.table('Z:/Victor/Documents/UBC Annotation Combined.txt', fill = TRUE, header = TRUE)

# 641 XY cross hybridizing probes
UBC_annoXY <- as.character(UBC_anno[grep("YES",UBC_anno$XY_Hits), 'IlmnID'])
# 2049 Autosomal probes
UBC_annoAuto <- as.character(UBC_anno[grep("YES",UBC_anno$Autosomal_Hits),'IlmnID'])

CH_mprice <- union(UBC_annoAuto, UBC_annoXY) # AND/OR

#' List of XY- and autosomal cross-hybridizing probes from m.price et al. (2013)
#'
#'
#' A character vector of 41937 CpG names that were identified to be able to cross
#' hybridize to other locations in the genome (both XY and autosomal).
#'
#'@format A character vector of CpG names with length 41937
#'@source Chen Y, Lemire M, Choufani S, et al. Discovery of cross-reactive 
#'probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 
#'microarray. Epigenetics. 2013;8(2):203-209. doi:10.4161/epi.23470.
#' 
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304
"mprice_CH"

#' List of cross-hybridizing probes from Chen et al. (2013)
#'
#'
#' A character vector of 29233 CpG names that were identified to be able to 
#' cross hybridize to other locations in the genome.
#'
#'@format A character vector of CpG names with length 29233
#'@source Price, E. M., Cotton, A. M., Lam, L. L., Farré, P., Emberly, E.,
#'Brown, C. J., … Kobor, M. S. (2013). Additional annotation enhances potential
#'for biologically-relevant analysis of the Illumina Infinium
#'HumanMethylation450 BeadChip array. Epigenetics & Chromatin, 6, 4.
#' http://doi.org/10.1186/1756-8935-6-4
#' 
#'http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/'
#'file: '48639-non-specific-probes-Illumina450k.xlsx'
"chen_CH"


#' Probes that are on both the 450K and EPIC arrays
#'
#' A character vector of 453152 probes names that are on both the 450k and EPIC 
#'Illummina DNA methylation microarrays. All of the snp probes on the EPIC (59)
#'are included - only 6 are missing from the 450k array. These names were taken
#'from minfiData and minfiDataEPIC, because the illumina annotations do not have
#'snp probe information.
#'
#'@format A character vector of probe names with length 453152
#'@source minfiData 0.24.0, minfiDataEPIC 1.4.0
"overlap_450k_EPIC"

#'Probes identified as invariable in indepedent buccal datasets
#'
#'A character vector of 120009 probe names identified as nonvariable between
#'multiple BUCCAL datasets.
#' 
#'@format A character vector of probe names with length 120009
#'@source Edgar R.D., Jones M.J., Robinson W.P., Kobor M.S. (2017) An 
#'empirically driven data reduction method on the human 450K methylation array 
#'to remove tissue specific non-variable CpGs. Clin. Epigenet ., 9, 11.
#'https://github.com/redgar598/Tissue_Nonvariable_450K_CpGs
"invar_Bu"

#'Probes identified as invariable in indepedent blood datasets
#' 
#'A character vector of 114204 probe names identified as nonvariable between
#'multiple BLOOD datasets.
#' 
#'@format A character vector of probe names with length 114204
#'@source Edgar R.D., Jones M.J., Robinson W.P., Kobor M.S. (2017) An 
#'empirically driven data reduction method on the human 450K methylation array 
#'to remove tissue specific non-variable CpGs. Clin. Epigenet ., 9, 11.
#'https://github.com/redgar598/Tissue_Nonvariable_450K_CpGs
"invar_Bl"

#'Probes identified as invariable in indepedent placenta datasets
#' 
#'A character vector of 101367 probe names identified as nonvariable between
#'multiple PLACENTA datasets.
#' 
#'@format A character vector of probe names with length 101367
#'@source Edgar R.D., Jones M.J., Robinson W.P., Kobor M.S. (2017) An 
#'empirically driven data reduction method on the human 450K methylation array 
#'to remove tissue specific non-variable CpGs. Clin. Epigenet ., 9, 11.
#'https://github.com/redgar598/Tissue_Nonvariable_450K_CpGs
"invar_Pl"

#' Analyze detection p-values
#'
#' \code{analyzeDetP} Applies given thresholds to detection p value matrix to 
#' identify poor quality probes and samples.
#'
#' @param detP A p by n dataframe matrix of detection p values. For example from
#' minfi::detectionP(rgset).
#' @param pDat A n by m dataframe/tibble of metadata.
#' @param cpg_threshold A numeric value x, where 0 < x < 1, that specifies the 
#' detection p value threshold for a poor quality probe.
#' @return List of length two:
#' \itemize{
#'   \item fp: p by n dataframe with elements where 1's are placed when an
#'   observation failed cpg_threshold, and 0's otherwise.
#'   \item pDat: n by m + 1 dataframe containing a column ndetp corresponding
#'   to the number of probes that failed in sample i, where i = 1, 2, ..., n.
#' }
#' @examples 
#' library(minfiData)
#' library(dplyr)
#' data(RGsetEx)
#' 
#' x <- analyzeDetP(detP = detectionP(RGsetEx), 
#'                  pDat = as_tibble(pData(RGsetEx), 
#'                  cpg_threshold = 0.01))
analyzeDetP <- function(detP, pDat, cpg_threshold = 0.01) {
  if (ncol(detP) != nrow(pDat)){
    stop('Number of columns in detP does not match the number of rows in pDat')
  }
  
  pDat <- pDat %>% mutate(ndetp = colSums(detP > 0.01, na.rm = F))
  
  if (sum(pDat$ndetp > 0.05*nrow(RGsetEx))>0) {
    x <- pDat %>% filter(ndetp > 0.05*nrow(RGsetEx)) %>% select(Sample_Name)
    print(paste('The following samples had > 5% of their probes failing a 
                detection p value of 0.05%', x$Sample_Name))  
  }                     
  
  fp <- detP > 0.01
  return(list(as_tibble(fp), as_tibble(pDat)))
}
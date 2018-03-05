#' Preprocess methylation array data
#'
#' \code{preprocess} returns a list of length two, containing (1) filtered,
#' normalized dataframe and (2) pData with QC information
#'
#' @param mset A minfi::methylset object.
#' @param fp An p by n matrix/dataframe, where n = the number of samples and p =
#'  number of features. This should be the same dimensions as mset, unless snp
#'  probes are included in mset.
#' @param threshold_p A numeric value between 0 and 1 representing the fraction 
#' of observations required for a probe to be removed.
#' @param threshold_s A numeric value between 0 and 1 representing the fraction
#' of observations required for a sample to be removed.
#' @param ch A character vector specificying the set of cross hybridizing probes 
#' to filter out. Options: 'mprice', 'chen', 'both'.
#' @param threshold_cor A numeric value between 0 and 1 specifying the minimum
#' mean interarray correlation to keep a sample. If not supplied, then 
#' correlations are not calculated.
#' @param invariable_probes A character vector of specifying the set of
#' invariable probes to filter out. Options: 'Blood', 'Buccal', 'Placenta'.
#' @param overlapping A logical vector indicating if 450k/EPIC overlapping 
#' probes should be kept.
#' @param wrongsex A character vector of samplenames to filter out.
#' @param seed Sets a seed for normalization.
#' @return List of length two:
#' \itemize{
#'   \item data - A filtered and normalized dataframe
#'   \item pDat - A dataframe containing updated phenotype data information
#'   \item num_fp - number of probes filtered based on quality of signal
#'   \item num_ch - number of cross hybridizing probes filtered
#'   \item num_iv - number of invariable probes filtered
#' }
#' @details This function takes in a methylset object, and does the following:
#' (1) Filters poor quality probes based on fp, (2) normalizes beta values with
#' wateRmelon::BMIQ, (3) filters overlapping 450k/EPIC probes and invariable
#' probes. FP needs to be precomputed. Any of these probes can be omitted.
#' 
#' The default for probe filtering is 0.01, and the default for sample filtering
#' is 0.05.
#' 
#' The reason of filtering overlapping / invariable probes is that normalization
#' can benefit from these probes as extra observations, whereas poor quality 
#' probes are removed prior to normalization.
#' 
preprocess <- function(mset, fp = NULL, threshold_p = 0.01, threshold_s = 0.05, 
                       ch = NULL, invariable_probes = NULL, threshold_cor = NULL,
                      overlapping = F, wrongsex = NULL, seed = 1){
  pData(mset)$rownames <- rownames(pData(mset))
  pDat <- as_tibble(pData(mset))

  #checks
  #if (all(rownames(fp) != rownames(mset))){
  #  stop('row names of failed probes matrix needs to match row names of mset')
  #}

  if (length(wrongsex) > 0){
    if (!all(wrongsex %in% colnames(mset))){
      warning(paste(setdiff(wrongsex, colnames(mset)), 'are not in the data'))
    }
  }

  # filter probes on fp and ch
  if (!is.null(fp)) {
    # calculate failed probes / samples
    bp <- rownames(fp)[rowSums(fp) > threshold_p*ncol(fp)] # failed probes
    bs <- colnames(fp)[colSums(fp) > threshold_s*nrow(fp)] # failed samples

    # remove badprobes
    print(paste('Removing', length(bp), 'poor quality probes.'))
    mset <- mset[setdiff(rownames(mset), bp),]
    
    # remove samples based on num of detect_p and bead count
    if (length(bs) > 0){
      print(paste('Removing', length(bs), 'samples that had > 5% poor quality probes.'))
      mset <- mset[,setdiff(colnames(mset), bs)]
      pDat <- as_tibble(pData(mset))
    }
  }
  
  # remove crosshybridizing
  if (!is.null(ch)){
    if (ch == 'both'){
      cp <- union(mprice_CH, chen_CH)
    }
    if (ch == 'mprice'){
      cp <- mprice_CH
    }
    if (ch == 'chen'){
      cp <- chen_CH
    } 
    cp <- intersect(cp, rownames(mset))
    print(paste('Removing', length(cp),'cross hybridizing probes.'))
    mset <- mset[setdiff(rownames(mset), cp),]
  } else {
      stop('ch must be either "both", "chen", "mprice", or NULL.')
  }
  
#TODO: for custom supplied ch
#    if ((length(ch) > 1) && (intersect(ch, rownames(mset)) > 1)){
#      cp <- intersect(ch, rownames(mset))
#      print(paste(length(cp), ' ch probes found in mset.', sep = ''))
#    }

  # mean interarray correlation
  if (!is.null(threshold_cor)) {
    c <- cor(getBeta(mset), use = 'pairwise.complete.obs')
    sample_c <- apply(c, 1, mean) # Calculate mean
    pDat$meanSScor <- sample_c[match(names(sample_c), pDat$rownames)] 
  
    # remove samples based on mean interarray cor
    if (sum(pDat$meanSScor > threshold_cor) > 0){
      bs_c <- pDat %>% filter(threshold_cor < threshold_cor) %>% pull(rownames)
      print(paste('Removing', length(bs_c), 
                  'samples that had < 0.95 mean interarray correlation.'))
      pData(mset) <- DataFrame(pDat)
      mset <- mset[,setdiff(colnames(mset), bs_c)]
      pDat <- as_tibble(pData(mset))
    }
  }
  
  # remove incorrect sex samples
  if ((length(wrongsex) > 0) & (all(wrongsex %in% colnames(mset)))) {
    print(paste('Removing', length(wrongsex), 'incorrectly labelled samples, based on sex inference.'))
    pData(mset) <- DataFrame(pDat)
    mset <- mset[,setdiff(colnames(mset), wrongsex)]
    pDat <- as_tibble(pData(mset))
  }

  # Normalization
  print('Normalizing data...')
  set.seed(seed)
  BMIQ <- as.data.frame(BMIQ(mset), nfit = 50000)

  # remove non-overlapping probes
  if (overlapping == T) {
    print(paste('Keeping', length(intersect(rownames(BMIQ), overlap_450k_EPIC)),
                'probes that are not present in 450k and EPIC'))
    BMIQ <- BMIQ[intersect(rownames(BMIQ), overlap_450k_EPIC),]
  }

  # remove invariable probes
  if (!is.null(invariable_probes)) {
    if (invariable_probes != 'Blood'){
      if (invariable_probes != 'Buccal') {
        if (invariable_probes != 'Placenta')
          stop(paste('Only blood, buccal, and placenta invariable probe', 
                     'references are available.'))
      }
    }
    # calculate range per probe
    print('Calculating your dataset-specific invariable probes...')
    Variation<-function(x) {quantile(x, c(0.9),na.rm=T)[[1]] -
        quantile(x, c(0.1), na.rm=T)[[1]]}
    
    rr <- lapply(1:nrow(BMIQ), function(x) Variation(BMIQ[x,])) 
    rr <- unlist(rr)
    
    # get cpg names of probes with <0.05 range
    ip <-rownames(BMIQ)[which(rr<0.05)] 
    print(paste('There are', length(ip), 'invariable probes in your data.'))
    
    if (invariable_probes == 'Blood') {
      print('Using blood reference set...')
      iv <- intersect(ip, invar_Bl)
    }
    if (invariable_probes == 'Buccal'){
      print('Using buccal reference set...')
      iv <- intersect(ip, invar_Bu)
    }
    if (invariable_probes == 'Placenta'){
      print('Using placenta reference set...')
      iv <- intersect(ip, invar_Pl)
    }
    
    print(paste(length(iv),
                ' of your dataset-specific probes overlap the reference, and',
                ' will be removed.', sep = ''))

    BMIQ <- BMIQ[setdiff(rownames(BMIQ), iv),]
  } else {
    if (!is.null(invariable_probes)) {
      warning('No invariable probes filtered. Check format.')
    }
  }
  
  list(data = as_tibble(BMIQ, rownames = 'rownames'), pData = pDat,
       num_fp = ifelse(!is.null(fp), length(bp), 0), 
       num_ch = ifelse(!is.null(ch), length(cp), 0),
       num_iv = ifelse(!is.null(invariable_probes), length(iv), 0))
}

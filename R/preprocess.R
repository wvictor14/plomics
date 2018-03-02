#' Preprocess methylation array data
#'
#' \code{preprocess} returns a list of length two, containing (1) filtered,
#' normalized dataframe and (2) pData with QC information
#'
#' @param mset A minfi::methylset object.
#' @param fp An p by n matrix/dataframe, where n = the number of samples and p =
#'  number of features. This should be the same dimensions as mset, unless snp
#'  probes are included in mset.
#' @param ch A character vector of cross hybridizing probes to be filtered out.
#' @param invariable_probes A character vector of reference invariable probes to
#' filter against.
#' @param overlapping A character vector of overlapping probes on 450k/EPIC
#' array to keep.
#' @param wrongsex A character vector of samplenames to filter out.
#' @param seed Sets a seed for normalization.
#' @return List of length two:
#' \itemize{
#'   \item data - A filtered and normalized dataframe
#'   \item pDat - A dataframe containing updated phenotype data information
#' }


preprocess <- function(mset, fp = NULL, ch = NULL, invariable_probes = NULL,
                       pkeep = NULL, overlapping = NULL, wrongsex, seed = 1){
  # This function filters an mset and BMIQ normalizes any mset.
  # fp must be same dimension as mset, and must contain 1 if a observation is failed, or 0 if not
  # ch is a character vector of cross hybridizing probe names to filter out
  # overlapping is a character vector of overlapping probes on both 450k/EPIC

  # Note that I filter non-overlapping probes after normalization because these
  # probes are still useful for QC.

  # Invariable probes are also calculated and filtered after normalization.
  # invariable_probes is a character vector of probe names that are invariable
  # in an independent dataset (from R. Edgar)

  # wrongsex is a character vector of sample names that are to be filtered out based on sex inference
  # if pkeep is specified then fp, ch, and invariable_probes are overridden.
  # pkeep is a list of probes to keep (everything else will be filtered out)

  # load libraries
  require(wateRmelon)
  require(minfi)
  require(IlluminaHumanMethylation450kmanifest)

  pData(mset)$rownames <- rownames(pData(mset))
  pDat <- as.data.frame(pData(mset))

  #checks
  #if (all(rownames(fp) != rownames(mset))){
  #  stop('row names of failed probes matrix needs to match row names of mset')
  #}
  if (is.null(pkeep)) {
    if (all(colnames(fp) != colnames(mset)) && !is.null(fp)){
      stop('column names of failed probes matrix needs to match column names of mset')
    }
    if (!all(ch %in% rownames(mset))){
      warning(paste('only', sum(ch %in% rownames(mset)),
                    'out of', length(ch), 'ch probes specified are in the data'))
    }
  } else {
    if (!is.character(pkeep)){
      stop('pkeep is not a character vector')
    }
    print('Filtering probes not in pkeep...')
  }

  if (length(wrongsex) > 0){
    if (!all(wrongsex %in% colnames(mset))){
      warning(paste(setdiff(wrongsex, colnames(mset)), 'are not in the data'))
    }
  }

  # filter probes on fp and ch
  if (is.null(pkeep)) {
    # calculate failed probes / samples
    bp <- rownames(fp)[rowSums(fp) > 0.01*ncol(fp)] # vector of probes names that had > 0.01 failed observations
    bs <- colnames(fp)[colSums(fp) > 0.05*nrow(fp)] # vector of samples had > 0.05 failed samples

    # remove badprobes
    print(paste('Removing', length(bp), 'probes that failed detection p, or had a bead count < 3.'))
    mset <- mset[setdiff(rownames(mset), bp),]

    # remove crosshybridizing
    print(paste('Removing', sum(ch %in% rownames(mset)), 'cross hybridizing probes.'))
    mset <- mset[setdiff(rownames(mset), ch),]

    # remove samples based on num of detect_p and bead count
    if (length(bs) > 0){
      print(paste('Removing', length(bs), 'samples that had > 5% poor quality probes.'))
      mset <- mset[,setdiff(colnames(mset), bs)]
      pDat <- as.data.frame(pData(mset))
    } else {
      print('No samples had > 5% poor quality probes.')
    }
  }

  #filter probes on pkeep
  if (!is.null(pkeep)){
    mset <- mset[pkeep, ]
  }

  # mean interarray correlation
  c <- cor(getBeta(mset), use = 'pairwise.complete.obs')
  sample_c <- apply(c, 1, mean) # Calculate mean
  pDat$meanSScor <- sample_c[match(names(sample_c), rownames(pDat))] # put into pDat
  bs_c <- rownames(pDat)[which(pDat$meanSScor < 0.95)] # extra samples with meanSS<0.95

  # remove samples based on mean interarray cor
  if (length(bs_c) > 0){
    print(paste('Removing', length(bs_c), 'samples that had < 0.95 mean interarray correlation.'))
    pData(mset) <- DataFrame(pDat)
    mset <- mset[,setdiff(colnames(mset), bs_c)]
    pDat <- as.data.frame(pData(mset))
  } else {
    print('No samples had < 0.95 mean interarray correlation.')
  }

  # remove incorrect sex samples
  if ((length(wrongsex) > 0) & (all(wrongsex %in% colnames(mset)))) {
    print(paste('Removing', length(wrongsex), 'incorrectly labelled samples, based on sex inference.'))
    pData(mset) <- DataFrame(pDat)
    mset <- mset[,setdiff(colnames(mset), wrongsex)]
    pDat <- as.data.frame(pData(mset))
  }

  # Normalization
  print('Normalizing data...')
  set.seed(seed)
  BMIQ <- as.data.frame(BMIQ(mset), nfit = 10000)

  # remove non-overlapping probes
  if (!is.null(overlapping)) {
    print(paste('Filtering', length(intersect(rownames(BMIQ), overlapping)),
                'probes that are not present in 450k and EPIC'))
    BMIQ <- BMIQ[intersect(rownames(BMIQ), overlapping),]
  }

  # remove invariable probes
  if ((!is.null(invariable_probes)) && is.null(pkeep)) {
    print(paste('There are', sum(invariable_probes %in%
                rownames(mset)), 'out of', length(invariable_probes),
                'pre-determined invariable probes in the dataset.'))
    print('Calculating your dataset-specific invariable probes.')
    Variation<-function(x) {quantile(x, c(0.9),na.rm=T)[[1]] -
        quantile(x, c(0.1), na.rm=T)[[1]]}
    rr <- lapply(1:nrow(BMIQ), function(x) Variation(BMIQ[x,])) # calculate range per probe
    rr <- unlist(rr)
    ip <-rownames(BMIQ)[which(rr<0.05)] # get cpg names of probes with <0.05 range

    print(paste('There are', length(ip), 'invariable probes in your data.'))
    print(paste(length(intersect(ip, invariable_probes)),
                'of these overlap with the independent list and will be removed.'))
    invariable <- intersect(ip, invariable_probes)

    BMIQ <- BMIQ[setdiff(rownames(BMIQ), invariable),]
  } else {
    print('No invariable probes supplied.')
  }

  list(data = BMIQ, pData = pDat) # return process matrix and pData
}

#' The cumulative amount of variance explained in each PC while adjusting for
#' other meta variables.
#'
#' \code{cumulative} progressively adds more and more PCs into a linear model,
#' while adjusting for other meta variables in \code{pDat}.
#' @param pc_matrix A rotated p x n dataframe, where p are the PC scores
#' (rotated values) and n are samples.
#' @param pDat A n x m dataframe, containing both the metavariables to be
#' adjusted and the response variable
#' @param response The metavariable to be used as the dependent variable.
#' @return returns a list of Rsquared values from each linear model
#' @details Creates progressively longer glm models of the form, while
#' adjusting for other columns in pDat, such as sex:
#' \itemize{
#'   \item response ~ pc1 + [Sex]
#'   \item response ~ pc1 + pc2 + [Sex]
#'   \item response ~ pc1 + pc2 + ...  + [Sex]
#' }
cumulative <- function(pc_matrix, pDat, response){


  require(rms)
  # Format checking...
  if(nrow(pc_matrix) != nrow(pDat))
    stop('Number of rows (samples) in rotated matrix must equal number of rows
         (samples) in metadata')

  # make sure all sample labels are the same
  stopifnot(rownames(pc_matrix) %in% rownames(pDat))
  stopifnot(rownames(pDat) %in% rownames(pc_matrix))
  stopifnot(!is.null(colnames(pc_matrix))) # colnames of scores need to be present

  data <- cbind(pDat[rownames(pc_matrix),,drop=F], pc_matrix)

  y <- colnames(pc_matrix) # take PC scores names
  x <- colnames(pDat)   # take out meta categorical labels

  results <- list()     # initialize empty vectors to fill in
  z <- colnames(pDat)[which(colnames(pDat) != response)]

  for (i in y){ #for every PC
    z <- c(z, i)
    if (length(z) != 1)
      z <- paste(z, collapse = '+')
    form <- paste(response, '~', z)
    fit <- lrm(as.formula(form), data, maxit = 25)
    results[form] <- fit$stats['R2']
  }
  return(results)
}

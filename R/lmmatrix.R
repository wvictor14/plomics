#' Creates a matrix of p values or rsquared between a matrix of independent
#' variables and a dataframe of dependent variables
#'
#' \code{heatPCA} tests association between independent variables and dependent
#' variables through simple linear regression (1 predictor variable)
#'
#' @param dep A d x n matrix of dependent variables, where the variables are
#' arranged in rows, and columns for samples or observations.
#' @param ind A i x n dataframe of independent variables, in same format as dep.
#' @param metric Can be 'Rsquared' or 'Pvalue'. Determines contents of output
#' matrix
#'
#' @details Each independent variable is tested for their association with each
#' dependent variable in simple linear regression, and a pvalue or rsquared is
#' extracted and returned as a matrix.
#' @return A i x d matrix of pvalues or rsquared values.
#'
#' @examples
#' ## to calculate PC association with covariates
#'
#' library(minfiData)
#' data(RGsetEx)
#' betas <- na.omit(getBeta(RGsetEx))
#' pDat <- as.data.frame(pData(RGsetEx))
#' pc_obj <- prcomp(t(betas))
#' pc_matrix <- t(pc_obj$x)
#' cov <- t(pDat[,c('Sample_Group', 'age', 'sex')])
#'
#' lmmatrix(dep = pc_matrix, ind = cov)
#' lmmatrix(dep = pc_matrix, ind = cov, metric = 'Pvalue')

lmmatrix <- function(dep, ind, metric = 'Rsquared'){
  # define output matrix
  n_dep <- nrow(dep)
  n_ind <- nrow(ind)

  out <- matrix(NA, ncol=n_dep, nrow=n_ind)

  # run linear models
  for(i in 1:n_ind) {
    for(j in 1:n_dep) {
      fit <- summary(lm(dep[j,] ~ ind[i,], na.action=na.omit))
      if (metric == 'Pvalue'){
        out[i,j] <- fit$coefficients[2,4]
      }
      if (metric == 'Rsquared'){
        out[i,j] <- fit$r.squared
      }
    }
  }
  rownames(out) <- rownames(ind)
  colnames(out) <- rownames(dep)

  out
}

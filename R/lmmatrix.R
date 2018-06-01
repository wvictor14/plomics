#' Creates a matrix of p values or rsquared between a matrix of independent
#' variables and a dataframe of dependent variables
#'
#' \code{heatPCA} tests association between independent variables and dependent
#' variables through simple linear regression (1 predictor variable)
#'
#' @param dep A p x n matrix of dependent variables, where p are the variables
#' and n are observations or samples.
#' @param ind A n x m dataframe, where m columns are the number of covariates
#' to be tested for association with the dependent variables
#' @param metric Can be 'Rsquared' or 'Pvalue'. Determines contents of output
#' matrix
#'
#' @details Each independent variable is tested for their association with each
#' dependent variable in simple linear regression, and a pvalue or rsquared is
#' extracted and returned as a matrix.
#' @return A p x m matrix of pvalues or rsquared values.
#'
#' @examples
#' ## to calculate PC association with covariates
#'
#' library(minfiData)
#' data(RGsetEx)
#' betas <- na.omit(getBeta(RGsetEx))
#' pDat <- pData(RGsetEx)
#' pc_obj <- prcomp(t(betas))
#' pc_matrix <- pc_obj$x
#'
#' lmmatrix(dep = pc_matrix, ind = pDat[,c('Sample_Group', 'age', 'sex')])
#' lmmatrix(dep = pc_matrix, ind = pDat[,c('Sample_Group', 'age', 'sex')],
#'          metric = 'Pvalue')

lmmatrix <- function(dep, ind, metric = 'Rsquared'){
  # define output matrix
  n_dep <- ncol(dep)
  n_ind <- ncol(ind)
  out <- matrix(NA,ncol=n_ind, nrow=n_dep)

  # run linear models
  for(j in 1:n_ind)
  {
    for(i in 1:n_dep)
    {
      fit <- summary(lm(dep[,i] ~ ind[,j], na.action=na.omit))

      if (metric == 'Pvalue'){
        out[i,j] <- fit$coefficients[2,4]
      }
      if (metric == 'Rsquared'){
        out[i,j] <- fit$r.squared
      }
    }
  }
  rownames(out) <- colnames(dep)
  colnames(out) <- names(ind)

  out
}

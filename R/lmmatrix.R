#' Tests associations between multiple independent and dependent variables
#'
#' \code{lmmatrix} tests association between independent variables and dependent
#' variables through simple linear regression (with 1 predictor variable)
#'
#' @param dep n x d matrix or dataframe of dependent variables, where the
#' variables are arranged in columns, and samples or observations in rows.
#' @param ind n x i matrix or dataframe of independent variables, where the
#' variables are arranged in columns, and samples or observations in rows.
#' @param metric Determines what is returned from the linear models. 'Rsquared' 
#' returns the proportion of variance explained, and 'Pvalue' returns the 
#' p-value.
#' @details Each independent variable is tested for their association with each
#' dependent variable in simple linear regression (ind ~ dep), and a pvalue or
#' rsquared is extracted and returned as a matrix. broom::glance() is used to 
#' calculate p value for the whole model.
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
#' pc_matrix <- pc_obj$x
#' cov <- pDat[,c('Sample_Group', 'age', 'sex')]
#'
#' lmmatrix(dep = pc_matrix, ind = cov)
#' lmmatrix(dep = pc_matrix, ind = cov, metric = 'Pvalue')
#' @export

lmmatrix <- function(dep, ind, metric = 'Rsquared'){
  # define output matrix
  n_dep <- ncol(dep)
  n_ind <- ncol(ind)
  out <- matrix(NA, ncol=n_dep, nrow=n_ind)

  # run linear models
  for(i in 1:n_ind) {
    for(j in 1:n_dep) {
      fit <- stats::lm(dep[,j,drop=TRUE] ~ ind[,i,drop=TRUE], na.action=stats::na.omit)
      if (metric == 'Pvalue'){
        out[i,j] <- broom::glance(fit)$p.value
      }
      if (metric == 'Rsquared'){
        out[i,j] <- summary(fit)$r.squared
      }
    }
  }
  rownames(out) <- colnames(ind)
  colnames(out) <- colnames(dep)

  out
}

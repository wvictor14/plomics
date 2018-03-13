#' Create a heatmap of association between PCs and clinical variables
#'
#' \code{heatPCA} tests association between principal components calculated on
#' gene expression or DNAm data and clinical variables, the plots these on a
#' heatmap.
#'
#' @param pc_matrix A rotated p x n dataframe, where p are the PC scores
#' (rotated values) and n are samples.
#' @param pDat A n x m dataframe, where m is the number of clinical / meta
#' variables to be tested for association with the PCs.
#' @param metric Defaults to 'Rsquared'. 'P-value' not implemented yet.
#' @param returnMat A logical indicating if a dataframe of association values
#' should be returned
#' @return A dataframe of association values between PCs and each meta variable
#' @details Each PC supplied (I recommend only doing the top 10-20 for
#' readability) is tested for their correlation / association with each PC by
#' computing a linear model between each PC (dependent variable) and each
#' meta variable (dependent variable). The results are then plotted using
#' pheatmap::pheatmap and RColorBrewer.
#' ##TODO:: implement % variance explained.
#' ##TODO:: implement return p-values
heatPCA <- function(pc_matrix, pDat, metric = 'Rsquared', labels, returnMat = F){

  require(RColorBrewer) # for colors
  require(pheatmap) # for heatmap
  require(reshape2) # for MELT and DCAST

  # Format checking...
  if(nrow(pc_matrix) != nrow(pDat))
    stop('Number of rows (samples) in rotated matrix must equal number of rows
         (samples) in metadata')

  # make sure all sample labels are the same
  stopifnot(!is.null(colnames(pc_matrix))) #colnames of scores must be present

  # make sure metrics are Pval or Rsqaured
  if ((metric != 'Rsquared') & (metric != 'Pval'))
    stop('metric must be "Rsquared" or "Pval".')

  if (!missing(labels) & length(labels) != ncol(pc_matrix)){
    stop('Labels need to be same length as number of PCs supplied.')
  }

  # Calculate % variance explained
  #if(length(eigs) > 0){
  #  print(paste('Found', ncol(eigs), 'columns in eigs. Make sure you supply the
  #            full number of eigenvalues. (should = to # of features of data)'))
  #  labels <- colnames(pc_matrix) # initial labels
  #  for (i in 1:length(labels)){
  #    pc_i <- colnames(pc_matrix)[i]
  #    labels[i] <- paste(pc_i, ' (', round(eigs[i]*100/sum(eigs), 2), '%)',
  #                       sep = '')
  #  }
  #}

  # Create linear model for each PC against each meta variable and pull out r
  # squared values. Results is a nested lists object that contains i number of
  # elements equal to the number of PCs supplied. Each of these inner elements
  # are lists of j length (# of meta variables supplied in pDat) with the j-th
  # element equal to the R^2 value of the linear model fit between i-th PC ~
  # j-th meta variable.

  # bind both dataframes together
  data <- cbind(pDat[,,drop=F], pc_matrix)

  y = colnames(pc_matrix) # take PC scores names
  x = colnames(pDat)   # take out meta categorical labels

  results = list()     # initialize empty vector to fill in

  for (i in y){        # for every PC
    inner <- list()    # create an inner list of j number of elements where each
    for (j in x){      # element is the r^2 value between PCi ~ elementj
      fit <- lm(as.formula(paste(i, '~', j)), data)
      if (metric == 'Rsquared')
        inner[j] <-summary(fit)$r.squared
      if (metric == 'Pval')
        inner[j] <- summary(aov(fit))$coefficients[1,4]
    }
    results[[i]] <- inner
  }
  # Convert results into a dataframe with rows as meta variables, and columns
  # as PCs. Then return a pheatmap object.
  df <- dcast(melt(results), L2 ~ L1)
  rownames(df) <- df$L2
  df <- df[,y]
  if (!missing(labels)){
    colnames(df) <- labels  # put % variance explained in labels
  }

  colorbreaks <- seq(0, 1, by = 0.10)
  pheatmap(df, cluster_rows = F, cluster_cols = F, display_numbers = T,
           color=colorRampPalette(brewer.pal(8, "YlGnBu"))(length(colorbreaks)),
	   breaks = colorbreaks, main = paste(metric))
  if (returnMat == T) {
    return(df)
  }
}

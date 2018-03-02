cumulative <- function(pc_matrix, pDat, response){
  # creates progressively longer logistic models of the form, while adjusting 
  # for other columns in pDat, such as sex
  # response ~ pc1 + [Sex]
  # response ~ pc1 + pc2 + [Sex]
  # response ~ pc1 + pc2 + ...  + [Sex]
  # this function puts the resulting R2 values into a list
  
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
#' Tests if every column is independent of each other
#'
#' \code{pairtest} tests association between columns of a dataframe.
#' @param df a dataframe of variables, where the variables are arranged in 
#' columns, and samples or observations in rows.
#' @details Each variable is tested for their association with each other 
#' variable in simple linear regression (if at least one numerical variable), 
#' or chi squared test (if two categorical).
#' @importFrom utils combn
#' @importFrom stats as.formula
#' @importFrom stats chisq.test
#' @importFrom stats lm
#' @examples
#' ## to calculate PC association with covariates
#'
#' library(plyr)
#'
#' variables <- data.frame(
#'   row = c(1, 1, 2, 2, 3, 3),
#'   column = c(1, 1, 1, 2, 2, 2),
#'   Sex = c('m', 'm', 'm', 'f', 'f', 'f'),
#'   age = c(18, 19, 18, 27, 30, 16),
#'   ethnicity = c('AF', 'AF', 'AF', 'EU', 'EU', 'AS'))
#'   
#' pairtest(variables)
#' @export
pairtest <- function(df){
  
  # generate all combinations
  combos <- combn(ncol(df),2)
  
  plyr::adply(combos, 2, function(x) {
    # 2 categorical variables, use chisq
    if ((class(df[,x[1]]) %in% c('factor', 'character')) &&
        (class(df[,x[2]]) %in% c('factor', 'character'))){
      test <- chisq.test(df[, x[1]], df[, x[2]])
      
      out <- data.frame("Row" = colnames(df)[x[1]]
                        , "Column" = colnames(df[x[2]])
                        , "Chi.Square" = round(test$statistic,3)
                        ,  "df"= test$parameter
                        ,  "p.value" = round(test$p.value, 3)
      )
      return(out)}
    
    # 2 numerics, use regression
    if (class(df[,x[1]]) == 'numeric' &&
        class(df[,x[1]]) == 'numeric'){
      test <- summary(lm(as.formula(paste0(colnames(df)[x[1]], ' ~ ',
                                           colnames(df[x[2]]))), data = df))
      out <- data.frame("Row" = colnames(df)[x[1]]
                        , "Column" = colnames(df[x[2]])
                        , "Fstat" = round(test$fstatistic[1],3)
                        ,  "df"= test$df[2]
                        ,  "p.value" = round(test$coefficients[2,4], 3))
      return(out)
    }
    
    # 1 numeric, (1 cat), use regression
    if (sum(c(class(df[, x[1]]), class(df[,x[2]])) == 'numeric') == 1) {
      if(class(df[, x[1]]) == 'numeric') {
        test <- summary(lm(as.formula(paste0(colnames(df)[x[1]], ' ~ ',
                                             colnames(df[x[2]]))), data = df))
        
      } else {
        test <- summary(lm(as.formula(paste0(colnames(df[x[2]]), ' ~ ',
                                             colnames(df)[x[1]])), data = df))
      }
      
      out <- data.frame("Row" = colnames(df)[x[1]]
                        , "Column" = colnames(df[x[2]])
                        , "Fstat" = round(test$fstatistic[1],3)
                        ,  "df"= test$df[2]
                        ,  "p.value" = round(test$coefficients[2,4], 3))
      return(out)
    }})
}


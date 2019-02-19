#' Creates a sample sheet out of the Robinson Lab's master sample sheet
#'
#' \code{makeSampleSheet} returns the relevant idats and filepaths for a set of sample names.
#' @details Takes a vector of sample IDs and finds relevant experiment information. For example, 
#' plate, batch, Chip ID (Sentrix_ID), Row and column (Sentrix Position), and full file path are 
#' returned.
#' 
#' First it finds the sentrix Id and position from the master sample sheet (currently 450k only), 
#' then it finds the associated idats within the 450k Raw data folder, and returns the path info.
#' 
#' If a sample has been ran more than once, a dataframe will be printed along with batch / tissue 
#' information.
#' 
#' Note that the shared drive needs to be mapped to Z:
#' @param samples a vector of sample IDs
#' @param path filepath to master sample sheet, defaults to 'Z:/ROBLAB6 Infinium450k John/450k 
#' Sample Sheets/MASTER SS.xlsx'
#' 
#' @examples
#' samples <- c('PM4', 'PM30', 'PM47', 'PM123', 'PM130', 'PM252')
#' ss <- makeSampleSheet(samples)
#' rgset <- read.metharray.exp(targets = ss, verbose = T)
#' 
#' @export
makeSampleSheet <- function(samples, path = NULL){
  
  if(is.null(path)){
    path <- 'Z:/ROBLAB6 Infinium450k John/450k Sample Sheets/MASTER SS.xlsx'
  }
  if (!class(samples) %in% c('factor', 'character')){
    stop('Sample_Names must be a character or factor')
  }
  ### load master sample sheet
  master_ss <- readxl::read_xlsx(path, sheet = 1, skip = 7) 
  
  # filter to those matching sample names vector
  samples <- paste0(samples, '(?!\\d)') # ensures that PM# does not match PM## (e.g. PM3 != PM30)
  ss <- master_ss %>% dplyr::filter(grepl(paste(samples, collapse = '|'), Sample_Name, perl = T))

  # create file names and find filepath for idats
  ss$file_name <- paste0(ss$Sentrix_ID, '_', ss$Sentrix_Position,
                         '_Grn.idat')
  
  # recursively list all files in 450k Raw Data directory
  base <- 'Z:/ROBLAB6 Infinium450k John/450k Raw Data'
  x <- list.files(path = base,
                  recursive = T)
  
  y <- grep(paste(ss$file_name, collapse = '|'), x, value = T) # grep based on sample names
  #y <- y[-grep('Extra', y)] # remove extra chip matches
  
  # add file paths to ss
  df <- tibble::tibble(path = y, file_name = basename(y))
  ss <- ss %>% dplyr::left_join(df)
  
  # add the full file path as 'basename' and remove other path columns
  ss <- ss %>% 
    dplyr::mutate(Basename = file.path(base, gsub('_Grn.idat', '', ss$path))) %>%
    dplyr::select(-file_name, -path)
  
  return(ss)
}


#' Creates a sample sheet out of the Robinson Lab's master sample sheet
#'
#' \code{makeSampleSheet} find idats and returnsa data frame for a given set of samples.
#' @details Takes a vector of sample IDs and finds relevant experiment information. For example, 
#' plate, batch, Chip ID (Sentrix_ID), Row and column (Sentrix Position), and full file path are 
#' returned.
#' 
#' First it finds the sentrix Id and position from the master sample sheet (currently 450k only), 
#' then it finds the associated idats within the 450k Raw data folder, and returns file path 
#' information. If the robinson lab shared drive is mapped to a drive other than "Z:/", then you can
#' specify change the file paths with `ss_path` and `data_path`.
#' 
#' If a sample has been ran more than once, a dataframe will be printed along with batch / tissue 
#' information.
#' 
#' If there are copies of the same idats located in different directories, both will be added to the
#' data frame. Be aware that multiple copies of the same idats won't be able to be loaded with 
#' minfi::read.metharray.sheet.
#' 
#' @param samples a vector of sample IDs
#' @param ss_path file path to master sample sheet. Defaults to 'Z:/ROBLAB6 Infinium450k John/450k 
#' Sample Sheets/MASTER SS.xlsx'
#' @param data_path file path to a parent directory containing idats. Idats can be in 
#' subdirectories within the parent directory. Defaults to 'Z:/ROBLAB6 Infinium450k John/450k Raw 
#' Data'
#' 
#' @examples
#' samples <- c('PM4', 'PM30', 'PM47', 'PM123', 'PM130', 'PM252')
#' ss <- makeSampleSheet(samples)
#' rgset <- read.metharray.exp(targets = ss, verbose = T)
#' 
#' @export
makeSampleSheet <- function(samples, ss_path = NULL, data_path = NULL){
  
  if(is.null(ss_path)){
    ss_path <- 'Z:/ROBLAB6 Infinium450k John/450k Sample Sheets/MASTER SS.xlsx'
  }
  if (!class(samples) %in% c('factor', 'character')){
    stop('Sample_Names must be a character or factor')
  }
  
  ### load master sample sheet
  master_ss <- readxl::read_xlsx(ss_path, sheet = 1, skip = 7) 
  
  # filter to those matching sample names vector
  samples <- paste0(samples, '(?!\\d)') # ensures that PM# does not match PM## (e.g. PM3 != PM30)
  ss <- master_ss %>% dplyr::filter(grepl(paste(samples, collapse = '|'), Sample_Name, perl = T))

  # create file names and find filepath for idats
  ss$file_name <- paste0(ss$Sentrix_ID, '_', ss$Sentrix_Position,
                         '_Grn.idat')
  
  # recursively list all files in 450k Raw Data directory
  if(is.null(data_path)){
    data_path <- 'Z:/ROBLAB6 Infinium450k John/450k Raw Data'
  }
  
  x <- list.files(path = data_path,
                  recursive = T)
  
  y <- grep(paste(ss$file_name, collapse = '|'), x, value = T) # grep based on sample names
  #y <- y[-grep('Extra', y)] # remove extra chip matches
  
  # add file paths to ss
  df <- tibble::tibble(path = y, file_name = basename(y))
  ss <- ss %>% dplyr::left_join(df)
  
  # add the full file path as 'basename' and remove other path columns
  ss <- ss %>% 
    dplyr::mutate(Basename = file.path(data_path, gsub('_Grn.idat', '', ss$path))) %>%
    dplyr::select(-file_name, -path)
  
  return(ss)
}


#' Find file paths for IDAT files
#' 
#' \code{findsentrix} Search a directory for IDAT files matching to given set of sentrix IDs.
#' 
#' This function takes a set of sentrix IDs, searches a given directory, and returns all unique IDAT 
#' filepaths that match to each sentrix ID. 
#' 
#' A sentrix ID is a chip number provided by Illumina, followed by an underscore, folloewd by the
#' row and column identifier. Example: "9285451000_R02C01".
#' 
#' @param sentrix a vector of sentrix IDs.
#' @param directory a directory to search for IDAT files.
#' @param recursive whether to search recursively, defaults to TRUE 
#' 
#' @return a [tibble][tibble::tibble-package] returns a tibble with the given sentrix IDs and each 
#' match. Will return multiple matched idats when present, and will not return unmatched sentrix 
#' IDs.
#' 
#' @examples
#' 
#' # Example for robinson lab 
#' 
#' ## Read in master sample sheet
#' ss <- readxl::read_xlsx('Z:/ROBLAB6 Infinium450k John/Master_Sample_Sheet.xlsx')
#' 
#' ## specify idat directory
#' idat_dir <- 'Z:/ROBLAB6 Infinium450k John/EPIC Raw data/'
#' 
#' ss <- ss %>% 
#' 
#'   # Take the first 6 EPIC samples
#'   dplyr::arrange(desc(Platform)) %>% 
#'   dplyr::slice(1:6) %>% 
#'   dplyr::select(Sample_Name, Sentrix_ID, Sentrix_Position) %>%
#'   
#'   # create sentrix column
#'   dplyr::mutate(Sentrix = paste0(Sentrix_ID, '_', Sentrix_Position))
#'   
#' #join all matches, retaining unmatched and multiple matched IDs
#' ss <- ss %>%
#'   dplyr::full_join(findsentrix(sentrix = ss$Sentrix, directory = idat_dir))
#'
#' ## Now you can load in these samples with minfi::read.metharray.exp
#' rgset <- minfi::read.metharray.exp(targets = ss, verbose = TRUE)
#' 
#' @export

findsentrix <- function(sentrix, directory, recursive = TRUE) {
  
  # get all files in directory
  all_files <- list.files(path = directory, recursive = recursive, full.names = TRUE)

  # get full filenames for each sentrix 
  sentrix_files <- paste0(sentrix, '_Grn.idat')
  
  # match sentrix filenames to directory files
  matched_files <- paste(sentrix_files, collapse = '|') %>%
    grep(all_files, value = T) %>%
    stringr::str_replace('_Grn.idat', '')
  
  # store results in data frame
  results <- tibble::tibble(Sentrix = basename(matched_files),
                            Basename = matched_files) 
  
  return(results)
}

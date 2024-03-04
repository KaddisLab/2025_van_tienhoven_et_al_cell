#' Read MD5 Sums From Files
#'
#' Searches the specified directory for all files matching *md5.txt and extracts the MD5 sum for each file to a dataframe with one column as the full file path, the other the MD5 sum.
#'
#' @param directory A character string specifying the path of the directory to search within.
#'
#' @return A dataframe with two columns: `file_path` containing the full file paths, and `md5sum` containing the corresponding MD5 sums.
#'
#' @examples
#' \dontrun{
#' library(here)
#' # Assuming your md5.txt files are stored in the 'data/hpapdata' directory
#' fastq_md5sums <- read_md5sums(here::here("data/hpapdata"))
#' }
#'
#' @export
#'
#' @importFrom here here
#' @importFrom utils readLines
#' @importFrom base dirname
read_md5sums <- function(directory) {
  # Find all md5.txt files in the directory
  md5_files <- list.files(path = directory, pattern = "*md5.txt$", full.names = TRUE, recursive = TRUE)
  
  # Initialize an empty dataframe to store results
  md5_data <- data.frame(file_path = character(), md5sum = character(), stringsAsFactors = FALSE)
  
  # Loop over each md5.txt file
  for(md5_file in md5_files) {
    # Read the contents of the md5.txt file
    contents <- readLines(md5_file)
    
    # Extract md5sum and the corresponding file name
    for(line in contents) {
      parts <- strsplit(line, " +")[[1]]
      if(length(parts) >= 2) {
        md5sum <- parts[1]
        file_name <- parts[2]
        full_path <- file.path(dirname(md5_file), file_name)
        # Append to the dataframe
        md5_data <- rbind(md5_data, data.frame(file_path = full_path, md5sum = md5sum, stringsAsFactors = FALSE))
      }
    }
  }
  
  return(md5_data)
}


#' Check File MD5 Checksum
#'
#' Compares the MD5 checksum of a specified file against an expected checksum and returns `TRUE` if they match or `FALSE` otherwise.
#'
#' @param file_path The full path to the file to be checked.
#' @param expected_md5 The expected MD5 checksum to compare against.
#'
#' @return Logical `TRUE` if the actual MD5 checksum of the file matches the expected checksum, `FALSE` otherwise.
#'
#' @examples
#' \dontrun{
#' library(purrr)
#' # Assuming fastq_md5sums is a dataframe with columns file_path and md5sum
#' results <- map2_lgl(fastq_md5sums$file_path, fastq_md5sums$md5sum, check_md5sum)
#' }
#'
#' @export
#' @importFrom tools md5sum
check_md5sum <- function(fastq_md5sums) {
    file_path <- fastq_md5sums$file_path
    expected_md5 <- fastq_md5sums$md5sum
  # Return NA for missing files to distinguish from FALSE (mismatch)
  if (!file.exists(file_path)) {
    warning("File does not exist: ", file_path)
    return(NA)  
  }
  actual_md5 <- tools::md5sum(file_path)
  actual_md5 <- actual_md5[file_path]
  
  actual_md5 == expected_md5
}



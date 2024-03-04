#' Check MD5 sums of files in a path
#'
#' This function checks the MD5 sums of files in a given path against the expected MD5 sums
#' provided in corresponding txt files (matching *-md5.txt).
#'
#' @param path The path where the MD5 files and corresponding files to be checked are located.
#'
#' @return A character vector containing the paths of the files that have failed the MD5 sum check.
#'
#' @examples
#' if_interactive({
#' failed_files <- check_md5sums("data/hpapdata/HPAP-055")
#' })
#' @export
#' @importFrom future future::plan
#' @importFrom future.apply future.apply::future_lapply
#' @importFrom future.apply future.apply::future_sapply
check_md5sums <- function(path) {
    # Initialize multisession plan for parallel processing
    hprcc::init_multisession()    
    
    # Find the md5 files in the path recursively
    md5_files <- list.files(path, pattern = "-md5\\.txt$", full.names = TRUE, recursive = TRUE)
    failed_files <- character()

    # Function to process each MD5 file
    process_md5_file <- function(md5_file_path) {
        md5_contents <- readLines(md5_file_path)
        md5_file_dir <- dirname(md5_file_path)
        failed_for_this_file <- character()

        # Iterate over each line in the md5 file
        for (line in md5_contents) {
            parts <- strsplit(line, " +")[[1]]
            expected_md5 <- parts[1]
            file_name <- parts[2]
            file_path <- file.path(md5_file_dir, file_name)
            
            # Run md5sum check
            actual_md5 <- system(paste("md5sum", shQuote(file_path)), intern = TRUE)
            actual_md5 <- strsplit(actual_md5, " +")[[1]][1]
            
            # Compare expected and actual md5
            if (expected_md5 != actual_md5) {
                failed_for_this_file <- c(failed_for_this_file, file_path)
            }
        }
        
        return(failed_for_this_file)
    }

    # Use future_lapply to process each MD5 file in parallel
    results <- future.apply::future_lapply(md5_files, process_md5_file)

    # Combine all failed files from each md5 check
    failed_files <- unlist(results, use.names = FALSE)
    
    return(failed_files)
}

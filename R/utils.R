

#' Save an Object to File
#'
#' Saves an object to a file in CSV, RDS, or QS format. Creates the directory path
#' if it doesn't exist, and handles overwriting of existing files. The file format
#' is determined by the file extension in the object_path.
#'
#' @param object The object to save
#' @param object_path Character string specifying the file path where the object will be saved
#' @param overwrite Logical indicating whether to overwrite an existing file (default: TRUE)
#'
#' @return Character string containing the path to the saved file
#'
#' @importFrom readr write_csv
#' @importFrom qs qsave
#' @importFrom tools file_ext
#' @export
save_results <- function(object, object_path, overwrite = TRUE) {
    dir.create(dirname(object_path), showWarnings = FALSE, recursive = TRUE)
    message("Saving object...")

    # Check if the file exists and handle the overwrite condition
    if (file.exists(object_path) && overwrite) {
        message("File exists and overwrite is TRUE. Deleting existing file...")
        unlink(object_path)
    }

    file_type <- tools::file_ext(object_path)

    if (file_type == "csv") {
        readr::write_csv(object, object_path)
    } else if (file_type == "rds") {
        saveRDS(object, object_path)
    } else if (file_type == "qs") {
        qs::qsave(object, file = object_path)
    } else {
        stop("Unsupported file type: ", file_type)
    }

    message("Done.")
    return(object_path)
}

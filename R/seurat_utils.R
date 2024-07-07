#' Check If Seurat Object Has Been Normalized
#'
#' This function tests whether a given Seurat object has undergone normalization.
#' Normalization is assumed to have occurred if the `data` slot of the Seurat
#' object contains values. The function checks for the existence and non-emptiness
#' of the `data` slot as an indicator of normalization.
#'
#' @param seuratObj A Seurat object.
#'
#' @return Logical indicating whether the Seurat object has been normalized (`TRUE`)
#' or not (`FALSE`). 
#'
#' @examples
#' # Assuming `seuratObj` is your Seurat object
#' if (isNormalizedSeurat(seuratObj)) {
#'   print("The Seurat object has been normalized.")
#' } else {
#'   print("The Seurat object has not been normalized.")
#' }
#'
#' @export
isNormalized <- function(seurat_object) {
  activeAssay <- Seurat::DefaultAssay(seurat_object)
  
  !is.null(seurat_object[[activeAssay]]$data) && length(seurat_object[[activeAssay]]$data) > 0
}

#' Load a Seurat object from a file or variable
#'
#' This function loads a Seurat object from a file (.qs or .rds) or a variable. If the input is not a Seurat object, an error message is displayed.
#'
#' @param seurat_object A Seurat object or the path to a file containing a Seurat object (.qs or .rds).
#'
#' @return A Seurat object.
#'
#' @importFrom qs qread
#' @importFrom glue glue
#' @importFrom Seurat Seurat
#'
#' @examples
#' # Load a Seurat object from a .qs file
#' seurat_obj <- load_seurat("path/to/seurat_object.qs")
#'
#' # Load a Seurat object from a variable
#' seurat_obj <- create_seurat_object()
#' seurat_obj <- load_seurat(seurat_obj)
#'
#' @export
load_seurat <- function(seurat_object) {
  input_name <- deparse(substitute(seurat_object))
  if (!inherits(seurat_object, "Seurat")) {
    library("Seurat") |> suppressPackageStartupMessages()
    if (grepl("\\.qs$", seurat_object, ignore.case = TRUE)) {
      seurat_object <- qs::qread(seurat_object)
    } else if (grepl("\\.rds$", seurat_object, ignore.case = TRUE)) {
      seurat_object <- readRDS(seurat_object)
    } else {
      stop("Invalid Seurat object file format. Please provide the path to a .qs or .rds file.")
    }
    #Test if the object is a Seurat object. If not, show the name of variable forvided to the function and the message "...is not a Seurat object."
    if (!inherits(seurat_object, "Seurat")) {
    stop(glue::glue("{input_name} is not a Seurat object."))
    }   
  }
  return(seurat_object)
}


#' Download Files from Zenodo
#'
#' Downloads files from Zenodo given a vector of URLs and saves them to a specified directory.
#' This function is designed to be used in a `targets` pipeline, returning a vector of file paths
#' where the files were saved, suitable for use with `tar_targets(format = "file")`.
#'
#' @param urls A character vector of URLs pointing to the files to be downloaded from Zenodo.
#' @param dest_dir A character string specifying the directory where files should be saved.
#'             This directory will be created if it does not exist.
#'
#' @return A character vector of the paths to the successfully downloaded files.
#' @export
#'
#' @examples
#' urls <- c("https://zenodo.org/records/4546926/files/idx.annoy?download=1",
#'           "https://zenodo.org/records/4546926/files/ref.Rds?download=1")
#' dest_dir <- "temp" # Specify your destination directory here
#' downloaded_files <- download_zenodo_files(urls, dest_dir)
#' # Now `downloaded_files` contains paths to the downloaded files
download_zenodo_files <- function(urls, dest_dir) {
  if (!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
  
  downloaded_files <- c() # Initialize an empty vector to store file paths
  
  for (url in urls) {
    # Remove query parameters and extract the filename
    file_name <- basename(sub("\\?.*$", "", url))
    
    dest_path <- file.path(dest_dir, file_name)
    
    response <- httr::GET(url, httr::write_disk(dest_path, overwrite = TRUE))
    
    if (httr::status_code(response) == 200) {
      downloaded_files <- c(downloaded_files, dest_path)
    } else {
      warning("Failed to download: ", file_name)
    }
  }
  
  return(downloaded_files)
}



find_expression_valley <- function(seurat_object, column_name) {
    require(Seurat)
    require(ggplot2)

    # Extract the data for the specified column
    expr_values <- seurat_object[[column_name]]

    # Remove NA and zero values
    expr_values_clean <- expr_values[!is.na(expr_values) & expr_values > 0]

    # Perform kernel density estimation
    d <- density(expr_values_clean)

    # Find the valley using optimize()
    valley <- optimize(approxfun(d$x, d$y), interval = range(expr_values_clean))$minimum

    # Create the density plot with the identified valley
    p <- ggplot(data.frame(x = d$x, y = d$y), aes(x = x, y = y)) +
        geom_line() +
        geom_area(fill = "blue", alpha = 0.5) +
        geom_vline(xintercept = valley, color = "red", linetype = "dashed") +
        annotate("text",
            x = valley, y = 0, label = sprintf("Valley: %.2f", valley),
            vjust = -0.5, hjust = -0.1, color = "red"
        ) +
        theme_bw() +
        labs(
            title = paste("Density Plot of", column_name, "Expression with Identified Valley"),
            x = paste(column_name, "Expression Values"),
            y = "Density"
        )

    # Create a list to return multiple objects
    result <- list(
        valley = valley,
        plot = p,
        na_count = sum(is.na(expr_values)),
        zero_count = sum(expr_values == 0, na.rm = TRUE)
    )

    return(result)
}

# Example usage:
# result <- find_expression_valley(seurat_object, "INS_hk")
#
# print(paste("The valley between the two main peaks is at:", round(result$valley, 2)))
# print(paste("Number of NA values removed:", result$na_count))
# print(paste("Number of zero values removed:", result$zero_count))
#
# print(result$plot)


seurat_sce <- function(seurat_object) {
    # Check if the input is a Seurat object
    if (!inherits(seurat_object, "Seurat")) {
        stop("Input must be a Seurat object")
    }

    # Extract the counts matrix
    counts <- as(seurat_object[["RNA"]]$counts, "dgCMatrix")

    # Create the SingleCellExperiment object
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))

    # Add metadata to the SingleCellExperiment object
    col_data <- seurat_object[[]] |> S4Vectors::DataFrame()
    SummarizedExperiment::colData(sce) <- col_data

    return(sce)
}
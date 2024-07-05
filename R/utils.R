Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}


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

make_plots_clickable <- function() {
    require(knitr)
# Define a new knitr hook to modify the output of plots
knit_hooks$set(plot = function(x, options) {
    # Ensure the plot filename (x) is not null
    if (!is.null(x)) {
        # Use the output filename directly
        return(sprintf('<a href="%s" target="_blank"><img src="%s" alt="Plot image" style="max-width: 100%%;"></a>', x, x))
    }
})
}

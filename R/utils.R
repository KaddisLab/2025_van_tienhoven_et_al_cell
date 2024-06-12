Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}


save_results<-function(object, object_path) {
    dir.create(dirname(object_path), showWarnings = FALSE, recursive = TRUE)
    message("Saving object...")
    file_type <- tools::file_ext(object_path)
    if (file_type == "csv") {
        readr::write_csv(object, object_path, append = FALSE)
    } else if (file_type == "rds") {
        saveRDS(object, object_path)
    } else if (file_type == "qs") {
        qs::qsave(object, file = object_path)
    } else {
        stop("Unsupported file type: ", file_type)
    }
    qs::qsave(object, file = object_path)
    message("Done.")
    return(object_path)
}
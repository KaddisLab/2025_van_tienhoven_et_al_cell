Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}


save_results<-function(object, object_path) {
      # Save results
    dir.create(dirname(object_path), showWarnings = FALSE, recursive = TRUE)
    message("Saving object...")
    qs::qsave(object, file = object_path)
    message("Done.")
    return(object_path)
}
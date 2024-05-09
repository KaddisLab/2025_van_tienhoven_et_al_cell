seurat_aggregate_cell_annot <- function(...) {
    # Collect all input arguments as a list of vectors
    paths_list <- list(...)

    # Initialize an empty list to store each dataframe
    df_list <- list()

    # Process each vector of file paths in the list
    for (paths in paths_list) {
        # Modify paths to end with .csv
        csv_paths <- sapply(paths, function(path) {
            sub("\\.[^\\.]*$", ".csv", path)
        })

        # Read data into a dataframe and handle errors
        df <- tryCatch(
            {
                readr::read_csv(file = csv_paths)
            },
            error = function(e) {
                warning("Error reading CSV files: ", e$message)
                NULL # Return NULL for any read errors to filter out later
            }
        )

        # Only add the dataframe to the list if it's not NULL
        if (!is.null(df)) {
            df_list[[length(df_list) + 1]] <- df
        }
    }

    # Check if 'cell' column exists in all dataframes
    if (!all(sapply(df_list, function(df) "cell" %in% names(df)))) {
        stop("One or more data frames are missing the 'cell' column")
    }

    # Reduce the list of data frames by performing left joins on the 'cell' column
    aggregated_df <- Reduce(function(x, y) {
        dplyr::left_join(x, y, by = "cell")
    }, df_list)

    return(aggregated_df)
}

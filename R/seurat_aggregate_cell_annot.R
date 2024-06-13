seurat_aggregate_cell_annot <- function(...) {
    # Collect all input arguments as a list
    paths_list <- list(...)

    # Initialize an empty list to store each dataframe
    df_list <- list()

    # Process each element in the list
    for (element in paths_list) {
        if (is.character(element)) {
            # If the element is a character vector (paths), read CSV files into a dataframe
            csv_paths <- sapply(element, function(path) {
                sub("\\.[^\\.]*$", ".csv", path)
            })

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
        } else if (is.data.frame(element)) {
            # If the element is already a dataframe, add it to the list
            df_list[[length(df_list) + 1]] <- element
        } else {
            warning("Unsupported element type: ", class(element))
        }
    }

    # Filter out any NULL dataframes (in case of read errors)
    df_list <- Filter(Negate(is.null), df_list)

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

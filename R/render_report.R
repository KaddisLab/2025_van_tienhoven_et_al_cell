#' Render Quarto report in a temporary directory
#'
#' This function renders a Seurat report in a temporary directory. It saves the rendered report in the "reports" folder of the project directory.
#'
#' @param template A character string specifying the name of the template to be used for rendering the report
#'
#' @return A character string specifying the path of the rendered report
#'
#' @author Denis O'Meally
#' @export
render_report <- function(template) {

    session_dir <- here::here(paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = ""))
    file_name <- tools::file_path_sans_ext(basename(template))
    #temp_dir <- file.path(session_dir, glue::glue("{file_name}_report"))
    temp_input_file <- glue::glue("{session_dir}/{basename(template)}")
    dir.create(session_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Copy the template to the temp folder
    file.copy(template, session_dir)
    
    # Use quarto::quarto_render() to render the document in the temp directory
    xfun::in_dir(
        session_dir, 
        quarto::quarto_render(debug = TRUE, quiet = FALSE, 
            input = temp_input_file, output_format = "html",
            execute_dir = session_dir
        )
    )

    # Define the output file path in the temp directory
    temp_output_file <- file.path(session_dir, glue::glue("{file_name}.html"))
    # Define the final output path
    final_output_path <- here::here(glue::glue("reports/{file_name}.html"))
    # Make the output directory if it doesn't exist
    dir.create(dirname(final_output_path), showWarnings = FALSE, recursive = TRUE)
    # Move the file from the temp directory to the final output location
    file.copy(temp_output_file, final_output_path, overwrite = TRUE)
    # Remove the temporary directory and its contents
    unlink(session_dir, recursive = TRUE)
    # Return the path of the rendered file
    return(final_output_path)
}

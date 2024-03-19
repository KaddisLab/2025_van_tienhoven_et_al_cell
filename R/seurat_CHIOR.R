
seurat_CHOIR <- function(obj) {
    # https://www.choirclustering.com/articles/CHOIR.html
    require(CHOIR)

    options(parallelly.availableCores.methods = "Slurm")
    hprcc::init_multisession()

    cpus <- hprcc::slurm_allocation()$CPUs

    obj <- load_seurat(obj)

    obj <- CHOIR::CHOIR(obj, n_cores = cpus, random_seed = 42, verbose = TRUE)

    sample_id <- obj[[]]$orig.ident[1]
    choir_csv_path <- glue::glue("{analysis_cache}/CHOIR_out/{sample_id}_CHOIR.csv")
    choir_csv <- obj[[]] |>
        tibble::rownames_to_column(var = "cell") |>
        select("cell", contains("CHOIR")) 
    
    dir.create(dirname(choir_csv_path), showWarnings = FALSE, recursive = TRUE)
    write.csv(choir_csv, choir_csv_path, row.names = FALSE, quote = FALSE)

    return(choir_csv_path)
}

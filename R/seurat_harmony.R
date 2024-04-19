seurat_sketch_harmony <- function(seurat_object, batch = "batch", pancdb_metadata) {

    options(parallelly.availableCores.methods = "Slurm")
    #hprcc::init_multisession()

    # Seurat object
    seurat_object <- seurat_object |> load_seurat()
    project_name <- Seurat::Project(seurat_object)

    # Add sample metadata
    # TODO - should be in get_pancdb_metadata()
    
    meta_data <- pancdb_metadata |> 
        dplyr::filter(str_detect(reagent_kit, "10X")) |>
        dplyr::select(c(sample_id = "donor_id", "sample_sex", "sample_age", "sample_ethnicity", "ab_positive", "diabetes_status", "reagent_kit")) |>
        dplyr::mutate(batch = as.integer(as.factor(reagent_kit))) |>
        dplyr::select(-reagent_kit)
    
    seurat_object<- scCustomize::Add_Sample_Meta(seurat_object = seurat_object, meta_data = meta_data, join_by_seurat = "orig.ident", join_by_meta = "sample_id")
    
    Seurat::Idents(seurat_object) <- batch
    Seurat::DefaultAssay(seurat_object) <- "sketch"
    seurat_object[["sketch"]] <- split(seurat_object[["sketch"]], f = seurat_object$batch) 

    # Run PCA
    seurat_object <- Seurat::NormalizeData(seurat_object) |>
        Seurat::FindVariableFeatures() |>
        Seurat::ScaleData() |>
        Seurat::RunPCA()
    # Integrate across batch
    seurat_object <- Seurat::IntegrateLayers(
        object = seurat_object,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "harmony",
        verbose = FALSE
    )
    seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", return.model = TRUE)
    seurat_object <- JoinLayers(seurat_object)
    Project(seurat_object) <- glue::glue("harmony_{batch}_{project_name}")
    #--------------------------------------------------------------------------------
    # Save results
    seurat_object_path <- (glue::glue("{analysis_cache}/harmony_out/harmony_{batch}_{project_name}.qs"))
    dir.create(dirname(seurat_object_path), showWarnings = FALSE, recursive = TRUE)

    qs::qsave(seurat_object, seurat_object_path)

    return(seurat_object_path)
}

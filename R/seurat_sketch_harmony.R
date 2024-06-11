seurat_sketch_harmony <- function(seurat_object, group_by_vars = "batch", pancdb_metadata) {

    options(parallelly.availableCores.methods = "Slurm")
    #hprcc::init_multisession()

    # Seurat object
    seurat_object <- seurat_object |> load_seurat()
    project_name <- Seurat::Project(seurat_object)

    # Add sample metadata
    # TODO - should be in get_pancdb_metadata()

    meta_data <- pancdb_metadata |>
        dplyr::filter(str_detect(reagent_kit, "10X")) |>
        dplyr::select(c(sample_id = "donor_id", "sample_sex", "sample_age", "sample_ethnicity", "ab_positive", "diabetes_status", "reagent_kit", "tissue_source"))
        
    seurat_object<- scCustomize::Add_Sample_Meta(seurat_object = seurat_object, meta_data = meta_data, join_by_seurat = "orig.ident", join_by_meta = "sample_id")
    
    Seurat::DefaultAssay(seurat_object) <- "sketch"

    # Run PCA
    #TODO use SCTransform, add vars_to_regress as a parameter
    seurat_object <- Seurat::NormalizeData(seurat_object) |>
        Seurat::FindVariableFeatures() |>
        Seurat::ScaleData() |>
        Seurat::RunPCA()
    # Integrate across group.by.vars
    seurat_object <- harmony::RunHarmony(seurat_object, group.by.vars = group_by_vars)
    seurat_object <- Seurat::RunUMAP(seurat_object, dims = 1:30, reduction = "harmony", reduction.name = "umap_harmony", return.model = TRUE)
    seurat_object <- JoinLayers(seurat_object)

    new_project_name <- glue::glue("harmony__{paste(group_by_vars, collapse = '_')}__{project_name}")

    Seurat::Project(seurat_object) <- new_project_name
    #--------------------------------------------------------------------------------
    # Save results
    object_path <- (glue::glue("{analysis_cache}/harmony_out/{new_project_name}.qs"))
    save_results(seurat_object, object_path)
}

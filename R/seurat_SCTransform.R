#' Perform SCTransform on a Seurat object
#'
#' This function applies the SCTransform normalization method to a Seurat object, regressing out specified variables and saving the transformed object to a file.
#'
#' @param seurat_object A Seurat object or the path to a Seurat object.
#' @param assay The assay to use for SCTransform. Default is "RNA".
#' @param vars_to_regress A vector of variables to regress out. Default is c("percent_mt", "percent_rb").
#'
#' @return The file path to the saved transformed Seurat object.
#' @export
#'
#' @examples
#' \dontrun{
#' seurat_path <- seurat_SCTransform(seurat_object, assay = "RNA", vars_to_regress = c("percent_mt", "percent_rb"))
#' }
seurat_SCTransform <- function(seurat_object, assay = "RNA", vars_to_regress = c("percent_mt", "percent_rb")) {
    seurat_object <- load_seurat(seurat_object)
    sample_id <- Seurat::Project(seurat_object)

    seurat_object <- Seurat::SCTransform(
        object = seurat_object,
        assay = assay,
        return.only.var.genes = FALSE,
        conserve.memory = FALSE,
        residual.features = NULL,
        vars.to.regress = vars_to_regress,
        vst.flavor = "v2",
        seed.use = 42
    )

    # seurat_object_path <- glue::glue("{analysis_cache}/sctransform_out/{sample_id}_sct.qs")
    # dir.create(dirname(seurat_object_path), recursive = TRUE, showWarnings = FALSE)

    # qs::qsave(seurat_object, file = seurat_object_path)
    return(seurat_object)
}

#' Add Cell Metrics to a Seurat Object
#'
#' This function enriches a Seurat object with additional metadata columns representing
#' percentages of mitochondrial genes (MT), ribosomal proteins (RP and RL), hemoglobin genes (HB),
#' platelet-specific genes (PECAM1 and PF4), the XIST gene, and human Y chromosome genes. It is designed
#' to help in the identification of potential cell types and quality control metrics.
#'
#' @param seurat_object A Seurat object or a file path to a saved Seurat object.
#'
#' @return Seurat object with added metadata columns for each of the specified gene set percentages:
#'         `percent_mt`, `percent_rb`, `percent_hb`, `percent_pl`, `percent_xist`, and `percent_chrY`.
#'
#' @details The function calculates the percentage of expression for mitochondrial genes, ribosomal proteins,
#'          hemoglobin genes (excluding pseudogenes), platelet-specific genes, XIST for identifying female cells,
#'          and Y chromosome genes for identifying male cells. This is particularly useful for downstream
#'          analyses that require cell type identification or quality control assessments based on the expression
#'          of these gene sets. Note: The function requires the `chrY_genes` vector to be defined in the global
#'          environment or within the function's scope for the Y chromosome gene percentage calculation.
#'
#' @examples
#' # Assuming 'seuratObj' is your existing Seurat object
#' seuratObj <- seurat_add_cell_metrics(seuratObj)
#' # Now, 'seuratObj' contains additional metadata columns with the calculated percentages.
#'
#' @importFrom qs qread
#' @importFrom Seurat AddMetaData PercentageFeatureSet
#' @importFrom Matrix colSums
#' @export
seurat_add_cell_metrics <- function(seurat_object) {
        
        seurat_object <- load_seurat(seurat_object) |>
            Seurat::AddMetaData(metadata = data.frame(
                Seurat::PercentageFeatureSet(seurat_object, pattern = "^MT-"),
                Seurat::PercentageFeatureSet(seurat_object, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"),
                Seurat::PercentageFeatureSet(seurat_object, pattern = "^HB[^(P)]"),
                Seurat::PercentageFeatureSet(seurat_object, pattern = "^PECAM1|^PF4"),
                Seurat::PercentageFeatureSet(seurat_object, pattern = "^XIST"),
                ifelse(
                    sum(rownames(seurat_object[["RNA"]]$counts) %in% chrY_genes) > 2,
                    (Matrix::colSums(seurat_object[["RNA"]]$counts[rownames(seurat_object[["RNA"]]$counts) %in% chrY_genes, ]) / Matrix::colSums(seurat_object[["RNA"]]$counts)) %>% as.data.frame() * 100,
                    ifelse(
                        sum(rownames(seurat_object[["RNA"]]$counts) %in% chrY_genes) == 1,
                        Seurat::PercentageFeatureSet(seurat_object, pattern = chrY_genes[chrY_genes %in% rownames(seurat_object[["RNA"]]$counts)]),
                        Seurat::PercentageFeatureSet(seurat_object, pattern = "ZERO EXPRESSION")
                    )
                )
            ),
            col.name = c(
                "percent_mt", "percent_rb", "percent_hb", "percent_pl",
                "percent_xist", "percent_chrY"
            )
        )

        return(seurat_object)
}


#' Quick QC Filtering for Seurat Objects
#'
#' This function performs quality control (QC) filtering on a Seurat object based on specified criteria. 
#' It supports both fixed-threshold and MAD-based (median absolute deviation) filtering. 
#' For each QC metric listed in `qcFilters`, if its criterion is `NULL`, the function applies MAD-based filtering.
#' Otherwise, it uses the specified lower and upper bounds for filtering.
#' The `madMultiplier` parameter defines how many MADs away from the median a cell can be before it is filtered out.
#'
#' @param seuratObj A Seurat object on which QC filtering will be performed.
#' @param qcFilters A named list where each name corresponds to a QC metric in `seuratObj`.
#'                  The value should be a numeric vector of length two specifying the lower and upper bounds for filtering, 
#'                  or `NULL` to apply MAD-based filtering to that metric.
#' @param madMultiplier A numeric value specifying the multiplier for MAD-based filtering. 
#'                      Default is 3, indicating that cells more than 3 MADs from the median of a given metric are filtered out.
#' @param subset Logical indicating whether to return a subset of the Seurat object (TRUE) or 
#'               a logical vector indicating cells to keep (FALSE). Default is FALSE.
#'
#' @return Depending on the value of `subset`, either a subset of the `seuratObj` including only the cells 
#'         that passed all filters, or a logical vector indicating cells that passed the filtering criteria.
#'
#' @examples
#' qcFilters <- list(
#'   "percent_mt" = NULL,         # Apply MAD-based filtering for mitochondrial content
#'   "nFeature_RNA" = c(200, 8000)  # Apply fixed-threshold filtering for detected features
#' )
#' # Returns a logical vector
#' seurat_quick_qc(test_seurat_object, qcFilters) |> table()
#'
#' # Returns a subset of the Seurat object
#' seurat_quick_qc(test_seurat_object, qcFilters, subset = TRUE)
#'
#' @export
seurat_quick_qc <- function(seurat_object, qcFilters, madMultiplier = 3, subset = FALSE) {
    # Check if the input is a Seurat object, else read it from qs file
    if (!inherits(seurat_object, "Seurat")) {
        seurat_object <- qs::qread(seurat_object)
    }

  # Create a logical vector to keep track of which cells to keep
  cellsToKeep <- rep(TRUE, ncol(seurat_object))
  
    # Create a logical vector to keep track of which cells to keep
  cellsToKeep <- rep(TRUE, ncol(seurat_object))
  
  # Iterate over each QC filter and update the logical vector accordingly
  for (filterName in names(qcFilters)) {
    criteria <- qcFilters[[filterName]]
    metricValues <- Seurat::FetchData(seurat_object, vars = filterName)[,1]
    
    if (!is.null(criteria)) {
      # Apply the filter based on specified criteria if any
      lowerBound <- criteria[1]
      upperBound <- criteria[2]
      cellsToKeep <- cellsToKeep & metricValues >= lowerBound & metricValues <= upperBound
    } else {
      # Apply MAD-based filtering
      medianValue <- median(metricValues, na.rm = TRUE)
      madValue <- mad(metricValues, constant = 1, na.rm = TRUE) * madMultiplier
      lowerBound <- medianValue - madValue
      upperBound <- medianValue + madValue
      cellsToKeep <- cellsToKeep & metricValues >= lowerBound & metricValues <= upperBound
    }
  }
  
  if (subset) {
    # Subset the Seurat object to include only the cells that passed all filters
    seurat_object_subset <- subset(x = seurat_object, cells = colnames(seurat_object)[cellsToKeep])
    return(seurat_object_subset)
  } else {
    return(cellsToKeep)
  }
}


#' Human Y Chromosome Gene Names
#'
#' A character vector containing unique gene names located on the human Y chromosome. 
#' These genes were identified from a GTF file of the GRCh38 primary assembly genome, 
#' filtered to exclude genes that also appear on the X chromosome.
#'
#' @details This data object was created by importing gene information from a GTF file 
#' specific to the GRCh38 primary assembly. The process involved:
#' 
#' 1. Importing the GTF file using `rtracklayer::import`, then converting it to a data frame.
#' 2. Filtering for genes located on the Y chromosome by matching `seqnames` with "^chrY" or "^Y".
#' 3. Similarly, genes on the X chromosome were identified by matching `seqnames` with "^chrX" or "^X".
#' 4. The list of Y chromosome genes was then refined to exclude any genes that also appear on the X chromosome.
#'
#' This ensures that `chrY_genes` only contains genes uniquely identified with the Y chromosome, 
#' useful for gender determination in genomic studies or specific analyses focused on Y-linked genes.
#'
#' @usage data(chrY_genes)
#'
#' @format A character vector with each element being a unique gene name.
#' 
#' @source GTF file from: \url{{cellranger_path}/mkgtf/GRCh38.primary_assembly.genome_genes.filtered.gtf}
#' 
#' @examples
#' data(chrY_genes)
#' head(chrY_genes)
#'
#' @export
chrY_genes <- c("XGY2", "RNU6-1334P", "SRY", "RNASEH2CP1", "TOMM22P2", "ENSG00000286130", 
"RPS4Y1", "HSFY3P", "NAP1L1P2", "ENSG00000278847", "ZFY", "ZFY-AS1", 
"EEF1A1P41", "LINC00278", "AGPAT5P1", "PRRC2CP1", "TGIF2LY", 
"USP12PY", "RNF19BPY", "ENSG00000286050", "ENSG00000218410", 
"UBE2V1P3", "ENSG00000229308", "RNU6-303P", "SERBP1P2", "ENSG00000277930", 
"PCDH11Y", "RNU2-57P", "VDAC1P6", "EIF4A1P2", "KRT18P10", "SNX3P1Y", 
"RPL26P37", "TUSC2P1", "DLGAP5P1", "TTTY23B", "TSPY2", "FAM197Y9", 
"TSPY11P", "ENSG00000275352", "TSPY19P", "ENSG00000235094", "RBMY2GP", 
"LINC00280", "TTTY1B", "TTTY2B", "TTTY21B", "TTTY7", "TTTY8B", 
"ENSG00000235895", "TSPY17P", "SRIP3", "ENSG00000274365", "GOT2P5", 
"AMELY", "ATP5PFP1", "TBL1Y", "ENSG00000273906", "GPR143P", "PRKY", 
"RN7SKP282", "RNU6-941P", "RNU6-521P", "ENSG00000275280", "RBMY2HP", 
"ENSG00000223422", "TSPY12P", "ENSG00000236690", "RFTN1P1", "TTTY12", 
"ENSG00000278854", "ZNF92P1Y", "ENSG00000273863", "ZNF736P8Y", 
"BPY2DP", "ZNF736P7Y", "ZNF736P9Y", "RBMY2JP", "RBMY2KP", "TSPY24P", 
"ZNF736P6Y", "MTND6P1", "MTCYBP1", "MTND1P1", "MTND2P3", "TRIM60P3Y", 
"ZNF736P10Y", "ENSG00000275866", "LINC00279", "TTTY18", "TTTY19", 
"TTTY11", "ENSG00000228207", "ENSG00000173357", "ENSG00000273731", 
"ENSG00000230025", "OFD1P3Y", "CDY3P", "USP9YP22", "USP9YP4", 
"ELOCP4", "RBMY1A3P", "ENSG00000275310", "TTTY20", "TSPY4", "FAM197Y8", 
"TSPY8", "FAM197Y7", "ENSG00000286120", "TSPY7P", "FAM197Y6", 
"ENSG00000286173", "TSPY3", "FAM197Y5", "TSPY1", "FAM197Y4", 
"TSPY9P", "FAM197Y3", "TSPY6P", "FAM197Y2", "TSPY10", "FAM197Y1", 
"TSPY15P", "RBMY3AP", "TSPY25P", "TSPY16P", "ENSG00000232617", 
"TTTY8", "TTTY7B", "TTTY21", "TTTY2", "TTTY1", "TTTY22", "ENSG00000228379", 
"RBMY2NP", "ENSG00000274445", "TSPY18P", "TSPY13P", "TTTY23", 
"RBMY2OP", "RBMY2QP", "TSPY20P", "TSPY5P", "ENSG00000278478", 
"RNA5SP518", "RNA5SP519", "DUX4L31", "PCMTD1P1", "CDC27P2", "ENSG00000225840", 
"RNA5-8SP6", "ENSG00000273858", "CDRT15P10", "ENSG00000271365", 
"ENSG00000271309", "ENSG00000270570", "ENSG00000274231", "ENSG00000278212", 
"ANKRD20A6P", "SNX18P1Y", "DUX4L16", "DUX4L17", "DUX4L18", "DUX4L19", 
"PABPC1P5", "SLC9B1P1", "ACTR3BP1", "CHEK2P1", "MTND1P12", "ENSG00000224567", 
"ENSG00000271375", "RCC2P1", "ASS1P6", "MXRA5Y", "GYG2P1", "RPS24P1", 
"ARSFP1", "RN7SL702P", "FAM8A4P", "ARSLP1", "ARSDP1", "XGY1", 
"USP9Y", "SHROOM2P1", "MED14P1", "ENSG00000286009", "CDY4P", 
"DDX3Y", "CASKP1", "SFPQP1", "TAB3P1", "DPPA2P1", "UTY", "PSMA6P1", 
"TMSB4Y", "ANOS2P", "VCY", "VCY1B", "PNPLA4P1", "ENSG00000289706", 
"ENSG00000223517", "NLGN4Y", "AGKP1", "NLGN4Y-AS1", "MED13P1", 
"CYCSP46", "ENSG00000224518", "PUDPP1", "STSP1", "ENSG00000252689", 
"RNU6-109P", "RNU6-184P", "ENSG00000286201", "SURF6P1", "ELOCP35", 
"FAM41AY1", "ENSG00000283076", "TUBB1P2", "FAM224B", "RNA5SP520", 
"RNA5SP521", "ENSG00000273966", "CLUHP1", "TAF9P1", "ENSG00000275578", 
"CDY5P", "ELOCP36", "PRYP1", "ACTG1P2", "XKRY", "USP9YP23", "USP9YP27", 
"RNU1-128P", "TRAPPC2P3", "OFD1P1Y", "ELOCP6", "CDY2B", "CDY6P", 
"USP9YP7", "USP9YP6", "CDY7P", "USP9YP34", "USP9YP32", "CDY8P", 
"CDY2A", "ELOCP12", "OFD1P2Y", "TRAPPC2P8", "USP9YP15", "RNU1-95P", 
"USP9YP16", "XKRY2", "ACTG1P11", "PRYP2", "ELOCP26", "CDY9P", 
"TAF9P2", "CLUHP2", "FAM224A", "ENSG00000278391", "RNA5SP522", 
"RNA5SP523", "TUBB1P1", "ENSG00000282909", "FAM41AY2", "ELOCP13", 
"OFD1P4Y", "USP9YP14", "RNU1-48P", "USP9YP5", "ENSG00000251510", 
"XKRYP1", "PRYP5", "ELOCP7", "USP9YP1", "HSFY1", "GPM6BP1", "TTTY9B", 
"OFD1P5Y", "RAB9AP4", "TRAPPC2P7", "OFD1P6Y", "HSFY2", "GPM6BP2", 
"TTTY14", "USP9YP2", "ELOCP14", "PRYP6", "XKRYP2", "USP9YP10", 
"USP9YP28", "RNU1-41P", "OFD1P7Y", "MTCYBP2", "ZNF839P1", "CD24P4", 
"RNU6-255P", "ENSG00000274282", "GAPDHP19", "BCORP1", "TXLNGY", 
"ENSG00000267793", "ENSG00000260197", "KDM5D", "ENSG00000288049", 
"ENSG00000274837", "RCC2P2", "ZNF886P", "ZNF885P", "TTTY10", 
"KDM5DP1", "EIF1AY", "ENSG00000286247", "TBL1YP1", "RPS4Y2", 
"HSFY4P", "ENSG00000216844", "GAPDHP17", "TMEM167AP1", "ENSG00000254488", 
"TOMM22P1", "ENSG00000289707", "NEFLP1", "ENSG00000226918", "PRORY", 
"RBMY2EP", "ENSG00000236615", "ENSG00000288057", "RBMY2TP", "TSPY14P", 
"RBMY1HP", "ENSG00000242393", "RBMY1B", "RBMY1A1", "ENSG00000228257", 
"TTTY13", "ELOCP5", "CDY10P", "USP9YP3", "USP9YP8", "CDY11P", 
"OFD1P16Y", "ENSG00000286187", "ENSG00000229725", "RBMY1D", "ENSG00000251618", 
"RBMY1E", "ENSG00000227444", "RBMY2AP", "ENSG00000237968", "OFD1P8Y", 
"HSFY5P", "USP9YP17", "CDY12P", "ELOCP15", "PRY2", "ENSG00000236951", 
"TTTY6B", "RBMY1F", "TSPY23P", "RBMY2UP", "RBMY1KP", "TTTY5", 
"TSPY22P", "RBMY2FP", "TTTY25P", "TSPY21P", "RBMY1J", "ENSG00000235059", 
"TTTY6", "PRY", "ELOCP8", "CDY13P", "USP9YP24", "HSFY7P", "OFD1P9Y", 
"ENSG00000273589", "RBMY2BP", "ENSG00000224917", "RBMY2WP", "TTTY17A", 
"ENSG00000276829", "ZNF736P11Y", "TTTY4", "BPY2", "ENSG00000244000", 
"TRIM60P8Y", "ZNF736P3Y", "TRIM60P9Y", "DAZ1", "DAZ2", "PPP1R12BP2", 
"ENSG00000286744", "REREP1Y", "ENSG00000230073", "RBMY2CP", "ENSG00000236647", 
"OFD1P10Y", "HSFY6P", "USP9YP18", "CDY14P", "ELOCP16", "PRYP3", 
"XKRYP3", "USP9YP25", "USP9YP29", "RNU1-97P", "RAB9AP1", "TRAPPC2P9", 
"OFD1P11Y", "ELOCP9", "CDY15P", "XKRYP4", "USP9YP13", "USP9YP11", 
"TTTY3B", "USP9YP12", "RNU1-86P", "RAB9AP5", "TRAPPC2P10", "OFD1P12Y", 
"ELOCP10", "CDY1B", "CDY17P", "USP9YP35", "USP9YP31", "CDY18P", 
"ENSG00000235981", "GOLGA6L11P", "DNM1P24", "ENSG00000230977", 
"CSPG4P2Y", "CSPG4P3Y", "ENSG00000284071", "GOLGA2P2Y", "RN7SL818P", 
"UBE2Q2P4Y", "LINC00265-2P", "CICP2", "LINC00266-2P", "ENSG00000279274", 
"RBMY2XP", "TTTY17B", "TRIM60P5Y", "ZNF736P12Y", "TTTY4B", "BPY2B", 
"ENSG00000240566", "TRIM60P10Y", "ZNF736P2Y", "TRIM60P11Y", "DAZ3", 
"DAZ4", "ENSG00000278602", "ZNF736P1Y", "ENSG00000278197", "ZNF736P5Y", 
"BPY2C", "TTTY4C", "ZNF736P4Y", "ENSG00000274899", "TTTY17C", 
"RBMY2YP", "ENSG00000279115", "LINC00266-4P", "CICP1", "LINC00265-3P", 
"UBE2Q2P5Y", "GOLGA2P3Y", "RN7SL725P", "ENSG00000284380", "CSPG4P4Y", 
"CSPG4P1Y", "DNM1P48", "ENSG00000233619", "GOLGA6L16P", "ENSG00000235583", 
"CDY19P", "USP9YP36", "USP9YP33", "CDY20P", "CDY1", "ELOCP34", 
"OFD1P13Y", "TRAPPC2P5", "RAB9AP2", "USP9YP30", "RNU1-107P", 
"TTTY3", "USP9YP9", "USP9YP19", "XKRYP5", "CDY22P", "ELOCP17", 
"OFD1P18Y", "TRAPPC2P4", "RAB9AP3", "USP9YP20", "RNU1-40P", "USP9YP21", 
"XKRYP6", "PRYP4", "ELOCP11", "CDY23P", "USP9YP26", "HSFY8P", 
"OFD1P15Y", "RBMY2DP", "REREP2Y", "ENSG00000289705", "ENSG00000225876", 
"PPP1R12BP1", "RNU6-1314P", "CYCSP48", "ANKRD36P1", "ENSG00000277146", 
"TPTE2P4", "CYCSP49", "SLC25A15P1", "PARP4P1", "CCNQP2", "CTBP2P1"
)


#' Quick QC Filtering for Seurat Objects with Batch Processing
#'
#' This function performs quality control (QC) filtering on a Seurat object based on specified criteria,
#' with the ability to apply filters batch-wise according to a metadata column. It supports both fixed-threshold
#' and MAD-based (median absolute deviation) filtering. For each QC metric listed in `qcFilters`, if its criterion
#' is `NULL`, the function applies MAD-based filtering. Otherwise, it uses the specified lower and upper bounds for filtering.
#' The `madMultiplier` parameter defines how many MADs away from the median a cell can be before it is filtered out.
#' Batch processing allows for applying these filters within groups defined by a metadata column (e.g., per cluster or cell type).
#'
#' @param seuratObj A Seurat object on which QC filtering will be performed.
#' @param qcFilters A named list where each name corresponds to a QC metric in `seuratObj`.
#'                  The value should be a numeric vector of length two specifying the lower and upper bounds for filtering,
#'                  or `NULL` to apply MAD-based filtering to that metric.
#' @param madMultiplier A numeric value specifying the multiplier for MAD-based filtering.
#'                      Default is 3, indicating that cells more than 3 MADs from the median of a given metric are filtered out.
#' @param subset Logical indicating whether to return a subset of the Seurat object (TRUE) or
#'               a logical vector indicating cells to keep (FALSE). Default is FALSE.
#' @param batchVar A character string naming the metadata column to use for batch-wise filtering.
#'                 If NULL (default), filtering is applied across the entire dataset.
#'
#' @return Depending on the value of `subset`, either a subset of the `seuratObj` including only the cells
#'         that passed all filters, or a logical vector indicating cells that passed the filtering criteria.
#'
#' @examples
#' qcFilters <- list(
#'   "percent_mt" = NULL,         # Apply MAD-based filtering for mitochondrial content
#'   "nFeature_RNA" = c(200, 8000)  # Apply fixed-threshold filtering for detected features
#' )
#' # Returns a logical vector
#' seurat_quick_qc(test_seurat_object, qcFilters, batchVar = "cell_type") |> table()
#'
#' # Returns a subset of the Seurat object
#' seurat_quick_qc(test_seurat_object, qcFilters, subset = TRUE, batchVar = "cell_type")
#'
#' @export


#' Identify Outliers in a Seurat Object
#'
#' This function identifies outlier cells based on a specified metric within the metadata of a Seurat object. 
#' It can perform this identification globally or within specified batches, and supports log transformation 
#' and minimum difference adjustments for outlier determination.
#'
#' @param seuratObj A `Seurat` object containing cell metadata where the outlier analysis will be performed.
#' @param metricName The name of the metadata column to analyze for outliers.
#' @param nmads The number of median absolute deviations (MADs) from the median a value must be to be considered an outlier.
#' @param type Determines which tail(s) to consider for outlier detection: "both", "lower", or "higher".
#' @param log Logical indicating whether to log-transform the metric values before analysis.
#' @param subset Optional vector indicating a subset of cells to include in the analysis.
#' @param batchVar Optional name of the metadata column to use for batch-specific analysis.
#' @param minDiff The minimum difference from the median for a value to be considered an outlier, providing a way to adjust sensitivity.
#'
#' @return A logical vector indicating which cells in the `seuratObj` are considered outliers based on the specified metric.
#'
#' @details The function calculates medians and MADs either globally or for each batch separately if a `batchVar` is specified. 
#' It then determines outliers as those values that are a specified number of MADs away from the median, optionally considering only one tail 
#' and adjusting for log transformation and minimum difference.
#'
#' @examples
#' # Assume 'seuratObj' is your Seurat object and 'percent_mt' is a column in its metadata:
#' outliers <- findOutliersInSeurat(seuratObj, metricName = "percent_mt", nmads = 3, type = "both", log = TRUE)
#'
#' # To perform batch-specific outlier detection based on a 'batch' column in metadata:
#' outliers <- findOutliersInSeurat(seuratObj, metricName = "percent_mt", nmads = 3, 
#'                                  type = "both", log = TRUE, batchVar = "batch")
#'
#' # Visualizing outlier cells based on 'percent_mt':
#' plot(seuratObj@meta.data$percent_mt, col = ifelse(outliers, "red", "black"))
#'
#' @export
findOutliersInSeurat <- function(seuratObj, metricName, nmads = 3, type = "both", 
                                 log = FALSE, subset = NULL, batchVar = NULL, 
                                 minDiff = NA) {

metric <- seuratObj@meta.data[[metricName]]
names(metric) <- rownames(seuratObj@meta.data)
  
  if (!is.null(subset)) {
    metric <- metric[subset]
  }
  
  if (log) {
    metric <- log2(metric + 1)
  }
  
  if (!is.null(batchVar)) {
    # Calculate medians and MADs by batch
    batches <- seuratObj[[batchVar]][,1]
    uniqueBatches <- unique(batches)
    
    medians <- setNames(object = numeric(length(uniqueBatches)), nm = uniqueBatches)
    mads <- setNames(object = numeric(length(uniqueBatches)), nm = uniqueBatches)
    
    for (batch in uniqueBatches) {
      batchIndices <- which(batches == batch)
      batchMetrics <- metric[batchIndices]
      
      medians[batch] <- median(batchMetrics, na.rm = TRUE)
      mads[batch] <- mad(batchMetrics, constant = 1, na.rm = TRUE)
    }
  } else {
    # Calculate global median and MAD
    medians <- median(metric, na.rm = TRUE)
    mads <- mad(metric, constant = 1, na.rm = TRUE)
  }
  
  # Identify outliers
  outliers <- numeric(length(metric))
  names(outliers) <- names(metric)
  
  for (i in seq_along(metric)) {
    batch <- if (!is.null(batchVar)) batches[i] else 1
    thresholdLower <- medians[batch] - nmads * mads[batch]
    thresholdUpper <- medians[batch] + nmads * mads[batch]
    
    if (type == "both") {
      outliers[i] <- metric[i] < thresholdLower | metric[i] > thresholdUpper
    } else if (type == "lower") {
      outliers[i] <- metric[i] < thresholdLower
    } else if (type == "higher") {
      outliers[i] <- metric[i] > thresholdUpper
    }
    
    # Adjust for minimum difference
    if (!is.na(minDiff)) {
      if (type %in% c("both", "lower")) {
        outliers[i] <- outliers[i] & (medians[batch] - metric[i]) > minDiff
      }
      if (type %in% c("both", "higher")) {
        outliers[i] <- outliers[i] & (metric[i] - medians[batch]) > minDiff
      }
    }
  }
  
  return(as.logical(outliers))
}

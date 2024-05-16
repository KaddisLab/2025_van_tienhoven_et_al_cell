gtf_path <- "/ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

# Load the GTF file
gtf_data <- import(gtf_path, format = "gtf")

# Filter for entries related to XBP1
xbp1_data <- gtf_data[gtf_data$gene_name == "XBP1", ]

# View the filtered data
print(xbp1_data)

STAR  --runMode genomeGenerate --runThreadN 8 --genomeDir /ref_genomes/cellranger/human/star_cr_ref --genomeFastaFiles /ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/fasta/genome.fa  --sjdbGTFfile /ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf 

STAR --runThreadN 8 \
     --genomeDir /ref_genomes/cellranger/human/star_cr_ref \
     --soloType CB_UMI_Simple \
     --readFilesIn /home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/cellranger_out/10xv3/nfcore-scrnaseq-v2.5.1-refdata-gex-GRCh38-2020-A/cellranger/count/HPAP-038/outs/possorted_genome_bam.bam \
     --readFilesCommand samtools view -F 0x100 \
     --readFilesType SAM SE \
     --soloInputSAMattrBarcodeSeq CR UR \
     --soloInputSAMattrBarcodeQual CY UY \
     --soloCBwhitelist /ref_genomes/cellranger/human/cellranger_whitelists/10x_version3_whitelist.txt \
     --soloFeatures Gene GeneFull SJ Velocyto \
     --soloBarcodeReadLength 0


cd /scratch/domeally/DCD.tienhoven_scRNAseq.2024/starsolo_out/HPAP-047/;$HOME/bin/STAR --runThreadN 12 --genomeDir /ref_genomes/cellranger/human/star_cr_ref --soloType CB_UMI_Simple --readFilesIn /scratch/domeally/DCD.tienhoven_scRNAseq.2024/cellranger_out/10xv3/nfcore-scrnaseq-v2.5.1-refdata-gex-GRCh38-2020-A/cellranger/count/HPAP-047/outs/possorted_genome_bam.bam --readFilesCommand 'samtools view -F 0x100' --readFilesType SAM SE --soloInputSAMattrBarcodeSeq CR UR --soloInputSAMattrBarcodeQual CY UY --soloCBwhitelist /ref_genomes/cellranger/human/cellranger_whitelists/10x_version3_whitelist.txt --soloFeatures Gene GeneFull SJ Velocyto --soloBarcodeReadLength 0 --outFileNamePrefix /scratch/domeally/DCD.tienhoven_scRNAseq.2024/starsolo_out/HPAP-047/

 run_starsolo(
     technology = "10XV2",
     ref_dir = "/ref_genomes/cellranger/human/star_cr_ref",
     cellranger_folder = cellranger_run_folders[8]
 )


###

junction_counts <- read_tsv(
  "/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/starsolo_out/HPAP-047/Solo.out/SJ/raw/features.tsv",
  show_col_types = FALSE,
  col_names = c("chromosome", "start", "end", "strand", "intron_motif", "annotated", "unique_mapper_count", "multi_mapper_count", "max_overhang")
)
# Chromosome (chr)
# Start position of the junction
# End position of the junction
# Strand (1 for + strand, 2 for - strand, or sometimes encoded as 0 for unknown, 1 for +, 2 for -)
# Intron motif (0: non-canonical; 1: GT/AG; 2: CT/AC; 3: GC/AG; 4: CT/GC; 5: AT/AC; 6: GT/AT)
# Annotated (0: unannotated; 1: annotated based on GTF)
# Unique mapper read count
# Multi-mapper read count
# Maximum overhang:  length of the longest continuous read alignment spanning the junction on either side
# Load necessary libraries
library(rtracklayer)
library(GenomicRanges)
library(dplyr)

library(readr)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

# Load the junction data with proper column names
junction_counts <- read_tsv(
    "/home/domeally/workspaces/DCD.tienhoven_scRNAseq.2024/analysis_cache/starsolo_out/HPAP-047/Solo.out/SJ/raw/features.tsv",
    show_col_types = FALSE,
    col_names = c("chromosome", "start", "end", "strand", "intron_motif", "annotated", "unique_mapper_count", "multi_mapper_count", "max_overhang")
)
# Assuming the junction_counts data frame is already loaded and set up correctly
library(GenomicRanges)
library(rtracklayer)

# Load the GTF file to get exon information for XBP1
gtf_path <- "/ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
all_annotations <- import(gtf_path)

# Filter for XBP1 exons
xbp1_exons <- subset(all_annotations, gene_name == "XBP1" & type == "exon")
seqlevels(xbp1_exons) <- seqlevelsInUse(xbp1_exons)


# Convert the junction data frame to a GRanges object, ensure column names are correct
junctions <- GRanges(
    seqnames = as.factor(junction_counts$chromosome),
    ranges = IRanges(start = junction_counts$start, end = junction_counts$end),
    strand = ifelse(junction_counts$strand == 1, "+", "-")
)

# Define the critical splicing junction for XBP1s
critical_junction <- GRanges(
    seqnames = "chr22",
    ranges = IRanges(start = 28796122, end = 28796147),
    strand = "-"
)

# Find overlaps with the critical splicing junction
critical_overlaps <- findOverlaps(junctions, critical_junction, type = "equal")

# Count reads for the critical splicing junction (spliced isoform XBP1s)
spliced_counts <- sum(junction_counts$unique_mapper_count[queryHits(critical_overlaps)])

# Identify all overlaps with XBP1 exons
all_xbp1_overlaps <- findOverlaps(junctions, xbp1_exons)

# Adjust to extract non-critical junctions by removing critical junction indices
non_critical_indices <- setdiff(queryHits(all_xbp1_overlaps), queryHits(critical_overlaps))
non_critical_overlaps <- junctions[non_critical_indices]

# Count reads for non-critical junctions (unspliced isoform XBP1u)
unspliced_counts <- sum(junction_counts$unique_mapper_count[non_critical_indices])

# Output results
print(paste("Spliced (XBP1s) read counts:", spliced_counts))
print(paste("Unspliced (XBP1u) read counts:", unspliced_counts))

#* --- Primers from Rene --- *#
#Primers for XBP1us (unspliced) are highlighted blue: 5’ GGAGTTAAGACAGCG and 5’ CTGCAGAGGTGCACG
#Primers for XBP1s (spliced) are underlined: 5’ CTGAGTCCG/CAGCAG and 5’ GAGATGTTCTGGAGGGGTGA
#The 26nt small red nucleotides (spliced out by IRE1): GGTCTGCTGAGTCCGCAGGTGCAGGCC

spliced_exon <- "GGTCTGCTGAGTCCGCAGGTGCAGGCC"

XBP1us_primers <- c("GGAGTTAAGACAGCG", "CTGCAGAGGTGCACG")
XBP1s_primers <- c("CTGAGTCCGCAGCAG", "GAGATGTTCTGGAGGGGTGA")
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(GenomicRanges)

# Initialize an empty GRanges object to store the primer coordinates
granges_spliced <- GRanges()

# Extract the chromosome 22 sequence
chr_seq <- Hsapiens[["chr22"]]

# Search for matches of each primer in the chromosome 22 sequence
for (primer in XBP1s_primers) {
    matches <- matchPattern(primer, chr_seq)
    granges_spliced <- c(granges_spliced, GRanges(
        seqnames = Rle(rep("chr22", length(matches))),
        ranges = IRanges(start = start(matches), width = nchar(primer)),
        strand = rep("*", length(matches))
    ))
}

# Display the coordinates
granges_spliced



##################
## Annotate the hits

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Your GRanges object (after finding matches)
# granges_objects <- ...

# Load the TxDb object for annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Find overlapping genes
overlaps <- findOverlaps(granges_objects, genes(txdb))

# Create a data frame with the annotations
annotated <- data.frame(
    queryHits = queryHits(overlaps),
    subjectHits = subjectHits(overlaps),
    gene_id = mapIds(txdb, keys = subjectHits(overlaps), column = "GENEID", keytype = "TXID"),
    gene_symbol = mapIds(txdb, keys = subjectHits(overlaps), column = "SYMBOL", keytype = "TXID"),
    stringsAsFactors = FALSE
)

# Add the annotations to the GRanges object
mcols(granges_objects) <- annotated

# View the annotated GRanges object
granges_objects


#  -----------------------------------------------------------├

# Load MARVEL package
library(MARVEL)

# Load adjunct packages for selected MARVEL features
# General data processing, plotting
library(ggnewscale)
library(ggrepel)
library(reshape2)
library(plyr)
library(stringr)
library(textclean)

# Gene ontology analysis
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# ad hoc gene candidate gene analysis
library(gtools)

# Visualising splice junction location
library(GenomicRanges)
library(IRanges)
library(S4Vectors)
library(wiggleplotr)

# Load adjunct packages for this tutorial
library(Matrix)
library(data.table)
library(ggplot2)
library(gridExtra)
library(Seurat)
library(Matrix)
library(dplyr)


tar_load(seurat_sketch_750)
tar_load(starsolo_alltech)
seurat_object <- seurat_sketch_750[1]
seurat_object <- load_seurat(seurat_object)
DefaultAssay(seurat_object) <- "RNA"

starsolo_alltech <- starsolo_alltech[1]
base_path <- gsub("Summary.csv", "", starsolo_alltech)
sj_matrix_path <- file.path(base_path, "raw/matrix.mtx")
sj_barcodes_path <- file.path(base_path, "raw/barcodes.tsv")
sj_features_path <- file.path(base_path, "raw/features.tsv")
sj_counts <- readMM(sj_matrix_path) # Load the sparse matrix from Matrix Market format
sj_barcodes <- read.table(sj_barcodes_path, header = FALSE, col.names = c("barcode"))
sj_features <- read.table(sj_features_path, header = FALSE, col.names = c(
    "chromosome",
    "start",
    "end",
    "strand",
    "intron_motif",
    "annotation_status",
    "unique_mapping_reads",
    "multi_mapping_reads",
    "max_spliced_alignment_overhang"
))


counts_matrix <- GetAssayData(seurat_object, layer = "counts", assay = "RNA")
counts_matrix <- as(counts_matrix, "dgCMatrix")
cell_sums <- as.vector(Matrix::colSums(counts_matrix))
cell_sums_diag <- Diagonal(x = 1 / cell_sums)
normalized_counts <- ( counts_matrix %*% cell_sums_diag ) * 10000
seurat_object <- .clusterData(seurat_object)
df.gene.norm.pheno <- as.data.frame(seurat_object@meta.data |> rownames_to_column(var = "cell.id") |> select(cell.id, everything()))
df.gene.norm.feature <- data.frame(gene_short_name = rownames(seurat_object[["RNA"]]$counts))

pca_data <- as.data.frame(Embeddings(seurat_object, "umap")) |> rownames_to_column(var = "cell.id") |> select(cell.id, x =2, y=3)

gtf_path <- "/ref_genomes/cellranger/human/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
gtf <- as.data.frame(data.table::fread(gtf_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE))


marvel <- CreateMarvelObject.10x(
    gene.norm.matrix = as.matrix(normalized_counts),
    gene.norm.pheno = df.gene.norm.pheno,
    gene.norm.feature = df.gene.norm.feature, # Make sure this is defined correctly
    gene.count.matrix = as.matrix(counts_matrix),
    gene.count.pheno = df.gene.norm.pheno, # Assuming same metadata applies
    gene.count.feature = df.gene.norm.feature, # Assuming same gene info applies
    sj.count.matrix = sj_counts,
    sj.count.pheno = sj_barcodes,
    sj.count.feature = sj_features,
    pca = pca_data,
    gtf = gtf
)

marvel <- AnnotateGenes.10x(MarvelObject = marvel)
#marvel <- AnnotateSJ.10x(MarvelObject = marvel) #??
marvel <- ValidateSJ.10x(MarvelObject = marvel)

marvel <- CheckAlignment.10x(MarvelObject = marvel)

sample.metadata <- marvel$sample.metadata
cell.ids <- sample.metadata[, "cell.id"]

marvel <- adhocGene.TabulateExpression.Gene.10x(
    MarvelObject = marvel,
    cell.group.list = list("all" = cell.ids),
    gene_short_name = "XBP1",
    min.pct.cells = 10,
    downsample = TRUE
)

marvel$adhocGene$Expression$Gene$Plot
#Generate heatmaps based on the user input for pathway of interest
#' Heatmap generation based on pathway of interest
#'
#' This function generates a heatmap, using pheatmap, with absolute vst or rlog values from the DESeq2Results object for
#' particular pathways of interest as outputted from GSEA results. All gene names must be ENSEMBL IDs.
#'
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @param gsea_res Data frame. Results from GSEA.
#' @param deseq2_object A DESeq2Results object.
#' @param transformation Character. "vst" or "rlog" depending on user intent & sample size. Blind is set to FALSE.
#' @param annotations Character or character vector. A single or multiple variables by which the user wishes to annotate samples by (e.g. c("Treatment", "Harvest_Time"))
#' @param pathway_of_interest Character. The pathway of interest outputted from GSEA analysis (e.g. "KRAS_SIGNALLING_UP".
#'
#' @returns Heatmap.
#' @export
#'
generate_heatmap_of_interest <- function(gsea_res, deseq2_object, transformation = c("vst", "rlog"), annotations, pathway_of_interest) {

  #Match arguments for transformation based on input from function calling
  transformation <- match.arg(transformation)

  #Ensures GSEA results are a data frame & prints description pathways enriched
  gsea_df <- as.data.frame(gsea_res)
  print(gsea_df$Description)

  #Subsetting data based on defined pathway of interest
  pathway_of_interest <- pathway_of_interest
  pathway_data <- gsea_df[gsea_df$Description == pathway_of_interest, ]

  #Core enrichment string is extracted & split by each pathway
  core_enrichment_string <- pathway_data$core_enrichment
  leading_edge_genes <- unlist(strsplit(core_enrichment_string, "/"))

  print(leading_edge_genes)

  #Perform vst or rlog transformation as defined
  if (transformation == "vst") {
    ht_mp <- DESeq2::vst(deseq2_object, blind = FALSE)
  }

  else {
    ht_mp <- DESeq2::rlog(deseq2_object, blind = FALSE)
  }

  #Extract data only relating to leading edge genes of pathway of interest & define annotations for heatmap
  heatmap_data <- SummarizedExperiment::assay(ht_mp)[leading_edge_genes, ]
  annotation_col <- as.data.frame(SummarizedExperiment::colData(deseq2_object)[, annotations, drop = FALSE])

  #Map ENSEMBL IDs to gene symbols for readability
  ensembl_ids <- rownames(heatmap_data)
  message("Mapping ENSEMBL IDs to Gene Symbols...")
  hs_db <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))

  gene_names <- AnnotationDbi::mapIds(
    x = hs_db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  gene_names[is.na(gene_names)] <- ensembl_ids[is.na(gene_names)]
  rownames(heatmap_data) <- gene_names

  #Generate heatmap using pheatmap
  heatmap_graph <- pheatmap::pheatmap(heatmap_data,
                                cluster_rows = TRUE,
                                cluster_cols = TRUE,
                                show_rownames = TRUE,
                                show_colnames = TRUE,
                                annotation_col = annotation_col,
                                scale = "row",
                                main = paste0("Enriched Genes in ", pathway_of_interest),
                                fontsize = 12,
                                fontsize_row = 9,
                                fontsize_col = 12
                                )

  return(heatmap_graph)
}

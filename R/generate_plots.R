#Generate volcano plot from DESeq2 DEGs
#' Generation of volcano plot from DESeq2 DEGs
#'
#' Generates a volcano plot using ggplot using the DEGs output from DESeq2.
#' Allows for the adjusted p-value cutoff & LFC cutoff to be adjusted, &
#' whether or not to include ENSEMBL IDs or not.
#'
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @param res_df Data frame. Results from DESeq2 pipeline.
#' @param padj_cutoff Numeric. Adjusted p-value threshold. Default is 0.05.
#' @param lfc_cutoff Numeric. Log-fold change threshold. Default is 1.
#' @param map_ensembl Logical. If TRUE, ENSEMBL IDs will be mapped to gene symbols using AnnotationDbi.
#'
#' @returns Volcano plot highlighting the significantly DEGs.
#' @export
#'
generate_volcano <- function(res_df, padj_cutoff = 0.05, lfc_cutoff = 1, map_ensembl = TRUE) {

  #Ensures that the provided results are in a data frame
  res_df <- as.data.frame(res_df)

  #Remove any N/A padj values
  res_df <- res_df[!is.na(res_df[["padj"]]), ]

  #Setting parameters for LFC & adjusted p-value cutoffs.
  res_df <- res_df |>
    dplyr::mutate(expression = dplyr::case_when(
      log2FoldChange > lfc_cutoff & padj < padj_cutoff ~ "Up-regulated",
      log2FoldChange < -lfc_cutoff & padj < padj_cutoff ~ "Down-regulated",
      TRUE ~ "Not Significant"
    ))

  #Mapping ENSEMBL IDs to gene symbols for readability of volcano plot
  if (map_ensembl == TRUE) {
    message("Mapping ENSEMBL IDs to Gene Symbols...")

    hs_db <- get("org.Hs.eg.db", envir = asNamespace("org.Hs.eg.db"))

    res_df$symbol <- AnnotationDbi::mapIds(
      x = hs_db,
      keys = rownames(res_df),
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )

    res_df$symbol <- ifelse(is.na(res_df$symbol), rownames(res_df), res_df$symbol)
  }

  else {
    message("Skipping ENSEMBL ID mapping...Using IDs provided in rownames...")

    res_df$symbol <- rownames(res_df)
  }

  #Extracting the names of the top 20 genes which are most differential expressed
  #so they are shown on the volcano plot
  top_genes <- res_df |>
    dplyr::filter(padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff) |>
    dplyr::arrange(padj) |>
    utils::head(20)

  message("Generating volcano plot...")

  #Generation of volcano plot
  volcano_plot <- ggplot2::ggplot(res_df, ggplot2::aes(x = log2FoldChange, y = -log10(pvalue), color = expression)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(values = c("Up-regulated" = "red",
                                           "Not Significant" = "grey",
                                           "Down-regulated" = "blue")) +

    ggrepel::geom_text_repel(data = top_genes,
                             ggplot2::aes(label= symbol),
                             size = 3.5,
                             max.overlaps = 15,
                             box.padding = 0.5,
                             point.padding = 0.3,
                             segment.color = "black",
                             show.legend = FALSE) +

    ggplot2::geom_vline(xintercept = c(-lfc_cutoff,lfc_cutoff), linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "black") +
    ggplot2::labs(title = "Volcano Plot of DESeq2 Differentially Expressed Genes",
                  x = "log2 Fold Change",
                  y = "-log10(p-value)",
                  color = "Adjusted P-Value") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

    return(volcano_plot)
}

#Generate PCA plot from DESeq2 DEGs
#' PCA plot generation from DESeq2 object using ggplot2
#'
#' This function generates a PCA plot from the DESeq2 object by first
#' applying either a regularised log (rlog, where N < 20) or
#' variance stabilising tranformation (vst, where N > 20) & then
#' generating the plot using ggplot2.
#'
#' @param deseq2_object A DESeq2Results object. Generated for the DESeq2 pipeline.
#' @param transformation Character. "vst" or "rlog" depending on user intent & sample size. Blind is set to FALSE.
#' @param color_by Character. Variable/annotation by which PCA plot points will be coloured by.
#' @param shape_by Character. Optional, variable/annotation by which PCA points will be shaped. Default is NULL.
#'
#' @returns PCA plot coloured & shaped by annotations provided.
#' @export
#'
generate_pca <- function(deseq2_object, transformation = c("vst", "rlog"), color_by, shape_by = NULL) {

  #Match transformation variable to the argument provided
  transformation <- match.arg(transformation)

  message(paste("Applying", transformation, "transformation..."))
  if (transformation == "vst") {
    transformed_data <- DESeq2::vst(deseq2_object, blind = FALSE)
  }

  else {
    transformed_data <- DESeq2::rlog(deseq2_object, blind = FALSE)
  }

  #If shape_by is not null, it will be included for PCA plot
  if (!is.null(shape_by)) {
    intgroups <- c(color_by, shape_by)
  }

  else {
    intgroups <- color_by
  }

  message("Now extracting PCA coordinates...")
  pca_data <- DESeq2::plotPCA(transformed_data, intgroup = intgroups, returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))

  #Generation of PCA plot based on whether shape_by is NULL or not
  message("Generating PCA plot...")
  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = .data[[color_by]]))

  if (!is.null(shape_by)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(shape = .data[[shape_by]]), size = 3, alpha = 0.8)
  }

  else {
    p <- p + ggplot2::geom_point(size = 3, alpha = 0.8)
  }


  p <- p +
    ggplot2::xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggplot2::ggtitle("Principal Component Analysis") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")) +
    ggplot2::theme_minimal()

  return(p)
}

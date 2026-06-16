#Run GSEA using HALLMARK, GO, KEGG or custom gene list
#' Run GSEA using HALLMARK, GO or KEGG gene sets
#'
#' This function runs GSEA using clusterProfiler using HALLMARK,
#' GO or KEGG gene sets/terms based on the database provided.
#' It is possible to use the log2FoldChange or Wald statistic to perform this.
#' There is also an option to set a seed to ensure reproducibility of results.
#'
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#'
#' @param res_df Data frame. Data frame of results generated from the DESeq2 pipeline. Rownames must be ENSEMBL IDs.
#' @param database Character. "HALLMARK", "GO" or "KEGG" defining what terms to use for GSEA.
#' @param stat Character. "log2FoldChange" or "stat" to define which values are used for GSEA.
#' @param seed Number. Any number which will be used as the seed to ensure reproducibility. Default is NULL.
#'
#' @returns GSEA results.
#' @export
#'
run_gsea <- function(res_df, database = c("HALLMARK", "GO", "KEGG"), stat = c("log2FoldChange", "stat"), seed = NULL) {

  #Set a seed to ensure reproducibility if seed is not NULL
  if (!is.null(seed)) {
    set.seed(seed)
  }

  #Matched arguments based on function calling
  database <- match.arg(database)
  stat_col <- match.arg(stat)

  #Ensure results are in a data frame & clean data frame column containing data of interest
  res_df <- as.data.frame(res_df)
  res_df <- res_df[!is.na(res_df[[stat_col]])]
  #Subset data & order in descending order for GSEA
  ranked_genes <- res_df[[stat_col]]
  names(ranked_genes) <- rownames(res_df)

  ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  message(paste("Now undertaking GSEA using", database, "gene sets..."))

  if (database == "HALLMARK") {

    #Download HALLMARK terms
    hallmark_set <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    hallmark_geneset <- hallmark_set[, c("gs_name", "ensembl_gene")]

    #Perform GSEA with the p-value cut off of 0.05, p-value adjustment method as Benjamini-Hochberg method for calculating FDR q-values
    #so that the expected rate of false positives does not exceed 5%, with minimum gene set size being 10 & maximum being 800
    #using HALLMARK terms
    gsea_res <- clusterProfiler::GSEA(geneList = ranked_genes,
                                      TERM2GENE = hallmark_geneset,
                                      pvalueCutoff = 0.05,
                                      pAdjustMethod = "BH",
                                      verbose = FALSE,
                                      minGSSize = 10,
                                      maxGSSize = 800)

    gsea_res@result$Description <- gsub("HALLMARK_", "", gsea_res@result$Description)
  }

  else if (database == "GO") {

    #Perform GSEA with the p-value cut off of 0.05, p-value adjustment method as Benjamini-Hochberg method for calculating FDR q-values
    #so that the expected rate of false positives does not exceed 5%, with minimum gene set size being 10 & maximum being 800
    #using GO terms, biological process
    gsea_res <- clusterProfiler::gseGO(geneList = ranked_genes,
                                       ont = "BP",
                                       keyType = "ENSEMBL",
                                       minGSSize = 10,
                                       maxGSSize = 800,
                                       pvalueCutoff = 0.05,
                                       verbose = FALSE,
                                       OrgDb = "org.Hs.eg.db",
                                       pAdjustMethod = "BH")

  }

  else if (database == "KEGG") {

    message("Preparing results data frame for GSEA using KEGG gene sets...")
    #Map ENSEMBL IDs to ENTREZ IDs for KEGG gene terms to be used, ensure no IDs are duplicated
    ids <- clusterProfiler::bitr(names(ranked_genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
    dedup_ids <- ids[!duplicated(ids[c("ENSEMBL")]), ]
    kegg_df <- merge(res_df, dedup_ids, by.x = "row.names", by.y = "ENSEMBL")

    #Subset the data based on log2FoldChange or stat, have names as ENTREZ IDs & omit any N/A values
    kegg_gene_list <- kegg_df[[stat_col]]
    names(kegg_gene_list) <- kegg_df$ENTREZID
    kegg_gene_list <- stats::na.omit(kegg_gene_list)

    #Average any duplicate ENTREZ IDs, ensure they are all numeric & sort in descending order for GSEA
    kegg_gene_list <- tapply(kegg_gene_list, names(kegg_gene_list), mean)
    kegg_gene_list <- stats::setNames(as.numeric(kegg_gene_list), names(kegg_gene_list))
    kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)

    #Perform GSEA with the p-value cut off of 0.05, p-value adjustment method as Benjamini-Hochberg method for calculating FDR q-values
    #so that the expected rate of false positives does not exceed 5%, with minimum gene set size being 10 & maximum being 800
    #using KEGG terms
    message("Undergoing GSEA using KEGG gene sets...")
    gsea_res <- clusterProfiler::gseKEGG(geneList = kegg_gene_list,
                                         organism = "hsa",
                                         minGSSize = 10,
                                         maxGSSize = 800,
                                         pvalueCutoff = 0.05,
                                         pAdjustMethod = "BH",
                                         keyType = "ncbi-geneid")
  }

  return(gsea_res)

}

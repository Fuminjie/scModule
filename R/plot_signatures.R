# plot_signatures

#' Plot Gene Signature with plot_gene_signature
#' @param seurat_obj A Seurat object with UMAP or TSNE computed.
#' @param gene_set_name Name of the  gene set score in the Seurat object metadata.
#' @param reduction Reduction pattern in the Seurat object.
#' @return A UMAP or TSNE plot with coexpression of two gene signatures.
#' @export
plot_gene_signature <- function(seurat_obj, gene_set_name = "gene_set_score", reduction = "umap") {
  umap_data <- Seurat::Embeddings(seurat_obj, reduction)
  gene_set1_score <- seurat_obj@meta.data[[gene_set1_name]]
  gene_set2_score <- seurat_obj@meta.data[[gene_set2_name]]

  # Create a data frame for plotting
  plot_data <- data.frame(UMAP1 = umap_data[, 1], UMAP2 = umap_data[, 2], GeneSet1 = gene_set1_score, GeneSet2 = gene_set2_score)

  # Normalize scores for color blending
  plot_data$GeneSet1_norm <- (plot_data$GeneSet1 - min(plot_data$GeneSet1)) / (max(plot_data$GeneSet1) - min(plot_data$GeneSet1))
  plot_data$GeneSet2_norm <- (plot_data$GeneSet2 - min(plot_data$GeneSet2)) / (max(plot_data$GeneSet2) - min(plot_data$GeneSet2))

  # Create a blended color based on the two normalized scores (red for GeneSet1, blue for GeneSet2, purple for coexpression)
  plot_data$blended_color <- rgb(plot_data$GeneSet1_norm, 0, plot_data$GeneSet2_norm, maxColorValue = 1)

  # Plot
  ggplot2::ggplot(plot_data, ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = blended_color), size = 1) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(title = "UMAP with Gene Set Coexpression", x = "UMAP1", y = "UMAP2") +
    ggplot2::theme_minimal()
}

#' Plot Gene Signature with plot_two_gene_signature
#' @param seurat_obj A Seurat object with UMAP or TSNE computed.
#' @param gene_set1_name Name of the first gene set score in the Seurat object metadata.
#' @param gene_set2_name Name of the second gene set score in the Seurat object metadata.
#' @param reduction Reduction pattern in the Seurat object.
#' @return A UMAP or TSNE plot with coexpression of two gene signatures.
#' @export
plot_two_gene_signature <- function(seurat_obj, gene_set1_name = "gene_set1_score",
                                    gene_set2_name = "gene_set2_score", reduction = "umap") {
  umap_data <- Seurat::Embeddings(seurat_obj, reduction)
  gene_set1_score <- seurat_obj@meta.data[[gene_set1_name]]
  gene_set2_score <- seurat_obj@meta.data[[gene_set2_name]]

  # Create a data frame for plotting
  plot_data <- data.frame(UMAP1 = umap_data[, 1], UMAP2 = umap_data[, 2], GeneSet1 = gene_set1_score, GeneSet2 = gene_set2_score)

  # Normalize scores for color blending
  plot_data$GeneSet1_norm <- (plot_data$GeneSet1 - min(plot_data$GeneSet1)) / (max(plot_data$GeneSet1) - min(plot_data$GeneSet1))
  plot_data$GeneSet2_norm <- (plot_data$GeneSet2 - min(plot_data$GeneSet2)) / (max(plot_data$GeneSet2) - min(plot_data$GeneSet2))

  # Create a blended color based on the two normalized scores (red for GeneSet1, blue for GeneSet2, purple for coexpression)
  plot_data$blended_color <- rgb(plot_data$GeneSet1_norm, 0, plot_data$GeneSet2_norm, maxColorValue = 1)

  # Plot
  ggplot2::ggplot(plot_data, ggplot2::aes(x = UMAP1, y = UMAP2)) +
    ggplot2::geom_point(ggplot2::aes(color = blended_color), size = 1) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(title = "UMAP with Gene Set Coexpression", x = "UMAP1", y = "UMAP2") +
    ggplot2::theme_minimal()
}

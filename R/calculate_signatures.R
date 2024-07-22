# calculate_signatures.R

#' Add Gene Signature Scores to Seurat Object
#'
#' @param seurat_obj A Seurat object.
#' @param gene_set A vector of genes of interest.
#' @param name A character object to describe the name of the gene_set.
#' @return A Seurat object with added gene signature scores.
#' @export
add_gene_signature_scores <- function(seurat_obj, gene_set, Module_name = "Module_name") {
  seurat_obj <- Seurat::AddModuleScore(seurat_obj, features = list(gene_set), name = Module_name)
  return(seurat_obj)
}

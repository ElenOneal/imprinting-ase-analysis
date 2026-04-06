suppressMessages(library(dplyr))
suppressMessages(library(edgeR))

#' Filter and normalise featureCounts output for edgeR analysis.
#'
#' @param x featureCounts output object.
#' @param cross Cross direction 1 sample name pattern.
#' @param reversecross Cross direction 2 sample name pattern.
#' @param group Group vector for edgeR DGEList.
#' @param cpm.filter Minimum CPM threshold for filtering.
#' @param n.cpm.filter Minimum number of samples exceeding cpm.filter.
#' @return Normalised DGEList object.

filter_snpcounts <- function(x, cross, reversecross, group, cpm.filter, n.cpm.filter) {

  counts <- as.data.frame(x$counts)

  # Clean up column names
  cn <- gsub("\\..*", "", colnames(counts))
  colnames(counts) <- cn

  # Clean up gene names
  counts$gene_name   <- rownames(counts)
  rownames(counts)   <- gsub('.v2.1', '', counts$gene_name)
  counts             <- counts %>% dplyr::select(-gene_name)

  # Select relevant crosses
  counts <- counts %>% dplyr::select(matches(cross), matches(reversecross))
  counts <- counts[, sort(colnames(counts))]

  # Build DGEList and filter by CPM
  dge <- DGEList(counts = counts, group = group, genes = rownames(counts))
  keep <- rowSums(cpm(dge) > cpm.filter) >= n.cpm.filter
  dge  <- dge[keep, , keep.lib.sizes = FALSE]

  # TMM normalisation
  dge <- calcNormFactors(dge, method = "TMM")

  return(dge)
}

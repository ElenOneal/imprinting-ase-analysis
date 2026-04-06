suppressMessages(library(dplyr))

#' Merge per-sample allele count files into a single wide data frame.
#'
#' @param directory_path Directory containing allele count files.
#' @param file_pattern File pattern to match (e.g. '.nudatus.allelecounts.txt').
#' @param parent1name Name of parent 1 (e.g. 'dhro22').
#' @param parent2name Name of parent 2 (e.g. 'gmr2').
#'
#' @return Wide data frame with gene_name and per-sample parent1/parent2 counts.

multiMerge <- function(directory_path, file_pattern, parent1name, parent2name) {
  
  file_list  <- list.files(path = directory_path, pattern = file_pattern, full.names = TRUE)
  shortnames <- gsub(directory_path, '', file_list)
  shortnames <- gsub('/', '', shortnames)
  shortnames <- gsub(file_pattern, '', shortnames)
  
  datalist <- lapply(seq_along(file_list), function(i) {
    df <- read.table(file = file_list[i], header = FALSE, sep = "\t")
    names(df) <- c(
      "gene_name",
      paste0(parent1name, '_', shortnames[i]),
      paste0(parent2name, '_', shortnames[i])
    )
    df
  })
  
  merged <- Reduce(function(x, y) merge(x, y, by = "gene_name", all = TRUE), datalist)
  merged[is.na(merged)] <- 0
  # rownames(merged) <- merged$gene_name
  # merged$gene_name <- NULL
  
  return(merged)
}
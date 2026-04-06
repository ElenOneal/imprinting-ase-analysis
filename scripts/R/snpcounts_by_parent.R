prepare_snpcounts <- function(x,crossname,parent1,parent2){
  x <- unique(x)
  rownames(x) <- x$gene_name
  x <- x %>% dplyr::select(-gene_name)
  x[is.na(x)] <- 0
  samples <- crossname
  assign_mat_pat_dir1 <- function(x) ifelse(grepl(parent1, x), sub(parent1, "mat", x), sub(parent2, "pat", x))
  assign_mat_pat_dir2 <- function(x) ifelse(grepl(parent2, x), sub(parent2, "mat", x), sub(parent1, "pat", x))
  is_cross_Sample <- sapply(colnames(x), FUN=function(x) {strsplit(x, fixed=T, split="_")[[1]][2]}) %in% samples
  colnames(x)[is_cross_Sample] <- sapply(colnames(x)[is_cross_Sample], assign_mat_pat_dir1)
  colnames(x)[!is_cross_Sample] <- sapply(colnames(x)[!is_cross_Sample], assign_mat_pat_dir2)
  x
}

#' Convert edgeR MDS to ggplot
#'
#' @param y DGEList object
#' @param color_var column in y$samples to color by (e.g. 'group')
#' @param shape_var optional column in y$samples to shape by (e.g. 'batch')
#' @param labels named vector of group labels for legend
#' @param title plot title
#' @return ggplot object

plot_mds_gg <- function(y, color_var, shape_var = NULL, labels = NULL, title = "") {
  
  # Extract MDS coordinates
  mds <- plotMDS(y, plot = FALSE)
  
  # Build data frame
  df <- data.frame(
    x      = mds$x,
    y      = mds$y,
    sample = colnames(y),
    color  = y$samples[[color_var]]
  )
  
  if (!is.null(shape_var)) {
    df$shape <- y$samples[[shape_var]]
  }
  
  # Replace group numbers with labels if provided
  if (!is.null(labels)) {
    df$color <- factor(df$color, levels = names(labels), labels = labels)
  } else {
    df$color <- factor(df$color)
  }
  
  p <- ggplot(df, aes(x = x, y = y, color = color)) +
    { if (!is.null(shape_var)) geom_point(aes(shape = factor(shape)), size = 3)
      else geom_point(size = 3) } +
    labs(
      x     = paste0("Leading logFC dim 1"),
      y     = paste0("Leading logFC dim 2"),
      color = NULL,
      shape = NULL,
      title = title
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
  
  return(p)
}
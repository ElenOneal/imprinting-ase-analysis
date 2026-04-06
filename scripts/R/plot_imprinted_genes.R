suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

#' Plot allelic imbalance for imprinted genes in hybrid endosperm.
#'
#' @param df Data frame with columns propmat1, propmat2, mean_maternal, gene_name.
#' @param sig.genes Data frame with column design.readSourcepaternal (-1=MEG, 0=biallelic, 1=PEG).
#' @param cross Name of cross direction 1 (e.g. 'DHRO22_GMR2').
#' @param reversecross Name of cross direction 2 (e.g. 'GMR2_DHRO22').
#' @param species Species name for plot title (italicised).
#' @return ggplot object.

imprinting_plot <- function(df, sig.genes, cross, reversecross, species) {

  # Add readable imprinting status as a factor
  df$status <- factor(
    sig.genes$design.readSourcepaternal,
    levels = c(-1, 0, 1),
    labels = c("MEG", "Biallelic", "PEG")
  )

  # Count MEGs and PEGs for annotation
  n_megs <- sum(df$status == "MEG", na.rm = TRUE)
  n_pegs <- sum(df$status == "PEG", na.rm = TRUE)

  p <- ggplot(df, aes(x = propmat1 * 100, y = propmat2 * 100, color = status)) +
    geom_point(aes(alpha = status), size = 2) +
    scale_alpha_manual(values = c("MEG" = 1, "Biallelic" = 0.3, "PEG" = 1)) +
    scale_color_manual(
      values = c("MEG" = "#2166AC", "Biallelic" = "grey60", "PEG" = "#D6604D"),
      name = NULL
    ) +
    # Expected 2:1 maternal:paternal ratio in triploid endosperm
    geom_vline(xintercept = 67, color = "darkgreen", linewidth = 0.8, linetype = "dashed") +
    geom_hline(yintercept = 67, color = "darkgreen", linewidth = 0.8, linetype = "dashed") +
    # MEG and PEG counts
    annotate("text", x = 5,  y = 97,
             label = paste0("MEGs: ", n_megs),
             color = "#2166AC", hjust = 0, size = 3.5) +
    annotate("text", x = 5, y = 90,
             label = paste0("PEGs: ", n_pegs),
             color = "#D6604D", hjust = 0, size = 3.5) +
    labs(
      x = paste0("% maternal expression in ", cross),
      y = paste0("% maternal expression in ", reversecross),
      title = bquote("Allelic expression in" ~ italic(.(species)) ~ "hybrid endosperm")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      panel.grid.minor = element_blank(),
      plot.title       = element_text(size = 12)
    ) +
    guides(alpha = "none")

  return(p)
}


#' Plot histogram of mean maternal proportion across all expressed genes.
#'
#' @param df Data frame with column mean_maternal.
#' @param species Species name for plot title (italicised).
#' @return ggplot object.

maternal_histogram <- function(df, species) {

  ggplot(df, aes(x = mean_maternal)) +
    geom_histogram(binwidth = 0.05, fill = "grey70", color = "white") +
    geom_vline(xintercept = c(0.33, 0.67),
               linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
    geom_vline(xintercept = 0.5,
               linetype = "dotted", color = "grey40", linewidth = 0.8) +
    labs(
      x     = "Mean maternal proportion",
      y     = "Number of genes",
      title = bquote(italic(.(species)) ~ "endosperm allelic expression")
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

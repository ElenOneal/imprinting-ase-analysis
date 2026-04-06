suppressMessages(library(edgeR))
suppressMessages(library(limma))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(scales))

#' edgeR pre-analysis with automatic alias handling and full MDS visualization
#'
#' @param df data.frame with gene counts + column "gene_name"
#' @param samples data.frame containing sample metadata (must include "sample" and "group";
#'        "batch" optional but recommended)
#' @param min.cpm numeric CPM threshold (default 1)
#' @param N integer minimum number of samples with CPM > min.cpm
#' @param name label for plots/output
#' @param design_formula model formula (default ~ group + batch)
#' @param plot_prefix file prefix for PDF outputs
#'
#' @return list of DGEList, fit, qlf, design matrices, etc.
#'
edger_preanalysis <- function(df,
                              samples,
                              min.cpm = 1,
                              N = 3,
                              name = "edgeR_analysis",
                              design_formula = ~ group + batch,
                              plot_prefix = NULL) {
  
  stopifnot("gene_name" %in% colnames(df))
  count_cols <- setdiff(colnames(df), "gene_name")
  
  # Match sample order to count matrix
  if (!all(count_cols %in% samples$sample)) {
    stop("All column names (excluding 'gene_name') must appear in samples$sample.")
  }
  samples <- samples %>% slice(match(count_cols, sample))
  
  # Make counts matrix
  count_mat <- as.matrix(df %>% select(all_of(count_cols)))
  mode(count_mat) <- "integer"
  
  # Build DGEList
  y <- DGEList(counts = count_mat, samples = samples)
  y$samples$lib.size <- colSums(y$counts)
  
  # Filter low-expression genes
  keep <- rowSums(cpm(y) > min.cpm) >= N
  y <- y[keep, , keep.lib.sizes = FALSE]
  df_kept <- df[keep, , drop = FALSE]
  
  # Normalize
  y <- calcNormFactors(y, method = "TMM")
  
  # Adjust design formula if variables are missing
  vars <- all.vars(design_formula)
  missing <- setdiff(vars, colnames(y$samples))
  if (length(missing) > 0) {
    rhs <- paste(setdiff(vars, missing), collapse = " + ")
    design_formula <- as.formula(paste("~", rhs))
    message("Adjusted design_formula to: ", deparse(design_formula))
  }
  
  design_request <- model.matrix(design_formula, data = y$samples)
  
  # --- drop aliased (non-estimable) columns via QR decomposition ---
  qrD <- qr(design_request)
  design_full <- design_request[, qrD$pivot[seq_len(qrD$rank)], drop = FALSE]
  dropped_coef <- setdiff(colnames(design_request),
                          colnames(design_full))
  if (length(dropped_coef)) {
    message("Dropped aliased coefficients: ",
            paste(dropped_coef, collapse = ", "))
  }
  
  # Estimate dispersions + fit GLM
  y <- estimateDisp(y, design_full, robust = TRUE)
  fit <- glmQLFit(y, design_full, robust = TRUE)
  qlf <- glmQLFTest(fit)
  message(sprintf("edgeR common dispersion (%s): %.4f",
                  name, y$common.dispersion))
  
  # ------------------------ PLOTS ------------------------
  if (is.null(plot_prefix))
    plot_prefix <- gsub("[^A-Za-z0-9._-]+", "_", name)
  
  ## 1. BCV
  pdf(paste0(plot_prefix, "_BCV.pdf"), 7, 5)
  plotBCV(y, main = paste("BCV:", name))
  dev.off()
  
  ## 2. Smear
  pdf(paste0(plot_prefix, "_Smear.pdf"), 7, 5)
  plotSmear(qlf, main = paste("Smear:", name))
  abline(h = c(-2, 2), col = "gray", lty = 2)
  dev.off()
  
  ## 3. MDS (color by group)
  group_colors <- hue_pal()(nlevels(factor(y$samples$group)))
  names(group_colors) <- levels(factor(y$samples$group))
  pdf(paste0(plot_prefix, "_MDS_group.pdf"), 7, 5)
  plotMDS(y,
          col = group_colors[y$samples$group],
          pch = 16,
          main = paste("MDS: color=group", name))
  legend("topright",
         legend = names(group_colors),
         col = group_colors,
         pch = 16, cex = 0.8, title = "Group")
  dev.off()
  
  ## 4. MDS (color = group, shape = batch)
  if ("batch" %in% colnames(y$samples)) {
    batch_shapes <- c(16, 17, 15, 18, 3, 4)[seq_len(nlevels(factor(y$samples$batch)))]
    pdf(paste0(plot_prefix, "_MDS_groupBatch.pdf"), 7, 5)
    plotMDS(y,
            col = group_colors[y$samples$group],
            pch = batch_shapes[as.numeric(factor(y$samples$batch))],
            main = paste("MDS: color=group, shape=batch", name))
    legend("topleft", legend = names(group_colors),
           col = group_colors, pch = 16, cex = 0.8, title = "Group")
    legend("bottomright", legend = paste("Batch", levels(factor(y$samples$batch))),
           pch = batch_shapes, cex = 0.8, title = "Batch")
    dev.off()
  }
  
  ## 5. MDS (batch removed for visualization)
  if ("batch" %in% colnames(y$samples)) {
    lcpm <- cpm(y, log = TRUE, prior.count = 3)
    lcpm_nobatch <- removeBatchEffect(lcpm,
                                      batch = y$samples$batch,
                                      design = model.matrix(~ y$samples$group))
    pdf(paste0(plot_prefix, "_MDS_batchRemoved.pdf"), 7, 5)
    plotMDS(lcpm_nobatch,
            col = group_colors[y$samples$group],
            pch = 16,
            main = paste("MDS: batch removed for visualization -", name))
    legend("topright", legend = names(group_colors),
           col = group_colors, pch = 16, cex = 0.8, title = "Group")
    dev.off()
  }
  
  # Return everything
  invisible(list(
    y = y,
    fit = fit,
    qlf = qlf,
    keep = keep,
    lib.size = y$samples$lib.size,
    design_request = design_request,
    design_full = design_full,
    dropped_coef = dropped_coef,
    samples = y$samples,
    df_kept = df_kept
  ))
}

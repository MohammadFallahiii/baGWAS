# Manhattan (ggmanh) and Q-Q plots. Returns genomic inflation factor lambda.
.bg_plot_gwas <- function(GWAS, trait, n_snps, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  GWAS <- as.data.frame(GWAS)
  if (!all(c("CHR", "BP", "P", "SNP") %in% names(GWAS))) {
    warning("Missing columns for plotting; skipped.", call. = FALSE)
    return(NA_real_)
  }
  d <- data.frame(
    chromosome = GWAS$CHR,
    position   = as.numeric(GWAS$BP),
    P.value    = as.numeric(GWAS$P),
    SNP        = as.character(GWAS$SNP),
    stringsAsFactors = FALSE)
  d <- d[stats::complete.cases(d), , drop = FALSE]
  d <- d[is.finite(d$P.value) & d$P.value > 0, , drop = FALSE]
  if (nrow(d) == 0L) {
    warning("No valid p-values to plot.", call. = FALSE)
    return(NA_real_)
  }
  d$chromosome <- factor(d$chromosome, levels = unique(d$chromosome[
    order(suppressWarnings(as.numeric(as.character(d$chromosome))),
          as.character(d$chromosome))]))

  bonferroni <- 0.05 / n_snps
  suggestive <- 1e-5
  d$label <- ifelse(d$P.value < bonferroni, d$SNP, "")

  mpdata <- ggmanh::manhattan_data_preprocess(
    x = d, pval.colname = "P.value", chr.colname = "chromosome",
    pos.colname = "position", signif = c(bonferroni, suggestive),
    signif.col = c("red", "blue"))
  p1 <- ggmanh::manhattan_plot(
    x = mpdata, label.colname = "label", point.size = 1,
    plot.title = paste("Manhattan Plot:", trait))
  ggplot2::ggsave(file.path(out_dir, "manhattan.png"), p1,
                  width = 14, height = 6, dpi = 300)

  # Q-Q plot with lambda and 95% CI
  observed <- sort(d$P.value)
  n <- length(observed)
  k <- seq_len(n)
  expected <- k / (n + 1)
  chisq_obs <- stats::qchisq(observed, 1, lower.tail = FALSE)
  lambda <- stats::median(chisq_obs, na.rm = TRUE) / stats::qchisq(0.5, 1)

  qq <- data.frame(
    expected = -log10(expected),
    observed = -log10(observed),
    upper_ci = -log10(stats::qbeta(0.025, k, n - k + 1)),
    lower_ci = -log10(stats::qbeta(0.975, k, n - k + 1)))

  p2 <- ggplot2::ggplot(qq, ggplot2::aes(x = expected, y = observed)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
                         fill = "gray80", alpha = 0.6) +
    ggplot2::geom_point(size = 0.8, color = "black", alpha = 0.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1, color = "red",
                         linewidth = 0.8) +
    ggplot2::annotate("text", x = -Inf, y = Inf,
                      label = paste0("lambda = ", round(lambda, 3)),
                      hjust = -0.2, vjust = 2.5, size = 4) +
    ggplot2::labs(x = expression("Expected -Log"[10] * "(p-value)"),
                  y = expression("Observed -Log"[10] * "(p-value)"),
                  title = paste("QQ Plot:", trait)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      plot.background  = ggplot2::element_rect(fill = "white"),
      axis.text  = ggplot2::element_text(color = "black", size = 10),
      axis.title = ggplot2::element_text(color = "black", size = 12),
      panel.border = ggplot2::element_rect(color = "black", fill = NA,
                                           linewidth = 1),
      axis.line = ggplot2::element_blank()) +
    ggplot2::coord_cartesian(expand = FALSE)
  ggplot2::ggsave(file.path(out_dir, "qqplot.png"), p2,
                  width = 7, height = 7, dpi = 300)

  .bg_stamp("Plots saved: manhattan & qqplot | lambda=%.3f", lambda)
  lambda
}

.bg_require_cols <- function(df, cols, label) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop(label, " is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

.bg_validate <- function(strategy, organism, pheno, covariate, trait, geno,
                         group_col) {
  if (missing(organism) || !is.character(organism) || length(organism) != 1L)
    stop("'organism' must be a single character string.", call. = FALSE)
  .bg_plink_flag(organism)

  if (!is.data.frame(pheno)) stop("'pheno' must be a data frame.", call. = FALSE)
  if (!is.character(trait) || length(trait) != 1L)
    stop("'trait' must be a single column name of 'pheno'.", call. = FALSE)
  .bg_require_cols(pheno, c("FID", "IID"), "pheno")
  if (!trait %in% names(pheno))
    stop("Trait '", trait, "' not found in pheno.", call. = FALSE)
  if (!is.numeric(pheno[[trait]]))
    stop("Trait '", trait, "' must be numeric.", call. = FALSE)

  if (!is.null(covariate)) {
    if (!is.data.frame(covariate))
      stop("'covariate' must be NULL or a data frame.", call. = FALSE)
    .bg_require_cols(covariate, c("FID", "IID"), "covariate")
  }

  if (strategy == "freqGWAS") {
    if (!is.character(group_col) || length(group_col) != 1L)
      stop("'group_col' must be a single column name.", call. = FALSE)
    .bg_require_cols(pheno, group_col, "pheno")
    if (!is.null(covariate) && !group_col %in% names(covariate))
      stop("covariate needs the grouping column '", group_col, "'.",
           call. = FALSE)
  }

  geno <- .bg_normalize_geno(geno)
  for (ext in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(geno, ext)))
      stop("Missing genotype file: ", paste0(geno, ext), call. = FALSE)
  }
  invisible(geno)
}

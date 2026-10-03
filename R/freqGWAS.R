.bg_freqGWAS <- function(trait, geno, organism, pheno, covariate, group_col,
                         output_dir, threads) {
  plink1 <- .bg_plink1()
  flag <- .bg_plink_flag(organism)
  start_time <- Sys.time()

  tmp <- tempfile("baGWAS_freq_")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  tp <- function(...) file.path(tmp, ...)

  fam <- fread(paste0(geno, ".fam"), header = FALSE)
  bim <- fread(paste0(geno, ".bim"), header = FALSE)
  total_individuals <- nrow(fam)
  total_snps <- nrow(bim)

  .bg_stamp("Reading phenotype")
  Ph <- as.data.frame(pheno)
  Ph <- Ph[Ph$IID %in% fam$V2, , drop = FALSE]
  if (nrow(Ph) == 0L)
    stop("No individuals of 'pheno' were found in the .fam file.", call. = FALSE)
  Ph$.grp <- as.character(Ph[[group_col]])
  ugroups <- unique(Ph$.grp)
  grp_f <- factor(Ph$.grp, levels = ugroups)

  # group-level phenotype = mean of the trait within the group
  y <- as.numeric(tapply(as.numeric(Ph[[trait]]), grp_f, mean, na.rm = TRUE))
  y[is.nan(y)] <- NA_real_
  names(y) <- ugroups

  breed <- data.frame(FID = Ph$FID, IID = Ph$IID, Group = Ph$.grp,
                      stringsAsFactors = FALSE)
  breed_file <- tp("breed.txt")
  fwrite(breed, breed_file, sep = "\t", col.names = FALSE, quote = FALSE)

  .bg_stamp("Reading covariates")
  Xc <- NULL
  n_covariates <- 0L
  if (is.null(covariate)) {
    cat("  Covariates : 0\n")
    cat("  No covariates provided: using intercept-only model\n")
  } else {
    Cv <- as.data.frame(covariate)
    Cv$.grp <- as.character(Cv[[group_col]])
    num <- setdiff(names(Cv)[vapply(Cv, is.numeric, logical(1))],
                   c("FID", "IID", group_col, trait))
    if (length(num)) {
      cf <- factor(Cv$.grp, levels = ugroups)
      Xc <- vapply(num, function(v)
        as.numeric(tapply(Cv[[v]], cf, mean, na.rm = TRUE)), numeric(length(ugroups)))
      Xc <- matrix(Xc, nrow = length(ugroups), dimnames = list(ugroups, num))
      Xc[is.nan(Xc)] <- NA_real_
      for (j in seq_len(ncol(Xc))) {
        na <- is.na(Xc[, j])
        if (any(na)) Xc[na, j] <- stats::median(Xc[, j], na.rm = TRUE)
      }
      ok <- apply(Xc, 2, function(v) isTRUE(stats::var(v) > 0))
      Xc <- Xc[, ok, drop = FALSE]
      if (ncol(Xc) == 0L) Xc <- NULL
    }
    n_covariates <- if (is.null(Xc)) 0L else ncol(Xc)
    cat(sprintf("  Covariates : %d\n", n_covariates))
  }

  .bg_stamp("Reading genotype")
  cat(sprintf("  Total individuals  : %d\n", total_individuals))
  cat(sprintf("  Total Group/Breed  : %d\n", length(ugroups)))
  cat(sprintf("  Total SNPs         : %d\n", total_snps))

  .bg_stamp("Building AF matrix")
  .bg_run(plink1, c("--bfile", geno, flag, "--freq", "--within", breed_file,
                    "--out", tp("BreedMAF")))
  frq <- fread(tp("BreedMAF.frq.strat"), showProgress = FALSE)
  frq <- as.data.frame(frq)[, c("SNP", "CLST", "MAF")]
  wide <- as.data.frame(dcast(as.data.table(frq), CLST ~ SNP, value.var = "MAF"))
  M <- as.matrix(wide[, -1, drop = FALSE])
  rownames(M) <- as.character(wide$CLST)

  .bg_stamp("Computing kinship")
  t0 <- Sys.time()
  .bg_run(plink1, c("--bfile", geno, flag, "--cluster", "--matrix",
                    "--out", tp("IBS")))
  ibs <- as.matrix(fread(tp("IBS.mibs"), header = FALSE))
  ids <- fread(tp("IBS.mibs.id"), header = FALSE)
  grp_of <- breed$Group[match(paste(ids$V1, ids$V2),
                              paste(breed$FID, breed$IID))]
  keep <- !is.na(grp_of)
  ibs <- ibs[keep, keep, drop = FALSE]
  grp_of <- grp_of[keep]
  G <- 1 * outer(grp_of, ugroups, "==")
  cnt <- colSums(G)
  K <- crossprod(G, ibs %*% G) / outer(cnt, cnt)   # mean IBS between groups
  dimnames(K) <- list(ugroups, ugroups)
  kinship_time <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  groups <- ugroups[ugroups %in% rownames(M) & !is.na(y[ugroups]) &
                      is.finite(diag(K)[ugroups])]
  if (length(groups) < 5L)
    stop("Too few groups with phenotype and genotype data (", length(groups),
         ").", call. = FALSE)
  M <- M[groups, , drop = FALSE]
  K <- K[groups, groups, drop = FALSE]
  yv <- unname(y[groups])
  if (!is.null(Xc)) Xc <- Xc[groups, , drop = FALSE]
  Z <- diag(length(groups))
  SNPs <- colnames(M)
  n_snps <- length(SNPs)

  .bg_stamp("Running freqGWAS scan")
  scan_start <- Sys.time()
  one_snp <- function(j) {
    tryCatch({
      x <- M[, j]
      if (anyNA(x)) x[is.na(x)] <- mean(x, na.rm = TRUE)
      if (!is.finite(stats::var(x)) || stats::var(x) <= 0) return(NA_real_)
      X <- cbind(1, x, Xc)
      fit <- EMMREML::emmreml(y = yv, X = X, K = K, Z = Z, test = TRUE,
                              varbetahat = TRUE, varuhat = TRUE, PEVuhat = FALSE)
      pv <- as.matrix(fit$pvalbeta)
      as.numeric(pv[2, ncol(pv)])
    }, error = function(e) NA_real_)
  }
  csize <- max(100L, ceiling(n_snps / 100))
  starts <- seq(1L, n_snps, by = csize)
  pvals <- numeric(n_snps)
  cat(sprintf("Processing %d SNPs in %d chunks...\n", n_snps, length(starts)))
  for (s in starts) {
    idx <- s:min(n_snps, s + csize - 1L)
    pvals[idx] <- unlist(parallel::mclapply(idx, one_snp, mc.cores = threads))
    .bg_progress(max(idx), n_snps, scan_start)
  }
  cat("\n")
  gwas_time <- as.numeric(difftime(Sys.time(), scan_start, units = "mins"))

  m <- match(SNPs, bim$V2)
  GWAS <- data.frame(CHR = bim$V1[m], SNP = SNPs, BP = bim$V4[m], P = pvals,
                     stringsAsFactors = FALSE)
  err <- is.na(GWAS$P)

  out_dir <- file.path(output_dir, "Result", "freqGWAS")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (any(err)) {
    fwrite(data.frame(SNP = GWAS$SNP[err]),
           file.path(out_dir, "SNPs_with_errors.txt"), sep = "\t")
    cat(sprintf("Skipped %d SNPs due to model errors\n", sum(err)))
  }
  GWAS <- GWAS[!err, , drop = FALSE]
  fwrite(GWAS, file.path(out_dir, "baGWAS.Res"), sep = "\t", quote = FALSE)
  lambda <- .bg_plot_gwas(GWAS, trait, nrow(GWAS), out_dir)

  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  .bg_write_log(file.path(out_dir, "baGWAS.log"), "freqGWAS", organism, trait,
                out_dir, start_time, total_individuals, length(groups),
                n_covariates, total_snps, nrow(GWAS), total_time, kinship_time,
                gwas_time, lambda)
  .bg_stamp("freqGWAS finished | %s | SNPs=%d | runtime=%.2f min",
            trait, nrow(GWAS), total_time)
  .bg_annotate(GWAS, "freqGWAS", trait, organism, out_dir, lambda, total_time)
}

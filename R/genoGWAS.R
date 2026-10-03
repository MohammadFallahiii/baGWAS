# Align covariates with the FAM order; returns matrix (intercept first).
.bg_cov_matrix <- function(covariate, fam) {
  n <- nrow(fam)
  if (is.null(covariate)) {
    return(list(mat = matrix(1, n, 1, dimnames = list(NULL, "Intercept")),
                n_cov = 0L))
  }
  Cv <- as.data.frame(covariate)
  idx <- match(paste(fam[[1]], fam[[2]]), paste(Cv$FID, Cv$IID))
  if (anyNA(idx))
    warning(sum(is.na(idx)), " individual(s) have no covariate row; ",
            "values imputed with the median.", call. = FALSE)
  keep <- setdiff(names(Cv)[vapply(Cv, is.numeric, logical(1))],
                  c("FID", "IID"))
  X <- as.matrix(Cv[idx, keep, drop = FALSE])
  if (ncol(X) > 0L) {
    for (j in seq_len(ncol(X))) {
      na <- is.na(X[, j])
      if (any(na)) X[na, j] <- stats::median(X[, j], na.rm = TRUE)
    }
    ok <- apply(X, 2, function(v) isTRUE(stats::var(v) > 0))
    X <- X[, ok, drop = FALSE]
  }
  list(mat = cbind(Intercept = 1, X), n_cov = ncol(X))
}

.bg_genoGWAS <- function(trait, geno, organism, pheno, covariate, output_dir,
                         threads) {
  plink <- .bg_plink()
  gemma <- .bg_gemma()
  start_time <- Sys.time()

  tmp <- tempfile("baGWAS_geno_")
  dir.create(file.path(tmp, "output"), recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  tp <- function(...) file.path(tmp, ...)

  Ph <- as.data.frame(pheno)
  had_sex <- "sex" %in% names(Ph)
  sex <- if (had_sex) Ph$sex else rep(0, nrow(Ph))
  out_prefix <- sprintf("Geno_GWAS.%s", gsub("[^A-Za-z0-9]+", "_", trait))

  .bg_stamp("Reading phenotype")
  geno_pheno <- data.frame(FID = Ph$FID, IID = Ph$IID, sex = sex,
                           y = as.numeric(Ph[[trait]]))
  names(geno_pheno)[4] <- trait
  fwrite(geno_pheno, tp("geno.pheno"), sep = "\t", na = "-9", quote = FALSE)
  sex_args <- character()
  if (had_sex) {
    fwrite(geno_pheno[, 1:3], tp("geno.sex"), sep = "\t", quote = FALSE,
           col.names = FALSE)
    sex_args <- c("--update-sex", tp("geno.sex"))
  }

  .bg_stamp("Reading genotype")
  fam_raw <- fread(paste0(geno, ".fam"), header = FALSE)
  total_individuals <- nrow(fam_raw)
  .bg_run(plink,
          c("--bfile", geno, .bg_plink_flag(organism), "--nonfounders",
            "--allow-no-sex", "--pheno", tp("geno.pheno"),
            "--pheno-name", trait, sex_args, "--make-bed",
            "--threads", threads, "--out", tp("geno.GEMMA")))
  if (!all(file.exists(tp(paste0("geno.GEMMA", c(".bed", ".bim", ".fam"))))))
    stop("PLINK failed to create GEMMA input files.", call. = FALSE)

  fam <- fread(tp("geno.GEMMA.fam"), header = FALSE)
  analyzed_individuals <- nrow(fam)
  total_snps <- nrow(fread(tp("geno.GEMMA.bim"), header = FALSE))
  cat(sprintf("  Total individuals    : %d\n", total_individuals))
  cat(sprintf("  Included individuals : %d\n", analyzed_individuals))
  cat(sprintf("  Total SNPs           : %d\n", total_snps))

  .bg_stamp("Reading covariates")
  cm <- .bg_cov_matrix(covariate, fam)
  n_covariates <- cm$n_cov
  fwrite(as.data.frame(cm$mat), tp("geno.cov"), sep = " ", col.names = FALSE)
  cat(sprintf("  Covariates           : %d\n", n_covariates))
  if (is.null(covariate))
    cat("  No covariates provided: using intercept-only model\n")

  # phenotype in .fam order; GEMMA uses "NA" for missing values
  yy <- Ph[[trait]][match(paste(fam$V1, fam$V2), paste(Ph$FID, Ph$IID))]
  fwrite(data.frame(y = as.numeric(yy)), tp("pheno_gemma.txt"),
         sep = " ", col.names = FALSE, na = "NA")

  .bg_stamp("Computing kinship")
  t0 <- Sys.time()
  .bg_run(gemma, c("-bfile", tp("geno.GEMMA"), "-gk", "1", "-o", "RelMat",
                   "-outdir", tp("output")),
          env = .bg_gemma_env(threads))
  kinship_time <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (!file.exists(tp("output", "RelMat.cXX.txt")))
    stop("GEMMA failed to create the relationship matrix.", call. = FALSE)

  .bg_stamp("Running genoGWAS scan")
  gwas_start <- Sys.time()
  cat(sprintf("Processing %d SNPs...\n", total_snps))
  gemma_stdout <- tempfile("gemma_stdout_", fileext = ".log")
  gemma_stderr <- tempfile("gemma_stderr_", fileext = ".log")
  on.exit(unlink(c(gemma_stdout, gemma_stderr)), add = TRUE)

  job <- parallel::mcparallel({
    status <- system2(
      gemma,
      shQuote(c("-bfile", tp("geno.GEMMA"), "-k", tp("output", "RelMat.cXX.txt"),
                "-p", tp("pheno_gemma.txt"), "-c", tp("geno.cov"),
                "-lmm", "2", "-o", out_prefix, "-outdir", tp("output"))),
      stdout = gemma_stdout, stderr = gemma_stderr,
      env = .bg_gemma_env(threads))
    if (status != 0)
      stop("GEMMA failed (", status, "): ",
           paste(readLines(gemma_stderr, warn = FALSE), collapse = "\n"))
    status
  })

  read_pct <- function(path) {
    tryCatch({
      sz <- file.info(path)$size
      if (is.na(sz) || sz == 0) return(0L)
      raw <- readChar(path, sz, useBytes = TRUE)
      chunks <- strsplit(raw, "\r", fixed = TRUE)[[1]]
      m <- regmatches(chunks, regexpr("[0-9]+(?=%)", chunks, perl = TRUE))
      v <- suppressWarnings(as.integer(m))
      v <- v[!is.na(v) & v >= 0L & v <= 100L]
      if (length(v)) utils::tail(v, 1L) else 0L
    }, error = function(e) 0L)
  }

  repeat {
    res <- parallel::mccollect(job, wait = FALSE)
    if (!is.null(res)) break
    pct <- read_pct(gemma_stdout)
    .bg_progress(min(total_snps, round(total_snps * pct / 100)),
                 total_snps, gwas_start)
    Sys.sleep(1)
  }
  if (inherits(res[[1]], "try-error"))
    stop(as.character(res[[1]]), call. = FALSE)
  .bg_progress(total_snps, total_snps, gwas_start)
  cat("\n")
  gwas_time <- as.numeric(difftime(Sys.time(), gwas_start, units = "mins"))

  assoc_file <- tp("output", paste0(out_prefix, ".assoc.txt"))
  if (!file.exists(assoc_file))
    stop("GEMMA did not produce an association file.", call. = FALSE)
  assoc <- as.data.frame(fread(assoc_file))
  names(assoc) <- tolower(names(assoc))
  pick <- function(c) intersect(c, names(assoc))[1]
  cols <- c(CHR = pick(c("chr", "chrom", "chromosome")),
            SNP = pick(c("rs", "snp", "marker", "snpid")),
            BP  = pick(c("ps", "pos", "position", "bp")),
            P   = pick(c("p_wald", "p_lrt", "p_score", "pvalue", "p", "p-value")))
  if (anyNA(cols))
    stop("Missing required GEMMA columns: ",
         paste(names(cols)[is.na(cols)], collapse = ", "),
         ". Available: ", paste(names(assoc), collapse = ", "), call. = FALSE)
  res <- data.frame(CHR = assoc[[cols["CHR"]]], SNP = assoc[[cols["SNP"]]],
                    BP = assoc[[cols["BP"]]], P = assoc[[cols["P"]]],
                    stringsAsFactors = FALSE)

  out_dir <- file.path(output_dir, "Result", "genoGWAS")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(res, file.path(out_dir, "baGWAS.Res"), sep = "\t")
  lambda <- .bg_plot_gwas(res, trait, nrow(res), out_dir)

  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  .bg_write_log(file.path(out_dir, "baGWAS.log"), "genoGWAS", organism, trait,
                out_dir, start_time, total_individuals, analyzed_individuals,
                n_covariates, total_snps, nrow(res), total_time, kinship_time,
                gwas_time, lambda)
  .bg_stamp("genoGWAS finished | %s | SNPs=%d | runtime=%.2f min",
            trait, nrow(res), total_time)
  .bg_annotate(res, "genoGWAS", trait, organism, out_dir, lambda, total_time)
}

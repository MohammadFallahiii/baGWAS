.bg_stamp <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
  utils::flush.console()
}

# Run `cmd args` and report whether it can be executed (exit 126/127 mean
# "not executable"/"not found") together with the first line of its output.
.bg_probe <- function(cmd, args) {
  out <- tempfile()
  on.exit(unlink(out), add = TRUE)
  status <- suppressWarnings(tryCatch(
    system2(cmd, args, stdout = out, stderr = out),
    error = function(e) 127L))
  status <- as.integer(status)
  txt <- tryCatch(readLines(out, warn = FALSE), error = function(e) character())
  list(status = status, runs = !(status %in% c(126L, 127L)),
       version = if (length(txt)) txt[1] else "")
}

# Try each candidate command name in order; return the first one that is on
# PATH AND actually runs (and, optionally, whose version line matches `accept`).
.bg_find_tool <- function(cands, envvar = NULL, probe_args = "--version",
                          accept = NULL) {
  if (!is.null(envvar)) {
    v <- Sys.getenv(envvar, "")
    if (nzchar(v)) {
      if (file.exists(v)) return(normalizePath(v))
      if (nzchar(unname(Sys.which(v)))) return(v)
      stop("Environment variable ", envvar, " points to '", v,
           "' which was not found.", call. = FALSE)
    }
  }
  tried <- character()
  for (cmd in cands) {
    path <- unname(Sys.which(cmd))
    if (!nzchar(path)) {
      tried <- c(tried, sprintf("  %s: not found on PATH", cmd)); next
    }
    pr <- .bg_probe(cmd, probe_args)
    if (!pr$runs) {
      tried <- c(tried, sprintf("  %s: found at %s but cannot be executed (exit %d)",
                                cmd, path, pr$status)); next
    }
    if (!is.null(accept) && !grepl(accept, pr$version)) {
      tried <- c(tried, sprintf("  %s: runs but is not the required version (%s)",
                                cmd, pr$version)); next
    }
    return(cmd)   # command name, resolved by the shell via PATH
  }
  stop("No usable executable found among: ", paste(cands, collapse = ", "),
       "\n", paste(tried, collapse = "\n"),
       if (!is.null(envvar)) paste0("\nYou can also set ", envvar, "."),
       call. = FALSE)
}

# PLINK: any installed version, in this order of preference.
.bg_plink <- function() {
  .bg_find_tool(c("plink2", "plink1.9", "plink"), "BA_GWAS_PLINK")
}

# freqGWAS needs --cluster --matrix and --within (PLINK 1.x). Prefer a working
# 1.x build; fall back to plink2 with a warning.
.bg_plink1 <- function() {
  one <- tryCatch(
    .bg_find_tool(c("plink1.9", "plink"), "BA_GWAS_PLINK1", accept = "v1\\."),
    error = function(e) e)
  if (!inherits(one, "error")) return(one)
  two <- tryCatch(.bg_find_tool("plink2", NULL), error = function(e) NULL)
  if (!is.null(two)) {
    warning("No working PLINK 1.x found; using plink2. freqGWAS needs ",
            "--cluster --matrix and may fail.", call. = FALSE)
    return(two)
  }
  stop(conditionMessage(one), call. = FALSE)
}

.bg_gemma <- function() .bg_find_tool("gemma", "BA_GWAS_GEMMA", probe_args = "-h")

#' Check that the external tools needed by baGWAS are available
#'
#' Reports every PLINK/GEMMA command name baGWAS knows about: whether it is on
#' the \code{PATH}, whether it can actually be executed, and its version line.
#' Environment variables \code{BA_GWAS_PLINK}, \code{BA_GWAS_PLINK1} and
#' \code{BA_GWAS_GEMMA} can override the search.
#'
#' @return Invisibly, a data frame with columns \code{command}, \code{path},
#'   \code{runs} and \code{version}.
#' @examples
#' baGWAS_check_tools()
#' @export
baGWAS_check_tools <- function() {
  cmds <- list(plink2 = "--version", plink1.9 = "--version",
               plink = "--version", gemma = "-h")
  rows <- lapply(names(cmds), function(cmd) {
    path <- unname(Sys.which(cmd))
    if (!nzchar(path))
      return(data.frame(command = cmd, path = NA_character_, runs = FALSE,
                        version = NA_character_, stringsAsFactors = FALSE))
    pr <- .bg_probe(cmd, cmds[[cmd]])
    data.frame(command = cmd, path = path, runs = pr$runs,
               version = if (pr$runs) pr$version else
                 sprintf("cannot be executed (exit %d)", pr$status),
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  print(out, row.names = FALSE)
  invisible(out)
}

.bg_run <- function(cmd, args, env = character()) {
  out <- tempfile(); err <- tempfile()
  on.exit(unlink(c(out, err)), add = TRUE)
  status <- system2(cmd, shQuote(as.character(args)),
                    stdout = out, stderr = err, env = env)
  if (!identical(as.integer(status), 0L)) {
    stop(sprintf("Command failed (status %s): %s %s\nSTDOUT:\n%s\nSTDERR:\n%s",
                 status, cmd, paste(args, collapse = " "),
                 paste(readLines(out, warn = FALSE), collapse = "\n"),
                 paste(readLines(err, warn = FALSE), collapse = "\n")),
         call. = FALSE)
  }
  invisible(status)
}

.bg_plink_flag <- function(organism) {
  flags <- c(dog = "--dog", sheep = "--sheep", cattle = "--cow", cow = "--cow",
             horse = "--horse", mouse = "--mouse", rice = "--rice")
  f <- flags[tolower(organism)]
  if (is.na(f)) {
    stop("Unsupported organism '", organism, "'. Supported: ",
         paste(names(flags), collapse = ", "), call. = FALSE)
  }
  unname(f)
}

.bg_gemma_env <- function(n) {
  n <- as.integer(max(1L, n))
  c(sprintf("OMP_NUM_THREADS=%d", n), sprintf("OPENBLAS_NUM_THREADS=%d", n),
    sprintf("MKL_NUM_THREADS=%d", n), sprintf("BLIS_NUM_THREADS=%d", n),
    "OMP_DYNAMIC=FALSE", "OMP_PROC_BIND=spread", "OMP_PLACES=cores")
}

.bg_default_threads <- function() {
  n <- parallel::detectCores(logical = TRUE)
  if (is.na(n)) n <- 2L
  max(1L, n - 1L)
}

.bg_normalize_geno <- function(geno) {
  geno <- path.expand(as.character(geno)[1])
  geno <- sub("/+$", "", geno)
  geno <- sub("\\.(bed|bim|fam)$", "", geno)
  file.path(normalizePath(dirname(geno), winslash = "/", mustWork = FALSE),
            basename(geno))
}

.bg_progress <- function(current, total, start_time, width = 20) {
  pct <- if (total > 0) current / total else 1
  filled <- round(width * pct)
  bar <- paste0("[\033[31m", strrep("\u2588", filled), "\033[0m",
                strrep("\u2591", width - filled), "]")
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  rate <- if (elapsed > 0) current / elapsed else 0
  eta <- if (rate > 0) (total - current) / rate else 0
  cat(sprintf("\r%s %3.0f%% | %d/%d SNPs | %.0f/sec | ETA %.0fs   ",
              bar, pct * 100, as.integer(current), as.integer(total), rate, eta))
  utils::flush.console()
}

.bg_banner <- function(strategy, organism, trait, output_dir) {
  cat("-------------------------------------------------------------------------\n")
  cat(" baGWAS v1.0.5 by Mohammad.Fallahi, Dani.Fayazi-Kia and team (C) 2026\n")
  cat("-------------------------------------------------------------------------\n")
  cat("Strategy   :", strategy, "\n")
  cat("Organism   :", organism, "\n")
  cat("Trait      :", trait, "\n")
  cat("Output     :", file.path(output_dir, "Result"), "\n")
  cat("Start Time :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("-------------------------------------------------------------------------\n")
  utils::flush.console()
}

.bg_write_log <- function(path, strategy, organism, trait, output_dir,
                          start_time, total_individuals, analyzed_individuals,
                          n_covariates, total_snps, analyzed_snps,
                          total_time_min, kinship_time_min, gwas_time_min,
                          lambda) {
  f <- function(x) formatC(x, format = "f", digits = 3)
  lines <- c(
    "",
    "baGWAS v1.0.5 by Mohammad.Fallahi, Dani.Fayazi-Kia and team (C) 2026",
    "------------------------------------------------------------",
    sprintf("Strategy             : %s", strategy),
    sprintf("Organism             : %s", organism),
    sprintf("Trait                : %s", trait),
    sprintf("Output               : %s", output_dir),
    sprintf("Start Time           : %s", format(start_time, "%Y-%m-%d %H:%M:%S")),
    "",
    "Summary Statistics:",
    "------------------------------------------------------------",
    sprintf("Total individuals    : %s", total_individuals),
    sprintf("Analyzed individuals : %s", analyzed_individuals),
    sprintf("Covariates           : %s", n_covariates),
    sprintf("Total SNPs           : %s", total_snps),
    sprintf("Analyzed SNPs        : %s", analyzed_snps),
    "",
    "Computation Time:",
    "------------------------------------------------------------",
    sprintf("Total time           : %s min", f(total_time_min)),
    sprintf("Kinship              : %s min", f(kinship_time_min)),
    sprintf("GWAS scan            : %s min", f(gwas_time_min)),
    "------------------------------------------------------------",
    "",
    sprintf("Plots saved: Manhattan & qqplot | lambda=%s", f(lambda)))
  writeLines(lines, con = path)
}

.bg_annotate <- function(res, strategy, trait, organism, output_dir, lambda,
                         runtime_min) {
  res <- as.data.frame(res)
  attr(res, "strategy") <- strategy
  attr(res, "trait") <- trait
  attr(res, "organism") <- organism
  attr(res, "output_dir") <- output_dir
  attr(res, "lambda") <- lambda
  attr(res, "n_snps") <- nrow(res)
  attr(res, "runtime_min") <- runtime_min
  res
}

#' Genome-wide association analysis (genoGWAS or freqGWAS)
#'
#' Runs a GWAS using either a genotype-based strategy (PLINK + GEMMA linear
#' mixed model) or a frequency-based strategy (group allele frequencies with
#' EMMREML).
#'
#' @param strategy \code{"genoGWAS"} or \code{"freqGWAS"}.
#' @param organism One of \code{"dog"}, \code{"sheep"}, \code{"cattle"},
#'   \code{"cow"}, \code{"horse"}, \code{"mouse"}, \code{"rice"}.
#' @param pheno Data frame with columns \code{FID}, \code{IID} and the trait.
#' @param covariate \code{NULL} (intercept-only model) or a data frame with
#'   \code{FID}, \code{IID} and numeric covariate columns. For
#'   \code{freqGWAS} it must also contain \code{group_col}.
#' @param trait Name of the trait column in \code{pheno}.
#' @param geno PLINK binary file prefix (without \code{.bed/.bim/.fam}).
#' @param output_dir Directory in which the \code{Result} folder is created.
#' @param group_col Grouping column (breed/population) for \code{freqGWAS}.
#' @param threads Number of threads; default is all cores minus one.
#'
#' @return A data frame with columns \code{CHR}, \code{SNP}, \code{BP},
#'   \code{P}, with attributes \code{strategy}, \code{trait}, \code{organism},
#'   \code{output_dir}, \code{lambda}, \code{n_snps}, \code{runtime_min}.
#'   Files are written to \code{<output_dir>/Result/<strategy>/}.
#' @examples
#' \dontrun{
#' GWAS.res <- baGWAS(strategy = "genoGWAS", organism = "dog",
#'                    pheno = pheno, covariate = cov, trait = "Size",
#'                    geno = "~/Desktop/tst/geno")
#' }
#' @export
baGWAS <- function(strategy = c("genoGWAS", "freqGWAS"), organism, pheno,
                   covariate = NULL, trait, geno, output_dir = getwd(),
                   group_col = "FID", threads = NULL) {
  strategy <- match.arg(strategy)
  geno <- .bg_validate(strategy, organism, pheno, covariate, trait, geno,
                       group_col)
  if (is.null(threads)) threads <- .bg_default_threads()
  threads <- as.integer(max(1L, threads))
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  .bg_banner(strategy, organism, trait, output_dir)

  res <- if (strategy == "genoGWAS") {
    .bg_genoGWAS(trait, geno, organism, pheno, covariate, output_dir, threads)
  } else {
    .bg_freqGWAS(trait, geno, organism, pheno, covariate, group_col,
                 output_dir, threads)
  }
  invisible(res)
}

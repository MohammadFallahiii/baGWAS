\name{baGWAS}
\alias{baGWAS}
\title{Genome-wide association analysis (genoGWAS or freqGWAS)}
\usage{
baGWAS(
  strategy = c("genoGWAS", "freqGWAS"),
  organism,
  pheno,
  covariate = NULL,
  trait,
  geno,
  output_dir = getwd(),
  group_col = "FID",
  threads = NULL
)
}
\arguments{
\item{strategy}{\code{"genoGWAS"} or \code{"freqGWAS"}.}

\item{organism}{One of \code{"dog"}, \code{"sheep"}, \code{"cattle"},
\code{"cow"}, \code{"horse"}, \code{"mouse"}, \code{"rice"}.}

\item{pheno}{Data frame with columns \code{FID}, \code{IID} and the trait.}

\item{covariate}{\code{NULL} (intercept-only model) or a data frame with
\code{FID}, \code{IID} and numeric covariate columns. For \code{freqGWAS} it
must also contain \code{group_col}.}

\item{trait}{Name of the trait column in \code{pheno}.}

\item{geno}{PLINK binary file prefix (without \code{.bed/.bim/.fam}).}

\item{output_dir}{Directory in which the \code{Result} folder is created.}

\item{group_col}{Grouping column (breed/population) for \code{freqGWAS}.}

\item{threads}{Number of threads; default is all cores minus one.}
}
\value{
A data frame with columns \code{CHR}, \code{SNP}, \code{BP}, \code{P}, with
attributes \code{strategy}, \code{trait}, \code{organism}, \code{output_dir},
\code{lambda}, \code{n_snps}, \code{runtime_min}. Files are written to
\code{<output_dir>/Result/<strategy>/}.
}
\description{
Runs a GWAS using either a genotype-based strategy (PLINK + GEMMA linear
mixed model) or a frequency-based strategy (group allele frequencies with
EMMREML).
}
\examples{
\dontrun{
GWAS.res <- baGWAS(strategy = "genoGWAS", organism = "dog",
                   pheno = pheno, covariate = cov, trait = "Size",
                   geno = "~/Desktop/tst/geno")
}
}

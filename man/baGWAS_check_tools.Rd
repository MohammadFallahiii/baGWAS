\name{baGWAS_check_tools}
\alias{baGWAS_check_tools}
\title{Check that the external tools needed by baGWAS are available}
\usage{
baGWAS_check_tools()
}
\value{
Invisibly, a data frame with columns \code{command}, \code{path},
\code{runs} and \code{version}.
}
\description{
Reports every PLINK/GEMMA command name baGWAS knows about: whether it is on
the \code{PATH}, whether it can actually be executed, and its version line.
Environment variables \code{BA_GWAS_PLINK}, \code{BA_GWAS_PLINK1} and
\code{BA_GWAS_GEMMA} can override the search.
}
\examples{
baGWAS_check_tools()
}

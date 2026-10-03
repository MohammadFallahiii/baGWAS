# CRAN comments

## baGWAS 1.0.5

This file records development notes for the baGWAS source tree and is not a CRAN submission at this time.

The package requires external command-line programs (PLINK and, for `genoGWAS`, GEMMA). Automated package tests exercise R-level functionality and validation; full GWAS execution additionally requires appropriate genotype data and external executables.

Before a future CRAN submission, review external-software requirements, Unix-specific process management, platform support, dependency availability, licensing, examples, and statistical validation of both workflows.

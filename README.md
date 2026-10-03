# baGWAS

**baGWAS** is an R package for genome-wide association analysis that integrates non-matched genotype and phenotypic data at the inter-breed level through two complementary workflows:

- **`genoGWAS`** — genotype-based GWAS using PLINK-prepared binary genotype data and GEMMA linear mixed models.
- **`freqGWAS`** — frequency-based GWAS using group-level allele-frequency matrices, an IBS-derived group relationship matrix, and EMMREML.

The package provides a single main entry point, `baGWAS()`, and returns a standardized GWAS result table together with reproducible output files, diagnostic plots, and a run log.

> **Status:** Development release 1.0.5. The current release is intended for Linux/macOS (Unix-like systems) and requires external command-line software for the GWAS workflows.

## Main function

```r
baGWAS.res <- baGWAS(
  strategy = "",
  organism = "",
  pheno = pheno,
  covariate = NULL,
  trait = "",
  geno = ""
)
```

Choose one of the two strategies:

```r
strategy = "genoGWAS"
```

or

```r
strategy = "freqGWAS"
```

## Installation

### Install dependencies

`baGWAS` imports the following R packages:

- `data.table`
- `EMMREML`
- `ggplot2`
- `ggmanh`

Install the CRAN dependencies:

```r
install.packages(c("data.table", "EMMREML", "ggplot2"))
```

`ggmanh` is installed from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ggmanh", ask = FALSE, update = FALSE)
```

### Install baGWAS from GitHub

Replace `YOUR_GITHUB_USERNAME` with the account that hosts this repository:

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("YOUR_GITHUB_USERNAME/baGWAS")
```

The package can also be installed from a local source checkout:

```bash
R CMD INSTALL .
```

## External software requirements

baGWAS calls external command-line programs. These programs are not bundled with the R package and must be installed separately.

### PLINK

The genotype-based workflow can use an executable detected as `plink2`, `plink1.9`, or `plink`. The frequency-based workflow prefers a working PLINK 1.x executable because it uses the PLINK 1.x commands `--freq --within` and `--cluster --matrix`.

### GEMMA

`genoGWAS` uses GEMMA to calculate a relationship matrix and perform the linear mixed-model association scan.

### Check external tools

After installing the external programs, run:

```r
library(baGWAS)
baGWAS_check_tools()
```

You can override executable discovery with environment variables:

```r
Sys.setenv(BA_GWAS_PLINK  = "/path/to/plink2")
Sys.setenv(BA_GWAS_PLINK1 = "/path/to/plink1.9")
Sys.setenv(BA_GWAS_GEMMA  = "/path/to/gemma")
```

## Quick start

### 1. genotype-based GWAS (`genoGWAS`)

```r
library(baGWAS)

GWAS.res <- baGWAS(
  strategy = "genoGWAS",
  organism = "sheep",
  pheno = pheno,
  covariate = covariate,
  trait = "weight",
  geno = "/data/sheep/genotypes"
)
```

### 2. frequency-based GWAS (`freqGWAS`)

For `freqGWAS`, `group_col` identifies the group/breed/population used to calculate allele frequencies and group-level relationships.

```r
GWAS.res <- baGWAS(
  strategy = "freqGWAS",
  organism = "sheep",
  pheno = pheno,
  covariate = covariate,
  trait = "weight",
  geno = "/data/sheep/genotypes",
  group_col = "FID"
)
```

The `covariate` argument may be `NULL`. In that case, the model uses an intercept only.

## Input data

### Genotype data

`geno` is a PLINK binary dataset prefix. The following files must exist:

```text
genotypes.bed
genotypes.bim
genotypes.fam
```

Example:

```r
geno = "/home/user/data/sheep/genotypes"
```

A path ending in `.bed`, `.bim`, or `.fam` is normalized to the dataset prefix.

### Phenotype data

`pheno` must be a data frame containing:

```text
FID   IID   <trait>
```

The trait column must be numeric.

Example:

```r
pheno <- data.frame(
  FID = c("1", "2", "3"),
  IID = c("1", "2", "3"),
  weight = c(42.1, 45.7, 47.3)
)
```

For `freqGWAS`, `pheno` must additionally contain the grouping column specified by `group_col`.

### Covariates

For `genoGWAS`, numeric covariates are aligned to the PLINK `.fam` individual order by `FID` and `IID`. Missing numeric covariates are median-imputed and zero-variance covariates are removed.

For `freqGWAS`, numeric covariates are aggregated at the group level before entering the fixed-effect design matrix.

## Supported organisms

The current PLINK species mapping implemented by baGWAS is:

| `organism` | PLINK flag |
|---|---|
| `dog` | `--dog` |
| `sheep` | `--sheep` |
| `cattle` | `--cow` |
| `cow` | `--cow` |
| `horse` | `--horse` |
| `mouse` | `--mouse` |
| `rice` | `--rice` |

These entries describe the organism keywords implemented by the package. Users should confirm that the selected PLINK option is appropriate for their dataset and analysis.

## Output

For a run with `output_dir = "<output_dir>"`, the package creates:

```text
<output_dir>/
└── Result/
    ├── genoGWAS/
    │   ├── baGWAS.Res
    │   ├── manhattan.png
    │   ├── qqplot.png
    │   └── baGWAS.log
    └── freqGWAS/
        ├── baGWAS.Res
        ├── manhattan.png
        ├── qqplot.png
        ├── baGWAS.log
        └── SNPs_with_errors.txt   # only when SNP models fail
```

`baGWAS.Res` contains the standardized columns:

```text
CHR   SNP   BP   P
```

The returned object is a data frame with the same core result columns and run metadata stored as attributes, including strategy, trait, organism, output directory, genomic inflation statistic (`lambda`), number of analyzed SNPs, and runtime.

## Diagnostic plots

Each successful strategy writes:

- `manhattan.png` — Manhattan plot of the GWAS p-values.
- `qqplot.png` — Q-Q plot with confidence band and the genomic inflation statistic (`lambda`).

## Reproducibility and logs

baGWAS records the strategy, organism, trait, sample and SNP counts, run timing, and diagnostic information in `baGWAS.log`.

Intermediate files are created in temporary working directories and are removed after the analysis finishes.

## Statistical scope

`genoGWAS` and `freqGWAS` are distinct analytical workflows and should not be interpreted as interchangeable methods. Before interpreting results, users should evaluate genotype quality control, phenotype definition, sample structure, covariates, group definitions, missing data, model assumptions, and multiple-testing procedures.

The package performs the computational workflow; it does not replace study-specific quality control, model diagnostics, or biological interpretation.

## Testing

The repository contains `testthat` tests for package-level validation and core R functionality. Full end-to-end GWAS execution additionally requires appropriate genotype data and the external PLINK/GEMMA executables.

Run local tests with:

```r
testthat::test_local()
```

Run an R package check with:

```bash
R CMD check .
```

## Citation

Please cite **baGWAS** and the external methods/software used by your analysis. The methodological and software references are listed in [`REFERENCES.md`](REFERENCES.md).

For reproducibility, record at least:

- baGWAS version
- R version
- PLINK version
- GEMMA version
- EMMREML version
- ggmanh version

You can record the R package versions with:

```r
sessionInfo()
packageVersion("baGWAS")
packageVersion("EMMREML")
packageVersion("ggmanh")
```

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) before opening an issue or pull request. Please do not upload private genotype or phenotype data to GitHub.

## License

baGWAS is released under the MIT License. See [`LICENSE`](LICENSE).

## Repository maintenance

The `.github/workflows/R-CMD-check.yaml` workflow is configured to run automated R package checks on GitHub Actions. The project also includes issue templates, contribution guidance, citation metadata, and a security policy.

## GitHub installation example

Once the repository is public, users can install the development version with:

```r
remotes::install_github("YOUR_GITHUB_USERNAME/baGWAS")
```

Replace `YOUR_GITHUB_USERNAME` with the actual GitHub account that owns the repository.

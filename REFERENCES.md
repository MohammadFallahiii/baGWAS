# References

The following publications and software resources underpin components used by **baGWAS**.

## GEMMA

Zhou X, Stephens M. (2012). Genome-wide efficient mixed-model analysis for association studies. *Nature Genetics*, 44, 821–824. https://doi.org/10.1038/ng.2310

## PLINK

Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MAR, Bender D, Maller J, Sklar P, de Bakker PIW, Daly MJ, Sham PC. (2007). PLINK: a tool set for whole-genome association and population-based linkage analyses. *The American Journal of Human Genetics*, 81(3), 559–575. https://doi.org/10.1086/519795

Chang CC, Chow CC, Tellier LCAM, Vattikuti S, Purcell SM, Lee JJ. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience*, 4, 7. https://doi.org/10.1186/s13742-015-0047-8

## EMMREML

Akdemir D, Godfrey OU. (2015). EMMREML: Fitting Mixed Models with Known Covariance Structures. R package version 3.1. https://doi.org/10.32614/CRAN.package.EMMREML

## ggmanh

Lee J, Zheng X. `ggmanh`: Visualization Tool for GWAS Result. Bioconductor. DOI: https://doi.org/10.18129/B9.bioc.ggmanh

## Reproducibility

When reporting a baGWAS analysis, record the exact versions of:

- baGWAS
- R
- PLINK
- GEMMA
- EMMREML
- ggmanh

Example:

```r
sessionInfo()
packageVersion("baGWAS")
packageVersion("EMMREML")
packageVersion("ggmanh")
```

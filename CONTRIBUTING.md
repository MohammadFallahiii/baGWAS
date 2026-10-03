# Contributing to baGWAS

Thank you for contributing to baGWAS.

## Before opening an issue

Please check the README and existing issues first. For reproducible bug reports, include:

- operating system;
- R version;
- baGWAS version;
- PLINK version, when relevant;
- GEMMA version, when relevant;
- a minimal reproducible example when possible; and
- the exact error message and relevant `baGWAS.log` output.

Do not upload private genotype or phenotype data. Use synthetic, public, or anonymized examples instead.

## Pull requests

Please keep changes focused and document changes to user-facing behavior. Changes to statistical behavior should include a clear explanation and, where practical, a test.

Before submitting a pull request, run:

```r
testthat::test_local()
```

and, when possible:

```bash
R CMD check .
```

## Coding style

Use clear function names, explicit package namespaces, informative error messages, and comments for non-obvious logic. Avoid changing statistical behavior without documenting its methodological consequence.

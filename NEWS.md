# baGWAS 1.0.5

* Fixed external-command invocation so PLINK/GEMMA command names are not
  shell-quoted twice by `system2()`.
* Improved external-tool detection and diagnostics.

# baGWAS 1.0.4

* Tool detection tests each candidate (`plink2`, `plink1.9`, `plink`, `gemma`)
  and skips candidates that cannot be executed.
* `baGWAS_check_tools()` reports candidate commands and version information.

# baGWAS 1.0.3

* PLINK detection tries `plink2`, `plink1.9`, and `plink` in sequence.
* `freqGWAS` prefers PLINK 1.x for its stratified-frequency and IBS workflow.

# baGWAS 1.0.2

* External tools are invoked by command name resolved from `PATH`.

# baGWAS 1.0.1

* Fixed package metadata and maintainer information.
* PLINK and GEMMA are located at run time or through environment variables.
* Temporary files are isolated from the project working directory.
* `freqGWAS` uses the requested `group_col` and group-level matrix algebra.
* `freqGWAS` SNP scans use forked parallel workers on Unix-like systems.
* `genoGWAS` passes missing phenotypes to GEMMA as `NA`.
* Q-Q confidence intervals and lambda calculation were revised.

# baGWAS 1.0.0

* Initial release.

# Upload baGWAS 1.0.5 to GitHub

This guide describes the simplest browser-based upload workflow for a new GitHub repository.

## 1. Create the repository

Create a GitHub repository named:

```text
baGWAS
```

For the browser-upload method, leave the repository initialization options empty because this source tree already contains its own README, license, and Git-related metadata.

## 2. Extract this source tree

The package root must contain `DESCRIPTION` directly:

```text
baGWAS/
├── DESCRIPTION
├── NAMESPACE
├── README.md
└── R/
```

Do not create `baGWAS/baGWAS/DESCRIPTION`.

## 3. Upload through GitHub

On the repository page:

1. Open **Add file**.
2. Choose **Upload files**.
3. Drag the **contents of the package root** into the upload area.
4. Confirm that `DESCRIPTION` is visible at the repository root.
5. Use commit message:

```text
Initial release of baGWAS 1.0.5
```

6. Commit the upload.

## 4. Verify the repository

The repository root should contain:

```text
DESCRIPTION
NAMESPACE
README.md
LICENSE
NEWS.md
CITATION.cff
R/
man/
tests/
.github/
```

## 5. Test GitHub installation

After the repository is public:

```r
install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/baGWAS")
```

## 6. Create a release

For the first public version, create a Git tag and GitHub release named:

```text
v1.0.5
```

Suggested release title:

```text
baGWAS 1.0.5
```

Suggested release notes:

```text
baGWAS 1.0.5

- genotype-based GWAS workflow (genoGWAS)
- frequency-based GWAS workflow (freqGWAS)
- PLINK and GEMMA integration
- EMMREML-based frequency analysis
- Manhattan and Q-Q plots
- runtime logging and standardized result objects
- external-tool detection and diagnostics
- improved input validation and tool handling
```

## 7. Update the package later

Increment the version in `DESCRIPTION` and add a corresponding entry to `NEWS.md`. Commit and tag the update, for example:

```text
v1.0.6
```

Keep genotype, phenotype, result, and other potentially sensitive study files out of the public repository.

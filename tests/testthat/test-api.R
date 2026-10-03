ph <- data.frame(FID = 1:3, IID = 1:3, w = c(1, 2, 3))

test_that("invalid strategy is rejected", {
  expect_error(baGWAS(strategy = "foo", organism = "dog", pheno = ph,
                      trait = "w", geno = "x"))
})

test_that("unknown organism is rejected", {
  expect_error(baGWAS(strategy = "genoGWAS", organism = "unicorn",
                      pheno = ph, trait = "w", geno = "x"), "Unsupported")
})

test_that("missing trait is rejected", {
  expect_error(baGWAS(strategy = "genoGWAS", organism = "dog", pheno = ph,
                      trait = "nope", geno = "x"), "not found")
})

test_that("missing genotype files are reported", {
  expect_error(baGWAS(strategy = "genoGWAS", organism = "dog", pheno = ph,
                      trait = "w", geno = tempfile()), "Missing genotype")
})

test_that("genotype prefix is normalised", {
  expect_false(grepl("\\.bed$|/$", baGWAS:::.bg_normalize_geno("a/b.bed/")))
})

test_that("covariates default to intercept only", {
  fam <- data.frame(V1 = 1:3, V2 = 1:3)
  cm <- baGWAS:::.bg_cov_matrix(NULL, fam)
  expect_equal(cm$n_cov, 0L)
  expect_equal(ncol(cm$mat), 1L)
})

test_that("apply_feature_scaling minmax returns values in [0, 1]", {
  mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
  res <- apply_feature_scaling(mat, method = "minmax")
  expect_true(all(res >= 0 & res <= 1))
})

test_that("apply_feature_scaling none returns original matrix", {
  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  res <- apply_feature_scaling(mat, method = "none")
  expect_equal(res, mat)
})

test_that("log_transform returns log10(x + 1)", {
  mat <- matrix(c(0, 9, 99, 999), nrow = 2)
  res <- log_transform(mat)
  expect_equal(res, log10(mat + 1))
})

test_that("zscore_scale centers and scales columns", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
  res <- zscore_scale(mat)

  expect_equal(as.numeric(colMeans(res)), c(0, 0), tolerance = 1e-8)
  expect_equal(
    as.numeric(apply(res, 2, sd)),
    c(1, 1),
    tolerance = 1e-8
  )
})

test_that("minmax_scale rejects constant columns", {
  mat <- cbind(
    a = c(1, 2, 3),
    b = c(5, 5, 5)
  )
  expect_error(minmax_scale(mat), "constant columns")
})

test_that("remove_constant_features removes zero-variance columns", {
  mat <- cbind(
    a = c(1, 2, 3),
    b = c(5, 5, 5),
    c = c(2, 4, 6)
  )
  res <- remove_constant_features(mat)
  expect_equal(ncol(res), 2)
  expect_true(all(colnames(res) %in% c("a", "c")))
})

test_that("remove_constant_features errors if all columns are constant", {
  mat <- cbind(
    a = c(1, 1, 1),
    b = c(2, 2, 2)
  )
  expect_error(remove_constant_features(mat), "All columns were removed")
})

test_that("check_spectra_matrix rejects non-matrix input", {
  df <- data.frame(a = 1:3, b = 4:6)
  expect_error(check_spectra_matrix(df), "must be a matrix")
})

test_that("check_spectra_matrix rejects non-numeric matrix", {
  mat <- matrix(c("a", "b", "c", "d"), nrow = 2)
  expect_error(check_spectra_matrix(mat), "numeric")
})

test_that("run_clustering returns cluster labels", {
  set.seed(1)
  emb <- data.frame(UMAP1 = rnorm(20), UMAP2 = rnorm(20))
  res <- run_clustering(emb, centers = 3)

  expect_equal(length(res$cluster), 20)
  expect_true(is.factor(res$cluster))
})

test_that("run_clustering rejects embedding with fewer than 2 columns", {
  emb <- data.frame(UMAP1 = rnorm(10))
  expect_error(run_clustering(emb), "at least 2 columns")
})

test_that("build_cluster_dataframe returns expected columns", {
  emb <- data.frame(UMAP1 = rnorm(5), UMAP2 = rnorm(5))
  cluster <- factor(c(1, 1, 2, 2, 3))
  pixel_info <- data.frame(pixel_ID = paste0("p", 1:5), x = 1:5)

  res <- build_cluster_dataframe(emb, cluster, pixel_info)

  expect_equal(nrow(res), 5)
  expect_true(all(c("pixel_ID", "UMAP1", "UMAP2", "cluster", "x") %in% colnames(res)))
})

test_that("build_cluster_dataframe rejects missing pixel_ID", {
  emb <- data.frame(UMAP1 = rnorm(5), UMAP2 = rnorm(5))
  cluster <- factor(c(1, 1, 2, 2, 3))
  pixel_info <- data.frame(x = 1:5)

  expect_error(
    build_cluster_dataframe(emb, cluster, pixel_info),
    "pixel_ID"
  )
})

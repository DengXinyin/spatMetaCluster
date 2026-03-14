test_that("apply_feature_scaling returns values in [0, 1]", {
  mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3)
  res <- apply_feature_scaling(mat, method = "minmax")
  expect_true(all(res >= 0 & res <= 1))
})

test_that("remove_constant_features removes zero-variance columns", {
  mat <- cbind(
    a = c(1, 2, 3),
    b = c(5, 5, 5),
    c = c(2, 4, 6)
  )
  res <- remove_constant_features(mat)
  expect_equal(ncol(res), 2)
})

test_that("run_clustering returns cluster labels", {
  emb <- data.frame(UMAP1 = rnorm(20), UMAP2 = rnorm(20))
  res <- run_clustering(emb, centers = 3)
  expect_equal(length(res$cluster), 20)
})

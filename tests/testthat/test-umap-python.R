test_that("check_umap_python_env works when python modules are available", {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("umap"))
  skip_if_not(reticulate::py_module_available("numpy"))

  expect_invisible(check_umap_python_env())
})

test_that("run_umap_py returns embedding with requested dimensions", {
  skip_if_not_installed("reticulate")
  skip_if_not(reticulate::py_module_available("umap"))
  skip_if_not(reticulate::py_module_available("numpy"))

  set.seed(1)
  mat <- matrix(runif(100), nrow = 10)

  res <- run_umap_py(
    mat,
    metric = "cosine",
    n_components = 2L,
    random_state = 42L,
    n_jobs = 1L,
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 10)
  expect_equal(ncol(res), 2)
  expect_equal(colnames(res), c("UMAP1", "UMAP2"))
})

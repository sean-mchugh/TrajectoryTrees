test_that("batch files round trip with filters recorded in the manifest", {
  constructor <- public_function("make_trajectory_batches")
  loader <- public_function("load_trajectory_batches")

  expect_true(is.function(constructor))
  expect_true(is.function(loader))
  if (!is.function(constructor) || !is.function(loader)) {
    return(invisible())
  }

  output_dir <- file.path(tempdir(), paste0("trajectory-batch-", Sys.getpid()))
  tree <- small_simmap()
  manifest <- constructor(
    list(a = tree, b = tree, c = tree),
    output_dir = output_dir,
    batch_size = 2L,
    include_tt = FALSE,
    progress = FALSE
  )
  loaded <- loader(manifest)

  expect_s3_class(manifest, "trajectory_batch_manifest")
  expect_length(manifest$files, 2L)
  expect_named(loaded, c("a", "b", "c"))
  expect_s3_class(loaded$a, "trajectory_tree")
})

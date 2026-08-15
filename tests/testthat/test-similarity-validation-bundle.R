test_that("the preserved 500-tree similarity audit is self-contained", {
  audit_directory <- testthat::test_path(
    "..",
    "validation",
    "scenariomats_v3_v21_500"
  )
  manifest_file <- file.path(audit_directory, "selected_inventory.csv")

  expect_true(dir.exists(audit_directory))
  expect_true(file.exists(manifest_file))
  if (!file.exists(manifest_file)) {
    return(invisible())
  }

  inventory <- utils::read.csv(manifest_file, stringsAsFactors = FALSE)
  expect_equal(nrow(inventory), 500L)
  expect_identical(inventory$selection_order, seq_len(500L))
  expect_identical(anyDuplicated(inventory$saved_tree_id), 0L)
  expect_true(all(c(
    "tree_path",
    "tip_count",
    "state_count",
    "is_ultrametric"
  ) %in% names(inventory)))

  tree_files <- file.path(
    audit_directory,
    "trees",
    basename(inventory$tree_path)
  )
  expect_identical(
    anyDuplicated(normalizePath(tree_files, mustWork = FALSE)),
    0L
  )
  expect_true(all(file.exists(tree_files)))
})

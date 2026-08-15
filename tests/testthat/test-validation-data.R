test_that("validation inventory can be filtered by size and label", {
  loader <- public_function("load_validation_trees")

  expect_true(is.function(loader))
  if (!is.function(loader)) {
    return(invisible())
  }

  selected <- loader(
    min_tips = 5L,
    max_tips = 10L,
    exclude_label = "^generated_",
    limit = 3L,
    load = FALSE
  )

  expect_s3_class(selected, "trajectory_tree_inventory")
  expect_lte(nrow(selected), 3L)
  expect_true(all(selected$n_tip >= 5L & selected$n_tip <= 10L))
  expect_false(any(grepl("^generated_", selected$case_id)))
})

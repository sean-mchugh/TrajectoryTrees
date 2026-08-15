trajectory_v35_custom_attribute_paths <- function(
    object,
    path = character()) {
  custom_names <- grep(
    "^pst_v35_",
    names(attributes(object)),
    value = TRUE
  )
  current <- vapply(
    custom_names,
    function(name) paste(c(path, name), collapse = "::"),
    character(1)
  )
  nested <- character()
  if (is.list(object)) {
    object_names <- names(object)
    for (index in seq_along(object)) {
      component <- if (
          !is.null(object_names) &&
          nzchar(object_names[[index]])
      ) {
        object_names[[index]]
      } else {
        paste0("[[", index, "]]")
      }
      nested <- c(
        nested,
        trajectory_v35_custom_attribute_paths(
          object[[index]],
          c(path, component)
        )
      )
    }
  }
  c(current, nested)
}

test_that("one simmap becomes a stable trajectory_tree object", {
  constructor <- public_function("make_trajectory_objects")

  expect_true(is.function(constructor))
  if (!is.function(constructor)) {
    return(invisible())
  }

  result <- constructor(small_simmap(), include_tt = FALSE)
  expect_s3_class(result, "trajectory_tree")
  expect_true(all(c(
    "phylo", "scenario", "trans", "ss", "scenario_mats", "root_policy"
  ) %in% names(result)))
  expect_identical(
    unname(sort(trajectory_v35_custom_attribute_paths(result))),
    sort(c(
      "pst_v35_path_lookup",
      "pst_v35_time_tolerance",
      "scenario::pst_v35_stable_scenario_ids"
    ))
  )
})

test_that("internal V35 attributes require an explicit request", {
  constructor <- public_function("make_trajectory_objects")

  expect_true(is.function(constructor))
  if (!is.function(constructor)) {
    return(invisible())
  }

  result <- constructor(
    small_simmap(),
    include_tt = FALSE,
    keep_internal_attributes = TRUE
  )
  expect_true(all(c(
    "pst_v35_cpp_traversal_summary",
    "pst_v35_traversal_scenario_ledger",
    "pst_v35_stage_timings"
  ) %in% names(attributes(result))))
})

test_that("input list shape is preserved", {
  constructor <- public_function("make_trajectory_objects")

  expect_true(is.function(constructor))
  if (!is.function(constructor)) {
    return(invisible())
  }

  tree <- small_simmap()
  flat <- constructor(list(first = tree, second = tree), include_tt = FALSE)
  nested <- constructor(
    list(group_a = list(first = tree), group_b = list(second = tree)),
    include_tt = FALSE
  )

  expect_named(flat, c("first", "second"))
  expect_s3_class(flat[[1]], "trajectory_tree")
  expect_named(nested, c("group_a", "group_b"))
  expect_s3_class(nested$group_a$first, "trajectory_tree")
})

test_that("selected fields keep the stable class", {
  constructor <- public_function("make_trajectory_objects")
  selector <- public_function("select_trajectory_fields")

  expect_true(is.function(constructor))
  expect_true(is.function(selector))
  if (!is.function(constructor) || !is.function(selector)) {
    return(invisible())
  }

  result <- constructor(small_simmap(), include_tt = FALSE)
  selected <- selector(result, include = c("phylo", "scenario", "trans"))

  expect_s3_class(selected, "trajectory_tree")
  expect_named(selected, c("phylo", "scenario", "trans"))
})

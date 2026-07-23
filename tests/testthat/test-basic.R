test_that("default loci order is available", {
  loci <- default_loci_order()

  expect_type(loci, "character")
  expect_true(length(loci) > 0)
  expect_equal(loci[1], "rs10495407")
  expect_true("rs1490413" %in% loci)
})


test_that("default thresholds by locus are available", {
  thresholds <- default_thresholds_by_locus()

  expect_type(thresholds, "double")
  expect_true("rs729172" %in% names(thresholds))
  expect_equal(unname(thresholds["rs729172"]), 0.1)
  expect_equal(unname(thresholds["rs338882"]), 0.2)
})


test_that("default review balance margins are available", {
  expect_equal(default_review_balance_lower_margin(), 0.15)
  expect_equal(default_review_balance_upper_margin(), 0.15)
})


test_that("bundled haplotype database can be loaded", {
  db <- load_haplotype_database()

  expect_s3_class(db, "data.frame")
  expect_true(nrow(db) > 0)

  expected_columns <- c(
    "TargetSNP",
    "TargetSNP_allele",
    "Haplotype",
    "Sequence",
    "Polymorphic_sites",
    "Position_GRCh38"
  )

  expect_true(all(expected_columns %in% names(db)))
})


test_that("bundled haplotype database has no duplicated keys", {
  db <- load_haplotype_database()

  duplicated_keys <- db |>
    dplyr::mutate(sequence_key = normalize_sequence(.data$Sequence)) |>
    dplyr::count(.data$TargetSNP, .data$sequence_key, name = "n") |>
    dplyr::filter(.data$n > 1)

  expect_equal(nrow(duplicated_keys), 0)
})


test_that("normalize_locus standardizes rs identifiers", {
  expect_equal(normalize_locus("rs 1490413"), "rs1490413")
  expect_equal(normalize_locus("RS1490413"), "rs1490413")
})

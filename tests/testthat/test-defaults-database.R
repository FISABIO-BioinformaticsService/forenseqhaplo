test_that("default loci order is available", {
  loci <- default_loci_order()

  expect_type(loci, "character")
  expect_true(length(loci) > 0)
  expect_equal(loci[1], "rs10495407")
  expect_true("rs1490413" %in% loci)
  expect_equal(length(unique(loci)), length(loci))
})


test_that("default thresholds are available", {
  thresholds <- default_thresholds_by_locus()

  expect_type(thresholds, "double")
  expect_true("rs729172" %in% names(thresholds))
  expect_true("rs338882" %in% names(thresholds))
  expect_equal(unname(thresholds["rs729172"]), 0.1)
  expect_equal(unname(thresholds["rs338882"]), 0.2)

  expect_equal(default_heterozygote_threshold(), 0.3)
  expect_equal(default_min_homozygote_reads(), 30)
  expect_equal(default_min_heterozygous_reads(), 11)
  expect_equal(default_extra_signal_threshold(), 0.25)
  expect_equal(default_review_balance_lower_margin(), 0.15)
  expect_equal(default_review_balance_upper_margin(), 0.15)
  expect_equal(default_review_low_homozygote_multiplier(), 1.5)
})


test_that("bundled haplotype database can be loaded and validated", {
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
  expect_true(validate_haplotype_database())
})

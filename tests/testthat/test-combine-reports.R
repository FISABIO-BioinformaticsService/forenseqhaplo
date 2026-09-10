test_that("replicate review reports can be combined", {
  input_dir <- make_synthetic_report_dir()
  output_dir <- tempfile("forenseqhaplo_combine_output_")

  process_forenseq_folder(
    input_dir = input_dir,
    output_dir = output_dir,
    write_parameters_xlsx = FALSE,
    write_log = FALSE
  )

  combined <- combine_reports(
    reports_dir = output_dir,
    output_txt = "combine_haplotypes.txt"
  )

  expect_s3_class(combined$consensus_txt, "data.frame")
  expect_equal(nrow(combined$consensus_txt), 1)
  expect_true("sample_id" %in% names(combined$consensus_txt))
  expect_true(file.exists(combined$output_txt))
  expect_s3_class(combined$conflicts, "data.frame")
})

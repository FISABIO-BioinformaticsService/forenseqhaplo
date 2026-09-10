test_that("target TXT is written when requested", {
  input_dir <- make_synthetic_report_dir()
  output_dir <- tempfile("forenseqhaplo_target_output_")

  result <- process_forenseq_folder(
    input_dir = input_dir,
    output_dir = output_dir,
    write_target_txt = TRUE,
    write_parameters_xlsx = FALSE,
    write_log = FALSE
  )

  expect_true(file.exists(result$output_txt))
  expect_true(file.exists(result$output_target_txt))

  target <- read.delim(
    result$output_target_txt,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expect_true("sample_id" %in% names(target))
  expect_true(nrow(target) > 0)
})

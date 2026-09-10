test_that("DVI mode writes relationship columns in the expected order", {
  input_dir <- make_synthetic_report_dir()
  output_dir <- tempfile("forenseqhaplo_dvi_output_")
  dvi_file <- get_synthetic_dvi_file()

  result <- process_forenseq_folder(
    input_dir = input_dir,
    output_dir = output_dir,
    dvi = TRUE,
    dvi_relationships = dvi_file,
    write_parameters_xlsx = FALSE,
    write_log = FALSE
  )

  expect_true(file.exists(result$output_txt))

  txt <- read.delim(
    result$output_txt,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  expect_equal(names(txt)[1:3], c("family_id", "relationship", "sample_id"))
  expect_true(all(!is.na(txt$family_id)))
  expect_true(all(!is.na(txt$relationship)))
  expect_true(all(!is.na(txt$sample_id)))
})

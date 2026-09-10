test_that("folder processing works with bundled synthetic reports", {
  input_dir <- make_synthetic_report_dir()
  output_dir <- tempfile("forenseqhaplo_output_")

  result <- process_forenseq_folder(
    input_dir = input_dir,
    output_dir = output_dir,
    write_parameters_xlsx = TRUE,
    write_log = TRUE
  )

  expect_s3_class(result$familias_txt, "data.frame")
  expect_true(nrow(result$familias_txt) > 0)
  expect_true("sample_id" %in% names(result$familias_txt))

  expect_true(file.exists(result$output_txt))
  expect_true(file.exists(result$output_parameters_xlsx))
  expect_true(file.exists(file.path(output_dir, "processing_log.tsv")))

  review_files <- list.files(output_dir, pattern = "^review_.*\\.xlsx$", full.names = TRUE)
  expect_true(length(review_files) > 0)

  sheets <- readxl::excel_sheets(review_files[[1]])
  expect_true(all(c("Final_calls", "Review", "All_sequences") %in% sheets))
})


test_that("single sample processing creates a review workbook", {
  input_dir <- make_synthetic_report_dir()
  report_file <- list.files(
    input_dir,
    pattern = "^synthetic_forenseq_report_.*\\.xlsx$",
    full.names = TRUE
  )[[1]]

  output_xlsx <- file.path(tempfile("forenseqhaplo_single_"), "review_single.xlsx")
  dir.create(dirname(output_xlsx), recursive = TRUE, showWarnings = FALSE)

  result <- process_forenseq_sample(
    report_xlsx = report_file,
    output_xlsx = output_xlsx,
    write_parameters_xlsx = FALSE
  )

  expect_true(file.exists(output_xlsx))
  expect_s3_class(result$familias_row, "data.frame")
  expect_true(nrow(result$familias_row) == 1)

  final_calls <- readxl::read_excel(
    output_xlsx,
    sheet = "Final_calls",
    col_types = "text",
    .name_repair = "unique_quiet"
  )

  expect_true(all(c("sample_id", "TargetSNP", "zygosity") %in% names(final_calls)))
})

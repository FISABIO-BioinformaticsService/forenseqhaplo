make_synthetic_report_dir <- function() {
  extdata_dir <- system.file("extdata", package = "forenseqhaplo", mustWork = TRUE)

  report_files <- list.files(
    extdata_dir,
    pattern = "^synthetic_forenseq_report_.*\\.xlsx$",
    full.names = TRUE
  )

  testthat::skip_if(
    length(report_files) == 0,
    "Synthetic ForenSeq reports are not installed."
  )

  output_dir <- tempfile("forenseqhaplo_input_")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  file.copy(
    from = report_files,
    to = file.path(output_dir, basename(report_files)),
    overwrite = TRUE
  )

  output_dir
}

get_synthetic_dvi_file <- function() {
  dvi_file <- system.file(
    "extdata",
    "synthetic_dvi_relationships.xlsx",
    package = "forenseqhaplo",
    mustWork = FALSE
  )

  testthat::skip_if(
    !nzchar(dvi_file) || !file.exists(dvi_file),
    "Synthetic DVI relationships file is not installed."
  )

  dvi_file
}

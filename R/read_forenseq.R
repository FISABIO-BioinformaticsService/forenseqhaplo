find_forenseq_header_row <- function(path, sheet) {
  raw <- readxl::read_excel(
    path = path,
    sheet = sheet,
    col_names = FALSE,
    col_types = "text"
  )

  for (i in seq_len(nrow(raw))) {
    values <- unlist(raw[i, ], use.names = FALSE)
    values_norm <- normalize_colname(values)

    has_locus <- any(
      values_norm %in% c(
        "isnp locus",
        "isnp & variant reference snp",
        "isnp and variant reference snp"
      ),
      na.rm = TRUE
    )

    has_detected <- any(
      values_norm %in% c("detected bases", "detected base"),
      na.rm = TRUE
    )

    has_read <- any(
      values_norm %in% c("read", "reads", "read count"),
      na.rm = TRUE
    )

    has_sequence <- any(values_norm == "sequence", na.rm = TRUE)

    if (has_locus && has_detected && has_read && has_sequence) {
      return(i)
    }
  }

  stop(
    "Could not find the table header row in sheet '",
    sheet,
    "'. Expected columns: iSNP Locus, Detected Bases, Read and Sequence."
  )
}


#' Read sample ID from a ForenSeq Flanking Region Report
#'
#' Reads the sample identifier from cell B5 of the selected ForenSeq report
#' sheet. In the expected ForenSeq report layout, this cell corresponds to
#' the Sample field.
#'
#' @param report_xlsx Path to the ForenSeq Flanking Region Report Excel file.
#' @param report_sheet Sheet containing the iSNP coverage table.
#'
#' @return A character scalar with the sample ID.
#' @noRd
read_forenseq_sample_id <- function(
    report_xlsx,
    report_sheet = "iSNP Coverage"
) {
  value <- readxl::read_excel(
    path = report_xlsx,
    sheet = report_sheet,
    range = "B5",
    col_names = FALSE,
    col_types = "text"
  )[[1]][1]

  value <- clean_text(value)

  if (is.na(value) || value == "") {
    warning("Could not read sample_id from cell B5.")
    return(NA_character_)
  }

  value
}


#' Read iSNP coverage data from a ForenSeq report
#'
#' Reads and standardizes the relevant columns from the `iSNP Coverage` sheet
#' of a ForenSeq Flanking Region Report.
#'
#' @param report_xlsx Path to the ForenSeq Flanking Region Report Excel file.
#' @param report_sheet Sheet containing the iSNP coverage table.
#'
#' @return A tibble with standardized columns.
#' @noRd
read_forenseq_isnp_coverage <- function(
    report_xlsx,
    report_sheet = "iSNP Coverage"
) {
  header_row <- find_forenseq_header_row(
    path = report_xlsx,
    sheet = report_sheet
  )

  report <- readxl::read_excel(
    path = report_xlsx,
    sheet = report_sheet,
    skip = header_row - 1,
    col_names = TRUE,
    col_types = "text"
  )

  col_locus <- get_column(
    report,
    c(
      "iSNP Locus",
      "iSNP & Variant Reference SNP",
      "iSNP and Variant Reference SNP"
    )
  )

  col_detected <- get_column(
    report,
    c("Detected Bases", "Detected Base")
  )

  col_read <- get_column(
    report,
    c("Read", "Reads", "Read Count")
  )

  col_sequence <- get_column(
    report,
    c("Sequence")
  )

  out <- tibble::tibble(
    row_id = seq_len(nrow(report)),
    isnp_locus_original = as.character(report[[col_locus]]),
    `iSNP Locus` = normalize_locus(report[[col_locus]]),
    `Detected Bases` = normalize_detected_bases(report[[col_detected]]),
    Read = parse_read_count(report[[col_read]]),
    Sequence = as.character(report[[col_sequence]]),
    sequence_key = normalize_sequence(report[[col_sequence]])
  )

  out <- out[
    !is.na(out$`iSNP Locus`) &
      out$`iSNP Locus` != "" &
      !is.na(out$sequence_key) &
      out$sequence_key != "",
  ]

  if (nrow(out) == 0) {
    stop("The ForenSeq report is empty after cleaning.")
  }

  out
}

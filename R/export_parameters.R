resolve_output_parameters_xlsx_path <- function(
    write_parameters_xlsx = TRUE,
    output_parameters_xlsx = NULL,
    output_dir = NULL,
    output_xlsx = NULL
) {
  if (!isTRUE(write_parameters_xlsx)) {
    return(NULL)
  }

  if (!is.null(output_parameters_xlsx) && output_parameters_xlsx != "") {
    if (!is.null(output_dir) && dirname(output_parameters_xlsx) == ".") {
      return(file.path(output_dir, output_parameters_xlsx))
    }

    return(output_parameters_xlsx)
  }

  if (!is.null(output_dir)) {
    return(file.path(output_dir, "parameters.xlsx"))
  }

  if (!is.null(output_xlsx) && output_xlsx != "" && dirname(output_xlsx) != ".") {
    return(file.path(dirname(output_xlsx), "parameters.xlsx"))
  }

  "parameters.xlsx"
}


format_parameter_value <- function(x) {
  if (length(x) == 0 || is.null(x)) {
    return(NA_character_)
  }

  if (length(x) > 1) {
    x <- paste(as.character(x), collapse = ",")
  }

  x <- as.character(x)

  if (is.na(x)) {
    return(NA_character_)
  }

  x
}


create_processing_parameters_table <- function(
    processing_mode,
    input_dir = NA_character_,
    report_xlsx = NA_character_,
    output_dir = NA_character_,
    output_xlsx = NA_character_,
    output_txt = NA_character_,
    write_target_txt = FALSE,
    output_target_txt = NA_character_,
    haplotype_db = NULL,
    report_sheet = "iSNP Coverage",
    db_sheet = 1,
    loci_order = default_loci_order(),
    heterozygote_threshold = default_heterozygote_threshold(),
    extra_signal_threshold = default_extra_signal_threshold(),
    thresholds_by_locus = default_thresholds_by_locus(),
    min_homozygote_reads = default_min_homozygote_reads(),
    min_heterozygous_reads = default_min_heterozygous_reads(),
    review_balance_lower_margin = default_review_balance_lower_margin(),
    review_balance_upper_margin = default_review_balance_upper_margin(),
    review_low_homozygote_multiplier = default_review_low_homozygote_multiplier(),
    write_log = NA,
    dvi = FALSE
) {
  thresholds_by_locus <- normalize_thresholds_by_locus(thresholds_by_locus)
  loci_order <- normalize_locus(loci_order)
  loci_order <- loci_order[!is.na(loci_order) & loci_order != ""]

  scalar_parameters <- tibble::tibble(
    section = c(
      "processing",
      "processing",
      "processing",
      "input",
      "input",
      "input",
      "output",
      "output",
      "output",
      "output",
      "output",
      "database",
      "database",
      "panel",
      "panel",
      "calling",
      "calling",
      "calling",
      "calling",
      "review",
      "review",
      "review",
      "logging",
      "dvi"
    ),
    parameter = c(
      "package",
      "processing_mode",
      "created_at",
      "input_dir",
      "report_xlsx",
      "report_sheet",
      "output_dir",
      "output_xlsx",
      "output_txt",
      "write_target_txt",
      "output_target_txt",
      "haplotype_db",
      "db_sheet",
      "n_loci_order",
      "loci_order",
      "heterozygote_threshold",
      "extra_signal_threshold",
      "min_homozygote_reads",
      "min_heterozygous_reads",
      "review_balance_lower_margin",
      "review_balance_upper_margin",
      "review_low_homozygote_multiplier",
      "write_log",
      "dvi"
    ),
    value = c(
      "forenseqhaplo",
      format_parameter_value(processing_mode),
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      format_parameter_value(input_dir),
      format_parameter_value(report_xlsx),
      format_parameter_value(report_sheet),
      format_parameter_value(output_dir),
      format_parameter_value(output_xlsx),
      format_parameter_value(output_txt),
      format_parameter_value(write_target_txt),
      if (isTRUE(write_target_txt)) format_parameter_value(output_target_txt) else "not_written",
      if (is.null(haplotype_db)) "internal" else format_parameter_value(haplotype_db),
      format_parameter_value(db_sheet),
      format_parameter_value(length(loci_order)),
      paste(loci_order, collapse = ","),
      format_parameter_value(heterozygote_threshold),
      format_parameter_value(extra_signal_threshold),
      format_parameter_value(min_homozygote_reads),
      format_parameter_value(min_heterozygous_reads),
      format_parameter_value(review_balance_lower_margin),
      format_parameter_value(review_balance_upper_margin),
      format_parameter_value(review_low_homozygote_multiplier),
      format_parameter_value(write_log),
      format_parameter_value(dvi)
    )
  )

  if (length(thresholds_by_locus) > 0) {
    locus_thresholds <- tibble::tibble(
      section = "locus_specific_heterozygote_threshold",
      parameter = names(thresholds_by_locus),
      value = as.character(unname(thresholds_by_locus))
    )

    dplyr::bind_rows(scalar_parameters, locus_thresholds)
  } else {
    scalar_parameters
  }
}


export_processing_parameters_xlsx <- function(parameters, output_parameters_xlsx) {
  writexl::write_xlsx(
    list(Parameters = parameters),
    path = output_parameters_xlsx
  )

  invisible(parameters)
}

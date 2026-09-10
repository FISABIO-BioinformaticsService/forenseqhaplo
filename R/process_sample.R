process_forenseq_sample_core <- function(
    report_xlsx,
    haplotype_data,
    output_xlsx,
    output_txt = NULL,
    write_target_txt = FALSE,
    output_target_txt = NULL,
    haplotype_db = NULL,
    report_sheet = "iSNP Coverage",
    db_sheet = 1,
    loci_order = default_loci_order(),
    original_sample_id = NULL,
    txt_sample_id = NULL,
    heterozygote_threshold = default_heterozygote_threshold(),
    extra_signal_threshold = default_extra_signal_threshold(),
    thresholds_by_locus = default_thresholds_by_locus(),
    min_homozygote_reads = default_min_homozygote_reads(),
    min_heterozygous_reads = default_min_heterozygous_reads(),
    review_balance_lower_margin = default_review_balance_lower_margin(),
    review_balance_upper_margin = default_review_balance_upper_margin(),
    review_low_homozygote_multiplier = default_review_low_homozygote_multiplier(),
    write_parameters_xlsx = TRUE,
    output_parameters_xlsx = NULL,
    dvi = FALSE,
    dvi_data = NULL
) {
  if (is.null(original_sample_id)) {
    original_sample_id <- read_forenseq_sample_id(
      report_xlsx = report_xlsx,
      report_sheet = report_sheet
    )
  }

  if (is.na(original_sample_id) || original_sample_id == "") {
    original_sample_id <- tools::file_path_sans_ext(basename(report_xlsx))
  }

  if (is.null(txt_sample_id) || is.na(txt_sample_id) || txt_sample_id == "") {
    txt_sample_id <- original_sample_id
  }

  report_data <- read_forenseq_isnp_coverage(
    report_xlsx = report_xlsx,
    report_sheet = report_sheet
  )

  classified <- classify_all_loci(
    report_data = report_data,
    loci_order = loci_order,
    heterozygote_threshold = heterozygote_threshold,
    extra_signal_threshold = extra_signal_threshold,
    thresholds_by_locus = thresholds_by_locus,
    min_homozygote_reads = min_homozygote_reads,
    min_heterozygous_reads = min_heterozygous_reads
  )

  selected_evidence <- annotate_with_haplotype_db(
    df = classified$selected_rows,
    haplotype_db = haplotype_data
  )

  all_sequences_annotated <- annotate_with_haplotype_db(
    df = classified$all_sequences,
    haplotype_db = haplotype_data
  )

  familias_outputs <- create_familias_outputs(
    selected_evidence = selected_evidence,
    locus_summary = classified$locus_summary,
    all_sequences_annotated = all_sequences_annotated,
    review_balance_lower_margin = review_balance_lower_margin,
    review_balance_upper_margin = review_balance_upper_margin,
    review_low_homozygote_multiplier = review_low_homozygote_multiplier
  )

  familias_input <- familias_outputs$familias_input
  familias_target_input <- familias_outputs$familias_target_input
  final_calls <- familias_outputs$final_calls
  review <- familias_outputs$review
  all_sequences <- familias_outputs$all_sequences

  dvi_info <- NULL

  if (isTRUE(dvi)) {
    dvi_info <- get_dvi_relationship_for_sample(
      dvi_data = dvi_data,
      sample_id = txt_sample_id
    )
  }

  if (!isTRUE(write_target_txt) && !is.null(output_target_txt)) {
    warning("`output_target_txt` was provided but `write_target_txt = FALSE`; it will be ignored.")
  }

  output_target_txt_path <- resolve_output_target_txt_path(
    write_target_txt = write_target_txt,
    output_target_txt = output_target_txt,
    output_dir = NULL,
    output_txt_path = output_txt,
    output_xlsx = output_xlsx
  )

  output_parameters_xlsx_path <- resolve_output_parameters_xlsx_path(
    write_parameters_xlsx = write_parameters_xlsx,
    output_parameters_xlsx = output_parameters_xlsx,
    output_dir = NULL,
    output_xlsx = output_xlsx
  )

  familias_row <- create_familias_row(
    familias_input = familias_input,
    sample_id = txt_sample_id,
    loci_order = loci_order,
    include_missing_loci = TRUE,
    missing_value = "",
    dvi_info = dvi_info
  )

  familias_target_row <- create_familias_row(
    familias_input = familias_target_input,
    sample_id = txt_sample_id,
    loci_order = loci_order,
    include_missing_loci = TRUE,
    missing_value = "",
    dvi_info = dvi_info
  )

  if (!is.null(output_txt)) {
    export_familias_txt(
      familias_rows = familias_row,
      output_txt = output_txt
    )
  }

  if (!is.null(output_target_txt_path)) {
    export_familias_txt(
      familias_rows = familias_target_row,
      output_txt = output_target_txt_path
    )
  }

  final_calls_excel <- final_calls |>
    add_sample_metadata(
      original_sample_id = original_sample_id,
      txt_sample_id = txt_sample_id,
      source_file = report_xlsx
    ) |>
    drop_excel_only_columns() |>
    order_final_calls_columns()

  review_excel <- review |>
    add_sample_metadata(
      original_sample_id = original_sample_id,
      txt_sample_id = txt_sample_id,
      source_file = report_xlsx
    ) |>
    drop_excel_only_columns() |>
    order_review_columns()

  all_sequences_excel <- all_sequences |>
    add_sample_metadata(
      original_sample_id = original_sample_id,
      txt_sample_id = txt_sample_id,
      source_file = report_xlsx
    ) |>
    drop_excel_only_columns() |>
    order_all_sequences_columns()

  writexl::write_xlsx(
    list(
      Final_calls = final_calls_excel,
      Review = review_excel,
      All_sequences = all_sequences_excel
    ),
    path = output_xlsx
  )

  processing_parameters <- NULL

  if (isTRUE(write_parameters_xlsx)) {
    processing_parameters <- create_processing_parameters_table(
      processing_mode = "sample",
      report_xlsx = report_xlsx,
      output_xlsx = output_xlsx,
      output_txt = output_txt,
      write_target_txt = write_target_txt,
      output_target_txt = output_target_txt_path,
      haplotype_db = haplotype_db,
      report_sheet = report_sheet,
      db_sheet = db_sheet,
      loci_order = loci_order,
      heterozygote_threshold = heterozygote_threshold,
      extra_signal_threshold = extra_signal_threshold,
      thresholds_by_locus = thresholds_by_locus,
      min_homozygote_reads = min_homozygote_reads,
      min_heterozygous_reads = min_heterozygous_reads,
      review_balance_lower_margin = review_balance_lower_margin,
      review_balance_upper_margin = review_balance_upper_margin,
      review_low_homozygote_multiplier = review_low_homozygote_multiplier,
      write_log = NA,
      dvi = dvi
    )

    export_processing_parameters_xlsx(
      parameters = processing_parameters,
      output_parameters_xlsx = output_parameters_xlsx_path
    )
  }

  process_summary <- tibble::tibble(
    source_file = basename(report_xlsx),
    original_sample_id = original_sample_id,
    txt_sample_id = txt_sample_id,
    status = "OK",
    error_message = NA_character_,
    output_xlsx = output_xlsx,
    n_loci_observed_in_report = dplyr::n_distinct(report_data$`iSNP Locus`),
    n_loci_processed = nrow(final_calls),
    n_loci_with_haplotypes_txt = sum(final_calls$included_in_txt, na.rm = TRUE),
    n_loci_with_target_txt = sum(final_calls$included_in_target_txt, na.rm = TRUE),
    n_homozygous_loci = sum(final_calls$zygosity == "homozygous", na.rm = TRUE),
    n_heterozygous_loci = sum(final_calls$zygosity == "heterozygous", na.rm = TRUE),
    n_no_call_loci = sum(final_calls$zygosity == "no_call", na.rm = TRUE),
    n_review_rows = nrow(review),
    n_database_not_found_loci = sum(
      final_calls$txt_action == "locus_excluded_from_txt_not_found_in_database",
      na.rm = TRUE
    ),
    dvi = isTRUE(dvi)
  )

  invisible(
    list(
      original_sample_id = original_sample_id,
      txt_sample_id = txt_sample_id,
      familias_row = familias_row,
      familias_target_row = familias_target_row,
      familias_input = familias_input,
      familias_target_input = familias_target_input,
      final_calls = final_calls,
      review = review,
      all_sequences = all_sequences,
      process_summary = process_summary,
      output_xlsx = output_xlsx,
      output_txt = output_txt,
      output_target_txt = output_target_txt_path,
      processing_parameters = processing_parameters,
      output_parameters_xlsx = output_parameters_xlsx_path
    )
  )
}


#' Process one ForenSeq Flanking Region Report
#'
#' Processes a single ForenSeq Flanking Region Report using sequence-based
#' haplotype calling.
#'
#' @param report_xlsx Path to the ForenSeq Flanking Region Report Excel file.
#' @param output_xlsx Path to the Excel review file to create.
#' @param output_txt Optional path to a Familias TXT file for this sample.
#' @param write_target_txt Logical. If `TRUE`, also writes a TargetSNP_allele TXT file.
#'   Defaults to `FALSE`.
#' @param output_target_txt Optional path to the TargetSNP_allele TXT file.
#'   Used only when `write_target_txt = TRUE`.
#' @param write_parameters_xlsx Logical. If `TRUE`, writes an Excel file with the
#'   processing parameters used in the run. Defaults to `TRUE`.
#' @param output_parameters_xlsx Optional path to the processing parameters Excel file.
#'   If `NULL`, `parameters.xlsx` is written next to `output_xlsx`.
#' @param haplotype_db Optional path to a user-provided haplotype database Excel file.
#' @param report_sheet Sheet containing the iSNP coverage table.
#' @param db_sheet Excel sheet to read when `haplotype_db` is provided.
#' @param loci_order Character vector with the loci order used in the TXT.
#' @param heterozygote_threshold General read-ratio threshold.
#' @param extra_signal_threshold Threshold used to flag relevant third sequences.
#' @param thresholds_by_locus Named numeric vector with locus-specific thresholds.
#' @param min_homozygote_reads Minimum reads required to duplicate one sequence as homozygous.
#' @param min_heterozygous_reads Minimum reads required for the second sequence
#'   to be accepted as the second allele in a heterozygous call.
#' @param review_balance_lower_margin Proportional lower margin below the allele
#'   balance threshold used to flag borderline calls for manual review.
#' @param review_balance_upper_margin Proportional upper margin above the allele
#'   balance threshold used to flag borderline calls for manual review.
#' @param review_low_homozygote_multiplier Multiplier applied to
#'   `min_homozygote_reads` to flag homozygous calls close to the minimum read
#'   threshold for manual review.
#' @param duplicates_output Optional path where duplicated haplotype database
#'   records are written if duplicated `TargetSNP + Sequence` keys are found.
#'   Mainly useful when using a custom haplotype database.
#' @param dvi Logical. If `TRUE`, DVI relationship columns are added.
#' @param dvi_relationships Path to the DVI relationships Excel file. The file must contain `sample_id`, `relationship`, and `family_id` columns.
#' @param dvi_sheet Excel sheet for the DVI relationships file.
#'
#' @return Invisibly returns a list with processed objects and output paths.
#'
#' @examples
#' \dontrun{
#' res <- process_forenseq_sample(
#'   report_xlsx = "path/to/sample Flanking Region Report.xlsx",
#'   output_xlsx = "path/to/review_sample.xlsx",
#'   output_txt = "path/to/haplotypes.txt"
#' )
#'
#' res$familias_row
#' head(res$final_calls)
#' res$review
#'
#' res_target <- process_forenseq_sample(
#'   report_xlsx = "path/to/sample Flanking Region Report.xlsx",
#'   output_xlsx = "path/to/review_sample_target.xlsx",
#'   output_txt = "path/to/haplotypes.txt",
#'   write_target_txt = TRUE,
#'   output_target_txt = "path/to/target.txt"
#' )
#' }
#' @export
process_forenseq_sample <- function(
    report_xlsx,
    output_xlsx,
    output_txt = NULL,
    write_target_txt = FALSE,
    output_target_txt = NULL,
    write_parameters_xlsx = TRUE,
    output_parameters_xlsx = NULL,
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
    duplicates_output = NULL,
    dvi = FALSE,
    dvi_relationships = NULL,
    dvi_sheet = 1
) {
  if (!isTRUE(dvi) && !is.null(dvi_relationships)) {
    warning("`dvi_relationships` was provided but `dvi = FALSE`; it will be ignored.")
  }

  dvi_data <- NULL

  if (isTRUE(dvi)) {
    dvi_data <- read_dvi_relationships(
      dvi_relationships = dvi_relationships,
      dvi_sheet = dvi_sheet
    )
  }

  haplotype_data <- read_haplotype_database(
    haplotype_db = haplotype_db,
    db_sheet = db_sheet,
    duplicates_output = duplicates_output
  )

  process_forenseq_sample_core(
    report_xlsx = report_xlsx,
    haplotype_data = haplotype_data,
    output_xlsx = output_xlsx,
    output_txt = output_txt,
    write_target_txt = write_target_txt,
    output_target_txt = output_target_txt,
    write_parameters_xlsx = write_parameters_xlsx,
    output_parameters_xlsx = output_parameters_xlsx,
    haplotype_db = haplotype_db,
    report_sheet = report_sheet,
    db_sheet = db_sheet,
    loci_order = loci_order,
    heterozygote_threshold = heterozygote_threshold,
    extra_signal_threshold = extra_signal_threshold,
    thresholds_by_locus = thresholds_by_locus,
    min_homozygote_reads = min_homozygote_reads,
    min_heterozygous_reads = min_heterozygous_reads,
    review_balance_lower_margin = review_balance_lower_margin,
    review_balance_upper_margin = review_balance_upper_margin,
    review_low_homozygote_multiplier = review_low_homozygote_multiplier,
    dvi = dvi,
    dvi_data = dvi_data
  )
}

resolve_output_txt_path <- function(output_txt, output_dir) {
  if (is.null(output_txt) || output_txt == "") {
    return(file.path(output_dir, "haplotypes.txt"))
  }

  if (dirname(output_txt) == ".") {
    return(file.path(output_dir, output_txt))
  }

  output_txt
}


resolve_output_target_txt_path <- function(
    write_target_txt = FALSE,
    output_target_txt = NULL,
    output_dir = NULL,
    output_txt_path = NULL,
    output_xlsx = NULL
) {
  if (!isTRUE(write_target_txt)) {
    return(NULL)
  }

  if (!is.null(output_target_txt) && output_target_txt != "") {
    if (!is.null(output_dir) && dirname(output_target_txt) == ".") {
      return(file.path(output_dir, output_target_txt))
    }

    return(output_target_txt)
  }

  if (!is.null(output_dir)) {
    return(file.path(output_dir, "target.txt"))
  }

  if (!is.null(output_txt_path) && output_txt_path != "" && dirname(output_txt_path) != ".") {
    return(file.path(dirname(output_txt_path), "target.txt"))
  }

  if (!is.null(output_xlsx) && output_xlsx != "" && dirname(output_xlsx) != ".") {
    return(file.path(dirname(output_xlsx), "target.txt"))
  }

  "target.txt"
}


list_forenseq_xlsx_files <- function(input_dir) {
  files <- list.files(
    path = input_dir,
    pattern = "\\.xlsx$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  files <- files[
    !startsWith(basename(files), "~$")
  ]

  files <- files[
    !grepl("^review_", basename(files), ignore.case = TRUE) &
      !grepl("^haplotypes", basename(files), ignore.case = TRUE) &
      !grepl("^target", basename(files), ignore.case = TRUE) &
      !grepl("^duplicated", basename(files), ignore.case = TRUE) &
      !grepl("^processing_log", basename(files), ignore.case = TRUE) &
      !grepl("^parameters", basename(files), ignore.case = TRUE) &
      !grepl("^log_", basename(files), ignore.case = TRUE)
  ]

  files
}


create_sample_id_with_suffix <- function(sample_id, sample_counter) {
  if (is.null(sample_counter[[sample_id]])) {
    sample_counter[[sample_id]] <- 1L
  } else {
    sample_counter[[sample_id]] <- sample_counter[[sample_id]] + 1L
  }

  occurrence <- sample_counter[[sample_id]]

  if (occurrence == 1L) {
    txt_sample_id <- sample_id
  } else {
    txt_sample_id <- paste0(sample_id, "_", occurrence)
  }

  list(
    sample_counter = sample_counter,
    occurrence = occurrence,
    txt_sample_id = txt_sample_id
  )
}


#' Process a folder of ForenSeq Flanking Region Reports
#'
#' Processes all valid `.xlsx` ForenSeq reports in a folder and creates one
#' combined haplotype TXT file plus one simplified Excel review file per sample.
#'
#' @param input_dir Folder containing ForenSeq report Excel files.
#' @param output_dir Folder where output files will be written.
#' @param output_txt Name or path of the combined haplotype TXT file.
#' @param write_target_txt Logical. If `TRUE`, also writes a TargetSNP_allele TXT file.
#'   Defaults to `FALSE`.
#' @param output_target_txt Optional name or path of the combined TargetSNP_allele
#'   TXT file. Used only when `write_target_txt = TRUE`.
#' @param write_parameters_xlsx Logical. If `TRUE`, writes `parameters.xlsx`
#'   with the processing parameters used in the run. Defaults to `TRUE`.
#' @param output_parameters_xlsx Optional name or path of the processing
#'   parameters Excel file. If `NULL`, `parameters.xlsx` is written in `output_dir`.
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
#' @param continue_on_error Logical. Continue processing if a sample fails.
#' @param write_log Logical. If `TRUE`, writes `processing_log.tsv`.
#'   Defaults to `TRUE`.
#' @param dvi Logical. If `TRUE`, DVI relationship columns are added.
#' @param dvi_relationships Path to the DVI relationships Excel file.
#' @param dvi_sheet Excel sheet for the DVI relationships file.
#'
#' @return Invisibly returns a list with the combined TXT tables, processing log and output paths.
#'
#' @examples
#' \dontrun{
#' res <- process_forenseq_folder(
#'   input_dir = "path/to/forenseq_reports",
#'   output_dir = "path/to/results"
#' )
#'
#' head(res$familias_txt)
#' res$processing_log
#'
#' res_target <- process_forenseq_folder(
#'   input_dir = "path/to/forenseq_reports",
#'   output_dir = "path/to/results_with_target",
#'   write_target_txt = TRUE
#' )
#' }
#' @export
process_forenseq_folder <- function(
    input_dir,
    output_dir,
    output_txt = "haplotypes.txt",
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
    continue_on_error = TRUE,
    write_log = TRUE,
    dvi = FALSE,
    dvi_relationships = NULL,
    dvi_sheet = 1
) {
  if (!dir.exists(input_dir)) {
    stop("Input folder does not exist: ", input_dir)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (!isTRUE(dvi) && !is.null(dvi_relationships)) {
    warning("`dvi_relationships` was provided but `dvi = FALSE`; it will be ignored.")
  }

  if (!isTRUE(write_target_txt) && !is.null(output_target_txt)) {
    warning("`output_target_txt` was provided but `write_target_txt = FALSE`; it will be ignored.")
  }

  dvi_data <- NULL

  if (isTRUE(dvi)) {
    dvi_data <- read_dvi_relationships(
      dvi_relationships = dvi_relationships,
      dvi_sheet = dvi_sheet
    )
  }

  output_txt_path <- resolve_output_txt_path(
    output_txt = output_txt,
    output_dir = output_dir
  )

  output_target_txt_path <- resolve_output_target_txt_path(
    write_target_txt = write_target_txt,
    output_target_txt = output_target_txt,
    output_dir = output_dir
  )

  output_parameters_xlsx_path <- resolve_output_parameters_xlsx_path(
    write_parameters_xlsx = write_parameters_xlsx,
    output_parameters_xlsx = output_parameters_xlsx,
    output_dir = output_dir
  )

  duplicates_output <- file.path(
    output_dir,
    "duplicated_haplotype_database_records.xlsx"
  )

  haplotype_data <- read_haplotype_database(
    haplotype_db = haplotype_db,
    db_sheet = db_sheet,
    duplicates_output = duplicates_output
  )

  files <- list_forenseq_xlsx_files(input_dir)

  if (length(files) == 0) {
    stop("No valid .xlsx files were found in input folder: ", input_dir)
  }

  familias_rows_list <- list()
  familias_target_rows_list <- list()
  log_list <- list()
  sample_counter <- list()

  for (file in files) {
    original_sample_id <- NA_character_
    txt_sample_id <- NA_character_
    output_xlsx_sample <- NA_character_

    iteration <- tryCatch(
      {
        original_sample_id <- read_forenseq_sample_id(
          report_xlsx = file,
          report_sheet = report_sheet
        )

        if (is.na(original_sample_id) || original_sample_id == "") {
          original_sample_id <- tools::file_path_sans_ext(basename(file))
        }

        suffix_info <- create_sample_id_with_suffix(
          sample_id = original_sample_id,
          sample_counter = sample_counter
        )

        sample_counter <- suffix_info$sample_counter
        occurrence <- suffix_info$occurrence
        txt_sample_id <- suffix_info$txt_sample_id

        output_xlsx_sample <- file.path(
          output_dir,
          paste0("review_", sanitize_filename(txt_sample_id), ".xlsx")
        )

        result <- process_forenseq_sample_core(
          report_xlsx = file,
          haplotype_data = haplotype_data,
          output_xlsx = output_xlsx_sample,
          output_txt = NULL,
          write_target_txt = FALSE,
          output_target_txt = NULL,
          write_parameters_xlsx = FALSE,
          output_parameters_xlsx = NULL,
          haplotype_db = haplotype_db,
          report_sheet = report_sheet,
          db_sheet = db_sheet,
          loci_order = loci_order,
          original_sample_id = original_sample_id,
          txt_sample_id = txt_sample_id,
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

        log_ok <- result$process_summary |>
          dplyr::mutate(
            duplicated_sample_id = occurrence > 1L,
            sample_id_occurrence = occurrence
          )

        message("OK: ", basename(file), " -> ", txt_sample_id)

        list(
          status = "OK",
          familias_row = result$familias_row,
          familias_target_row = result$familias_target_row,
          log = log_ok
        )
      },
      error = function(e) {
        error_message <- conditionMessage(e)

        log_error <- tibble::tibble(
          source_file = basename(file),
          original_sample_id = original_sample_id,
          txt_sample_id = txt_sample_id,
          status = "ERROR",
          error_message = error_message,
          output_xlsx = output_xlsx_sample,
          n_loci_observed_in_report = NA_integer_,
          n_loci_processed = NA_integer_,
          n_loci_with_haplotypes_txt = NA_integer_,
          n_loci_with_target_txt = NA_integer_,
          n_homozygous_loci = NA_integer_,
          n_heterozygous_loci = NA_integer_,
          n_no_call_loci = NA_integer_,
          n_review_rows = NA_integer_,
          n_database_not_found_loci = NA_integer_,
          dvi = isTRUE(dvi),
          duplicated_sample_id = NA,
          sample_id_occurrence = NA_integer_
        )

        message("ERROR: ", basename(file), " -> ", error_message)

        if (!continue_on_error) {
          stop(e)
        }

        list(
          status = "ERROR",
          familias_row = NULL,
          familias_target_row = NULL,
          log = log_error
        )
      }
    )

    log_list[[length(log_list) + 1L]] <- iteration$log

    if (identical(iteration$status, "OK")) {
      familias_rows_list[[length(familias_rows_list) + 1L]] <- iteration$familias_row
      familias_target_rows_list[[length(familias_target_rows_list) + 1L]] <- iteration$familias_target_row
    }
  }

  familias_rows <- dplyr::bind_rows(familias_rows_list)
  processing_log <- dplyr::bind_rows(log_list)

  if (nrow(familias_rows) == 0) {
    stop("No samples were processed successfully. The combined TXT was not created.")
  }

  export_familias_txt(
    familias_rows = familias_rows,
    output_txt = output_txt_path
  )

  familias_target_rows <- NULL

  if (isTRUE(write_target_txt)) {
    familias_target_rows <- dplyr::bind_rows(familias_target_rows_list)

    export_familias_txt(
      familias_rows = familias_target_rows,
      output_txt = output_target_txt_path
    )
  }

  if (isTRUE(write_log)) {
    log_path <- file.path(output_dir, "processing_log.tsv")

    utils::write.table(
      processing_log,
      file = log_path,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE,
      na = "",
      fileEncoding = "UTF-8"
    )
  }

  processing_parameters <- NULL

  if (isTRUE(write_parameters_xlsx)) {
    processing_parameters <- create_processing_parameters_table(
      processing_mode = "folder",
      input_dir = input_dir,
      output_dir = output_dir,
      output_txt = output_txt_path,
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
      write_log = write_log,
      dvi = dvi
    )

    export_processing_parameters_xlsx(
      parameters = processing_parameters,
      output_parameters_xlsx = output_parameters_xlsx_path
    )
  }

  message("--------------------------------------------------")
  message("Processing finished.")
  message("Samples OK: ", sum(processing_log$status == "OK", na.rm = TRUE))
  message("Samples with ERROR: ", sum(processing_log$status == "ERROR", na.rm = TRUE))
  message("Combined haplotype TXT created: ", output_txt_path)
  if (isTRUE(write_target_txt)) {
    message("Combined target TXT created: ", output_target_txt_path)
  }
  if (isTRUE(write_parameters_xlsx)) {
    message("Parameters Excel created: ", output_parameters_xlsx_path)
  }
  message("Output folder: ", output_dir)
  message("--------------------------------------------------")

  invisible(
    list(
      familias_txt = familias_rows,
      familias_target_txt = familias_target_rows,
      processing_log = processing_log,
      output_txt = output_txt_path,
      output_target_txt = output_target_txt_path,
      processing_parameters = processing_parameters,
      output_parameters_xlsx = output_parameters_xlsx_path,
      output_dir = output_dir
    )
  )
}

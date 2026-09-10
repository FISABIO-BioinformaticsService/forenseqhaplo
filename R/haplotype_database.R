#' Export bundled haplotype database to Excel
#'
#' Exports the internal haplotype database bundled with the package to an
#' editable Excel file. Users can edit this exported file and pass it to
#' `process_forenseq_folder()` through the `haplotype_db` argument.
#'
#' @param output_path Path to the Excel file to create.
#' @param overwrite Logical. If `FALSE`, the function stops when
#'   `output_path` already exists.
#'
#' @return Invisibly returns `output_path`.
#'
#' @examples
#' \dontrun{
#' # Export the bundled haplotype database to an editable Excel file
#' export_haplotype_database(
#'   output_path = "haplotype_database_export.xlsx",
#'   overwrite = TRUE
#' )
#'
#' # The exported Excel file can be edited and then used as a custom database
#' res <- process_forenseq_folder(
#'   input_dir = "path/to/forenseq_reports",
#'   output_dir = "path/to/results_custom_db",
#'   haplotype_db = "haplotype_database_export.xlsx"
#' )
#' }
#' @export
export_haplotype_database <- function(
    output_path = "haplotype_database.xlsx",
    overwrite = FALSE
) {
  if (file.exists(output_path) && !isTRUE(overwrite)) {
    stop(
      "File already exists: ", output_path, "\n",
      "Use overwrite = TRUE to replace it."
    )
  }

  writexl::write_xlsx(
    list(haplotype_database = haplotype_database),
    path = output_path
  )

  invisible(output_path)
}


#' Load bundled haplotype database
#'
#' Returns the internal haplotype database bundled with the package.
#'
#' @return A tibble/data frame with haplotype database records.
#'
#' @examples
#' db <- load_haplotype_database()
#' head(db)
#' dim(db)
#' @export
load_haplotype_database <- function() {
  haplotype_database
}

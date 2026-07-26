standardize_haplotype_database <- function(db) {
  col_target_snp <- get_column(db, c("TargetSNP"))
  col_target_snp_allele <- get_column(db, c("TargetSNP_allele"))
  col_haplotype <- get_column(db, c("Haplotype"))
  col_sequence <- get_column(db, c("Sequence"))
  col_polymorphic_sites <- get_column(db, c("Polymorphic_sites"))
  col_position_grch38 <- get_column(db, c("Position_GRCh38"))

  out <- tibble::tibble(
    TargetSNP = normalize_locus(db[[col_target_snp]]),
    TargetSNP_allele = clean_text(db[[col_target_snp_allele]]) |>
      stringr::str_to_upper(),
    Haplotype = clean_text(db[[col_haplotype]]),
    Sequence = clean_text(db[[col_sequence]]),
    Polymorphic_sites = clean_text(db[[col_polymorphic_sites]]),
    Position_GRCh38 = clean_text(db[[col_position_grch38]])
  )

  out <- out[
    !is.na(out$TargetSNP) &
      out$TargetSNP != "" &
      !is.na(out$TargetSNP_allele) &
      out$TargetSNP_allele != "" &
      !is.na(out$Haplotype) &
      out$Haplotype != "" &
      !is.na(out$Sequence) &
      out$Sequence != "" &
      !is.na(out$Polymorphic_sites) &
      out$Polymorphic_sites != "" &
      !is.na(out$Position_GRCh38) &
      out$Position_GRCh38 != "",
  ]

  dplyr::distinct(out)
}


add_sequence_key <- function(db) {
  db |>
    dplyr::mutate(
      sequence_key = normalize_sequence(.data$Sequence)
    )
}


find_haplotype_database_duplicates <- function(db) {
  db |>
    add_sequence_key() |>
    dplyr::count(.data$TargetSNP, .data$sequence_key, name = "n") |>
    dplyr::filter(.data$n > 1)
}


#' Validate a haplotype database
#'
#' Validates either the bundled haplotype database or a user-provided Excel
#' haplotype database. The required fields are target SNP, target SNP allele,
#' haplotype, sequence, polymorphic sites and GRCh38 position.
#'
#' The effective matching key is `TargetSNP + Sequence`, with `Sequence`
#' normalized internally during validation and matching.
#'
#' @param haplotype_db Path to a user Excel database. If `NULL`, the bundled
#'   internal database is validated.
#' @param db_sheet Excel sheet to read when `haplotype_db` is provided.
#'
#' @return Invisibly returns `TRUE` if validation succeeds.
#' @export
validate_haplotype_database <- function(
    haplotype_db = NULL,
    db_sheet = 1
) {
  db <- read_haplotype_database(
    haplotype_db = haplotype_db,
    db_sheet = db_sheet
  )

  required_columns <- c(
    "TargetSNP",
    "TargetSNP_allele",
    "Haplotype",
    "Sequence",
    "Polymorphic_sites",
    "Position_GRCh38"
  )

  missing_columns <- setdiff(required_columns, names(db))

  if (length(missing_columns) > 0) {
    stop(
      "The haplotype database is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  empty_haplotypes <- db[
    is.na(db$Haplotype) | db$Haplotype == "",
  ]

  if (nrow(empty_haplotypes) > 0) {
    stop("The haplotype database contains empty Haplotype values.")
  }

  duplicates <- find_haplotype_database_duplicates(db)

  if (nrow(duplicates) > 0) {
    stop("The haplotype database contains duplicated TargetSNP + Sequence keys.")
  }

  invisible(TRUE)
}


read_haplotype_database <- function(
    haplotype_db = NULL,
    db_sheet = 1,
    duplicates_output = NULL
) {
  if (is.null(haplotype_db)) {
    db <- standardize_haplotype_database(load_haplotype_database())
  } else {
    raw_db <- readxl::read_excel(
      path = haplotype_db,
      sheet = db_sheet,
      col_types = "text",
      .name_repair = "unique_quiet"
    )

    db <- standardize_haplotype_database(raw_db)
  }

  if (nrow(db) == 0) {
    stop("The haplotype database is empty after cleaning.")
  }

  duplicates <- find_haplotype_database_duplicates(db)

  if (nrow(duplicates) > 0) {
    duplicated_rows <- db |>
      add_sequence_key() |>
      dplyr::semi_join(
        duplicates,
        by = c("TargetSNP", "sequence_key")
      ) |>
      dplyr::arrange(.data$TargetSNP, .data$sequence_key) |>
      dplyr::select(-.data$sequence_key)

    if (!is.null(duplicates_output)) {
      writexl::write_xlsx(
        list(duplicated_records = duplicated_rows),
        path = duplicates_output
      )
    }

    stop(
      "The haplotype database contains duplicated TargetSNP + Sequence keys.",
      if (!is.null(duplicates_output)) {
        paste0("\nDuplicated records were written to: ", duplicates_output)
      } else {
        ""
      }
    )
  }

  db |>
    dplyr::select(
      dplyr::all_of(c(
        "TargetSNP",
        "TargetSNP_allele",
        "Haplotype",
        "Sequence",
        "Polymorphic_sites",
        "Position_GRCh38"
      ))
    )
}

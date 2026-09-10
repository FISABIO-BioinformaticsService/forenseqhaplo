read_dvi_relationships <- function(
    dvi_relationships,
    dvi_sheet = 1
) {
  if (is.null(dvi_relationships) || is.na(dvi_relationships) || dvi_relationships == "") {
    stop("`dvi_relationships` must be provided when `dvi = TRUE`.")
  }

  relationships_raw <- readxl::read_excel(
    path = dvi_relationships,
    sheet = dvi_sheet,
    col_types = "text",
    .name_repair = "unique_quiet"
  )

  col_sample_id <- get_column(
    relationships_raw,
    c("sample_id", "sample", "id_muestra")
  )

  col_relationship <- get_column(
    relationships_raw,
    c("relationship", "relation", "relacion", "relaci\u00f3n")
  )

  col_family_id <- get_column(
    relationships_raw,
    c("family_id")
  )

  relationships <- tibble::tibble(
    family_id = clean_text(relationships_raw[[col_family_id]]),
    relationship = clean_text(relationships_raw[[col_relationship]]),
    sample_id = clean_text(relationships_raw[[col_sample_id]])
  )

  relationships <- relationships |>
    dplyr::filter(
      !is.na(.data$sample_id),
      .data$sample_id != "",
      !is.na(.data$relationship),
      .data$relationship != "",
      !is.na(.data$family_id),
      .data$family_id != ""
    )

  if (nrow(relationships) == 0) {
    stop("The DVI relationships file is empty after cleaning.")
  }

  duplicated_samples <- relationships |>
    dplyr::count(.data$sample_id, name = "n") |>
    dplyr::filter(.data$n > 1)

  if (nrow(duplicated_samples) > 0) {
    stop(
      "The DVI relationships file contains duplicated sample_id values: ",
      paste(duplicated_samples$sample_id, collapse = ", ")
    )
  }

  relationships
}


get_dvi_relationship_for_sample <- function(
    dvi_data,
    sample_id
) {
  if (is.null(dvi_data)) {
    stop("DVI data are required when `dvi = TRUE`.")
  }

  sample_id_query <- as.character(sample_id)

  hit <- dvi_data |>
    dplyr::filter(.data$sample_id == sample_id_query)

  if (nrow(hit) == 0) {
    stop(
      "No DVI relationship found for sample_id: ",
      sample_id,
      ". The `sample_id` column in the DVI relationships file must match the final TXT sample_id."
    )
  }

  if (nrow(hit) > 1) {
    stop("More than one DVI relationship found for sample_id: ", sample_id)
  }

  hit[1, ]
}

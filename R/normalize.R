normalize_sequence <- function(x) {
  out <- as.character(x)
  out <- stringr::str_replace_all(out, "\u00A0", "")
  out <- stringr::str_replace_all(out, "\\s+", "")
  out <- stringr::str_to_upper(out)
  out[out == ""] <- NA_character_
  out
}


normalize_detected_bases <- function(x) {
  out <- as.character(x)
  out <- stringr::str_replace_all(out, "\u00A0", "")
  out <- stringr::str_replace_all(out, "[\\s,/;]+", "")
  out <- stringr::str_to_upper(out)
  out[out == ""] <- NA_character_
  out
}


normalize_locus <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\u00A0", " ")
  x <- stringr::str_squish(x)

  rs_list <- stringr::str_extract_all(
    stringr::str_to_lower(x),
    "rs\\s*\\d+",
    simplify = FALSE
  )

  out <- vapply(
    rs_list,
    function(z) {
      z <- z[!is.na(z)]

      if (length(z) == 0) {
        NA_character_
      } else {
        z <- stringr::str_replace_all(z, "\\s+", "")
        tail(z, 1)
      }
    },
    character(1)
  )

  out <- ifelse(is.na(out), stringr::str_to_lower(x), out)
  out[out == ""] <- NA_character_
  out
}


normalize_colname <- function(x) {
  out <- as.character(x)
  out <- stringr::str_replace_all(out, "\u00A0", " ")
  out <- stringr::str_squish(out)
  stringr::str_to_lower(out)
}


clean_text <- function(x) {
  out <- as.character(x)
  out <- stringr::str_replace_all(out, "\u00A0", " ")
  out <- stringr::str_squish(out)
  out[out == ""] <- NA_character_
  out
}


parse_read_count <- function(x) {
  out <- as.character(x)
  out <- stringr::str_replace_all(out, "\u00A0", "")
  out <- stringr::str_replace_all(out, "\\s+", "")
  out <- stringr::str_replace_all(out, "[^0-9]", "")
  out[out == ""] <- NA_character_
  as.numeric(out)
}


append_reason <- function(x, reason) {
  x <- ifelse(is.na(x), "", x)
  ifelse(x == "", reason, paste(x, reason, sep = "; "))
}


sanitize_filename <- function(x) {
  x <- clean_text(x)

  if (is.na(x) || x == "") {
    x <- "missing_sample_id"
  }

  x <- stringr::str_replace_all(x, "[/\\\\:*?\"<>|]", "_")
  x <- stringr::str_replace_all(x, "\\s+", "_")
  x
}


build_path_with_suffix <- function(output_xlsx, suffix, ext) {
  file.path(
    dirname(output_xlsx),
    paste0(
      tools::file_path_sans_ext(basename(output_xlsx)),
      suffix,
      ext
    )
  )
}


get_column <- function(df, candidates) {
  names_norm <- normalize_colname(names(df))
  candidates_norm <- normalize_colname(candidates)

  pos <- match(candidates_norm, names_norm)
  pos <- pos[!is.na(pos)]

  if (length(pos) == 0) {
    stop(
      "None of these columns were found: ",
      paste(candidates, collapse = ", ")
    )
  }

  names(df)[pos[1]]
}

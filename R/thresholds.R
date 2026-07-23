normalize_thresholds_by_locus <- function(thresholds_by_locus) {
  if (is.null(thresholds_by_locus) || length(thresholds_by_locus) == 0) {
    return(numeric(0))
  }

  if (is.null(names(thresholds_by_locus)) || any(names(thresholds_by_locus) == "")) {
    stop(
      "`thresholds_by_locus` must be a named numeric vector, ",
      "for example c(rs729172 = 0.1)."
    )
  }

  values <- as.numeric(thresholds_by_locus)
  names_norm <- normalize_locus(names(thresholds_by_locus))

  out <- values
  names(out) <- names_norm

  out[!is.na(names(out)) & names(out) != ""]
}


get_heterozygote_threshold <- function(
    locus,
    heterozygote_threshold = default_heterozygote_threshold(),
    thresholds_by_locus = default_thresholds_by_locus()
) {
  locus <- normalize_locus(locus)
  thresholds_by_locus <- normalize_thresholds_by_locus(thresholds_by_locus)

  if (!is.na(locus) && locus %in% names(thresholds_by_locus)) {
    return(as.numeric(thresholds_by_locus[[locus]]))
  }

  as.numeric(heterozygote_threshold)
}

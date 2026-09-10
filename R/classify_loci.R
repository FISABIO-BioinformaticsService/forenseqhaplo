classify_locus <- function(
    locus_data,
    locus = NULL,
    heterozygote_threshold = default_heterozygote_threshold(),
    extra_signal_threshold = default_extra_signal_threshold(),
    thresholds_by_locus = default_thresholds_by_locus(),
    min_homozygote_reads = default_min_homozygote_reads(),
    min_heterozygous_reads = default_min_heterozygous_reads()
) {
  if (is.null(locus) || is.na(locus) || locus == "") {
    if (nrow(locus_data) > 0 && "iSNP Locus" %in% names(locus_data)) {
      locus <- unique(locus_data$`iSNP Locus`)[1]
    } else {
      locus <- NA_character_
    }
  }

  locus <- normalize_locus(locus)

  locus_data <- locus_data |>
    dplyr::arrange(
      dplyr::desc(.data$Read),
      .data$row_id
    ) |>
    dplyr::mutate(
      sequence_rank = dplyr::row_number()
    )

  original_rows <- nrow(locus_data)
  n_sequences <- dplyr::n_distinct(locus_data$sequence_key)
  n_detected_bases <- dplyr::n_distinct(stats::na.omit(locus_data$`Detected Bases`))

  read_1 <- if (nrow(locus_data) >= 1) locus_data$Read[1] else NA_real_
  read_2 <- if (nrow(locus_data) >= 2) locus_data$Read[2] else NA_real_
  read_3 <- if (nrow(locus_data) >= 3) locus_data$Read[3] else NA_real_

  top_sequence_1 <- if (nrow(locus_data) >= 1) locus_data$Sequence[1] else NA_character_
  top_sequence_2 <- if (nrow(locus_data) >= 2) locus_data$Sequence[2] else NA_character_
  top_sequence_3 <- if (nrow(locus_data) >= 3) locus_data$Sequence[3] else NA_character_

  top_detected_bases_1 <- if (nrow(locus_data) >= 1) locus_data$`Detected Bases`[1] else NA_character_
  top_detected_bases_2 <- if (nrow(locus_data) >= 2) locus_data$`Detected Bases`[2] else NA_character_
  top_detected_bases_3 <- if (nrow(locus_data) >= 3) locus_data$`Detected Bases`[3] else NA_character_

  allele_balance <- ifelse(
    is.na(read_1) || is.na(read_2) || read_1 <= 0,
    NA_real_,
    read_2 / read_1
  )

  top2_total <- ifelse(
    is.na(read_1) || is.na(read_2),
    NA_real_,
    read_1 + read_2
  )

  read_1_top2_fraction <- ifelse(
    is.na(top2_total) || top2_total <= 0,
    NA_real_,
    read_1 / top2_total
  )

  read_2_top2_fraction <- ifelse(
    is.na(top2_total) || top2_total <= 0,
    NA_real_,
    read_2 / top2_total
  )

  extra_signal_ratio <- ifelse(
    is.na(read_1) || is.na(read_3) || read_1 <= 0,
    NA_real_,
    read_3 / read_1
  )

  extra_signal_above_threshold <- !is.na(extra_signal_ratio) &&
    extra_signal_ratio >= extra_signal_threshold

  threshold_used <- get_heterozygote_threshold(
    locus = locus,
    heterozygote_threshold = heterozygote_threshold,
    thresholds_by_locus = thresholds_by_locus
  )

  heterozygote_balance_pass <- !is.na(allele_balance) &&
    allele_balance >= threshold_used

  heterozygote_read_count_pass <- !is.na(read_2) &&
    read_2 >= min_heterozygous_reads

  call_type <- NA_character_
  decision <- NA_character_

  warning_vec <- character(0)
  review_reason_vec <- character(0)
  review_required <- FALSE

  if (nrow(locus_data) == 0) {
    selected_rows <- locus_data[0, ]

    call_type <- "no_call_missing_locus_in_report"
    decision <- "no_call"

    warning_vec <- c(warning_vec, "locus_missing_from_forenseq_report")
    review_reason_vec <- c(review_reason_vec, "locus_missing_from_forenseq_report")
    review_required <- TRUE

  } else if (nrow(locus_data) == 1) {
    if (!is.na(read_1) && read_1 >= min_homozygote_reads) {
      selected_rows <- dplyr::bind_rows(locus_data[1, ], locus_data[1, ]) |>
        dplyr::mutate(allele_copy = dplyr::row_number())

      call_type <- "single_sequence_duplicated"
      decision <- "duplicate_single_observed_sequence"

    } else {
      selected_rows <- locus_data[0, ] |>
        dplyr::mutate(allele_copy = integer())

      call_type <- "not_interpretable_low_homozygote_reads"
      decision <- "not_interpretable"

      warning_vec <- c(warning_vec, "single_sequence_below_min_homozygote_reads")
      review_reason_vec <- c(review_reason_vec, "single_sequence_below_min_homozygote_reads")
      review_required <- TRUE
    }

  } else {
    if (isTRUE(heterozygote_balance_pass) && isTRUE(heterozygote_read_count_pass)) {
      selected_rows <- dplyr::bind_rows(locus_data[1, ], locus_data[2, ]) |>
        dplyr::mutate(allele_copy = dplyr::row_number())

      call_type <- "two_sequences_accepted"
      decision <- "keep_two_most_abundant_sequences"

    } else if (isTRUE(heterozygote_balance_pass) && !isTRUE(heterozygote_read_count_pass)) {
      selected_rows <- locus_data[0, ] |>
        dplyr::mutate(allele_copy = integer())

      call_type <- "not_interpretable_low_heterozygous_reads"
      decision <- "not_interpretable"

      warning_vec <- c(warning_vec, "second_sequence_below_min_heterozygous_reads")
      review_reason_vec <- c(review_reason_vec, "second_sequence_below_min_heterozygous_reads")
      review_required <- TRUE

    } else {
      if (!is.na(read_1) && read_1 >= min_homozygote_reads) {
        selected_rows <- dplyr::bind_rows(locus_data[1, ], locus_data[1, ]) |>
          dplyr::mutate(allele_copy = dplyr::row_number())

        call_type <- "homozygous_collapsed_low_balance"
        decision <- "duplicate_primary_sequence"

        warning_vec <- c(warning_vec, "second_sequence_below_balance_threshold")
        review_reason_vec <- c(review_reason_vec, "second_sequence_collapsed_low_balance")
        review_required <- TRUE

      } else {
        selected_rows <- locus_data[0, ] |>
          dplyr::mutate(allele_copy = integer())

        call_type <- "not_interpretable_low_homozygote_reads"
        decision <- "not_interpretable"

        warning_vec <- c(
          warning_vec,
          "second_sequence_below_balance_threshold",
          "primary_sequence_below_min_homozygote_reads"
        )

        review_reason_vec <- c(
          review_reason_vec,
          "second_sequence_collapsed_low_balance",
          "primary_sequence_below_min_homozygote_reads"
        )

        review_required <- TRUE
      }
    }
  }

  if (nrow(locus_data) >= 3 && isTRUE(extra_signal_above_threshold)) {
    warning_vec <- c(warning_vec, "third_sequence_above_extra_signal_threshold")
    review_reason_vec <- c(review_reason_vec, "third_sequence_above_extra_signal_threshold")
    review_required <- TRUE
  }

  warning <- paste(unique(warning_vec), collapse = "; ")
  warning <- ifelse(warning == "", NA_character_, warning)

  review_reason <- paste(unique(review_reason_vec), collapse = "; ")
  review_reason <- ifelse(review_reason == "", NA_character_, review_reason)

  if (nrow(selected_rows) > 0) {
    selected_rows <- selected_rows |>
      dplyr::mutate(
        `iSNP Locus` = locus,
        allele_copy = dplyr::row_number(),
        call_type = call_type,
        decision = decision,
        allele_balance = allele_balance,
        threshold_used = threshold_used,
        read_1 = read_1,
        read_2 = read_2,
        read_3 = read_3,
        read_1_top2_fraction = read_1_top2_fraction,
        read_2_top2_fraction = read_2_top2_fraction,
        extra_signal_ratio = extra_signal_ratio,
        extra_signal_above_threshold = extra_signal_above_threshold,
        min_homozygote_reads = min_homozygote_reads,
        min_heterozygous_reads = min_heterozygous_reads,
        original_rows = original_rows,
        n_sequences = n_sequences,
        n_detected_bases = n_detected_bases,
        warning = warning,
        review_required = review_required,
        review_reason = review_reason
      )
  }

  selected_copy_counts <- selected_rows |>
    dplyr::count(.data$row_id, name = "selected_copies")

  all_sequences <- locus_data |>
    dplyr::left_join(
      selected_copy_counts,
      by = "row_id"
    ) |>
    dplyr::mutate(
      selected_copies = dplyr::coalesce(.data$selected_copies, 0L),
      selected_for_call = .data$selected_copies > 0,
      discard_reason = dplyr::case_when(
        .data$selected_for_call ~ NA_character_,
        call_type == "not_interpretable_low_heterozygous_reads" ~ "not_selected_low_heterozygous_reads",
        call_type == "not_interpretable_low_homozygote_reads" ~ "not_selected_low_homozygote_reads",
        call_type == "not_interpretable_no_sequence" ~ "no_sequence_observed",
        TRUE ~ "not_selected_for_final_call"
      ),
      call_type = call_type,
      decision = decision,
      allele_balance = allele_balance,
      threshold_used = threshold_used,
      read_1 = read_1,
      read_2 = read_2,
      read_3 = read_3,
      read_1_top2_fraction = read_1_top2_fraction,
      read_2_top2_fraction = read_2_top2_fraction,
      extra_signal_ratio = extra_signal_ratio,
      extra_signal_above_threshold = extra_signal_above_threshold,
      min_homozygote_reads = min_homozygote_reads,
      min_heterozygous_reads = min_heterozygous_reads,
      original_rows = original_rows,
      n_sequences = n_sequences,
      n_detected_bases = n_detected_bases,
      warning = warning,
      review_required = review_required,
      review_reason = review_reason
    )

  locus_summary <- tibble::tibble(
    `iSNP Locus` = locus,
    call_type = call_type,
    decision = decision,
    allele_balance = allele_balance,
    threshold_used = threshold_used,
    read_1 = read_1,
    read_2 = read_2,
    read_3 = read_3,
    read_1_top2_fraction = read_1_top2_fraction,
    read_2_top2_fraction = read_2_top2_fraction,
    extra_signal_ratio = extra_signal_ratio,
    extra_signal_above_threshold = extra_signal_above_threshold,
    min_homozygote_reads = min_homozygote_reads,
    min_heterozygous_reads = min_heterozygous_reads,
    original_rows = original_rows,
    n_sequences = n_sequences,
    n_detected_bases = n_detected_bases,
    top_sequence_1 = top_sequence_1,
    top_sequence_2 = top_sequence_2,
    top_sequence_3 = top_sequence_3,
    top_detected_bases_1 = top_detected_bases_1,
    top_detected_bases_2 = top_detected_bases_2,
    top_detected_bases_3 = top_detected_bases_3,
    review_required = review_required,
    review_reason = review_reason,
    warning = warning
  )

  list(
    selected_rows = selected_rows,
    locus_summary = locus_summary,
    all_sequences = all_sequences
  )
}


classify_all_loci <- function(
    report_data,
    loci_order = NULL,
    heterozygote_threshold = default_heterozygote_threshold(),
    extra_signal_threshold = default_extra_signal_threshold(),
    thresholds_by_locus = default_thresholds_by_locus(),
    min_homozygote_reads = default_min_homozygote_reads(),
    min_heterozygous_reads = default_min_heterozygous_reads()
) {
  duplicated_sequences <- report_data |>
    dplyr::count(.data$`iSNP Locus`, .data$sequence_key, name = "n") |>
    dplyr::filter(.data$n > 1)

  if (nrow(duplicated_sequences) > 0) {
    warning(
      "Duplicated Sequence values were found within the same locus. ",
      "The sequence-based algorithm assumes one row per Sequence."
    )
  }

  observed_loci <- unique(report_data$`iSNP Locus`)
  observed_loci <- observed_loci[!is.na(observed_loci) & observed_loci != ""]

  if (is.null(loci_order)) {
    loci <- observed_loci
  } else {
    expected_loci <- normalize_locus(loci_order)
    expected_loci <- expected_loci[!is.na(expected_loci) & expected_loci != ""]
    expected_loci <- unique(expected_loci)

    loci <- unique(c(
      expected_loci,
      setdiff(observed_loci, expected_loci)
    ))
  }

  results <- lapply(
    loci,
    function(locus) {
      locus_data <- report_data |>
        dplyr::filter(.data$`iSNP Locus` == locus)

      classify_locus(
        locus_data = locus_data,
        locus = locus,
        heterozygote_threshold = heterozygote_threshold,
        extra_signal_threshold = extra_signal_threshold,
        thresholds_by_locus = thresholds_by_locus,
        min_homozygote_reads = min_homozygote_reads,
        min_heterozygous_reads = min_heterozygous_reads
      )
    }
  )

  selected_rows <- dplyr::bind_rows(
    lapply(results, function(x) x$selected_rows)
  )

  locus_summary <- dplyr::bind_rows(
    lapply(results, function(x) x$locus_summary)
  )

  all_sequences <- dplyr::bind_rows(
    lapply(results, function(x) x$all_sequences)
  )

  if (nrow(selected_rows) > 0) {
    check <- selected_rows |>
      dplyr::count(.data$`iSNP Locus`, name = "n") |>
      dplyr::filter(.data$n != 2)

    if (nrow(check) > 0) {
      stop(
        "Some loci do not have exactly two selected rows after sequence-based classification."
      )
    }
  }

  list(
    selected_rows = selected_rows,
    locus_summary = locus_summary,
    all_sequences = all_sequences
  )
}

annotate_with_haplotype_db <- function(df, haplotype_db) {
  haplotype_db_join <- haplotype_db |>
    dplyr::mutate(
      sequence_key = normalize_sequence(.data$Sequence)
    )

  df |>
    dplyr::left_join(
      haplotype_db_join |>
        dplyr::select(
          dplyr::all_of(c(
            "TargetSNP",
            "sequence_key",
            DB_Sequence = "Sequence",
            "Haplotype",
            "TargetSNP_allele",
            "Polymorphic_sites",
            "Position_GRCh38"
          ))
        ),
      by = c(
        "iSNP Locus" = "TargetSNP",
        "sequence_key" = "sequence_key"
      )
    ) |>
    dplyr::mutate(
      TargetSNP = .data[["iSNP Locus"]],
      match_status = dplyr::if_else(
        !is.na(.data$Haplotype),
        "matched",
        "not_found"
      ),
      DB_Sequence = dplyr::coalesce(.data$DB_Sequence, "NA"),
      Haplotype = dplyr::coalesce(.data$Haplotype, "NA"),
      TargetSNP_allele = dplyr::coalesce(.data$TargetSNP_allele, "NA"),
      Polymorphic_sites = dplyr::coalesce(.data$Polymorphic_sites, "NA"),
      Position_GRCh38 = dplyr::coalesce(.data$Position_GRCh38, "NA")
    )
}


is_valid_target_allele <- function(x) {
  !is.na(x) &
    x != "" &
    x != "NA"
}


first_valid <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "NA"]

  if (length(x) == 0) {
    NA_character_
  } else {
    x[[1]]
  }
}


create_familias_outputs <- function(
    selected_evidence,
    locus_summary,
    all_sequences_annotated,
    review_balance_lower_margin = default_review_balance_lower_margin(),
    review_balance_upper_margin = default_review_balance_upper_margin(),
    review_low_homozygote_multiplier = default_review_low_homozygote_multiplier()
) {
  final_calls_list <- list()

  for (i in seq_len(nrow(locus_summary))) {
    summary_row <- locus_summary[i, ]
    locus <- summary_row$`iSNP Locus`

    selected <- selected_evidence |>
      dplyr::filter(.data$`iSNP Locus` == locus) |>
      dplyr::arrange(.data$allele_copy)

    observed <- all_sequences_annotated |>
      dplyr::filter(.data$`iSNP Locus` == locus) |>
      dplyr::arrange(.data$sequence_rank)

    has_two_selected <- nrow(selected) == 2
    all_matched <- has_two_selected && all(selected$match_status == "matched")

    final_haplotype_1 <- if (has_two_selected) selected$Haplotype[1] else NA_character_
    final_haplotype_2 <- if (has_two_selected) selected$Haplotype[2] else NA_character_

    final_target_allele_1 <- if (has_two_selected) selected$TargetSNP_allele[1] else NA_character_
    final_target_allele_2 <- if (has_two_selected) selected$TargetSNP_allele[2] else NA_character_

    final_sequence_1 <- if (has_two_selected) selected$Sequence[1] else NA_character_
    final_sequence_2 <- if (has_two_selected) selected$Sequence[2] else NA_character_

    final_db_sequence_1 <- if (has_two_selected) selected$DB_Sequence[1] else NA_character_
    final_db_sequence_2 <- if (has_two_selected) selected$DB_Sequence[2] else NA_character_

    report_haplotype_1 <- if (nrow(observed) >= 1) observed$Haplotype[1] else NA_character_
    report_haplotype_2 <- if (nrow(observed) >= 2) observed$Haplotype[2] else NA_character_

    report_target_allele_1 <- if (nrow(observed) >= 1) observed$TargetSNP_allele[1] else NA_character_
    report_target_allele_2 <- if (nrow(observed) >= 2) observed$TargetSNP_allele[2] else NA_character_

    report_sequence_1 <- if (nrow(observed) >= 1) observed$Sequence[1] else NA_character_
    report_sequence_2 <- if (nrow(observed) >= 2) observed$Sequence[2] else NA_character_

    report_db_sequence_1 <- if (nrow(observed) >= 1) observed$DB_Sequence[1] else NA_character_
    report_db_sequence_2 <- if (nrow(observed) >= 2) observed$DB_Sequence[2] else NA_character_

    report_match_status_1 <- if (nrow(observed) >= 1) observed$match_status[1] else NA_character_
    report_match_status_2 <- if (nrow(observed) >= 2) observed$match_status[2] else NA_character_

    polymorphic_sites <- first_valid(c(selected$Polymorphic_sites, observed$Polymorphic_sites))
    position_grch38 <- first_valid(c(selected$Position_GRCh38, observed$Position_GRCh38))

    included_in_txt <- isTRUE(all_matched)

    txt_value <- if (included_in_txt) {
      paste(final_haplotype_1, final_haplotype_2, sep = ",")
    } else {
      NA_character_
    }

    included_in_target_txt <- included_in_txt &&
      isTRUE(is_valid_target_allele(final_target_allele_1)) &&
      isTRUE(is_valid_target_allele(final_target_allele_2))

    txt_target_value <- if (included_in_target_txt) {
      paste(final_target_allele_1, final_target_allele_2, sep = ",")
    } else {
      NA_character_
    }

    zygosity <- if (included_in_txt) {
      if (identical(as.character(final_haplotype_1), as.character(final_haplotype_2))) {
        "homozygous"
      } else {
        "heterozygous"
      }
    } else if (identical(as.character(summary_row$call_type), "no_call_missing_locus_in_report")) {
      "no_call"
    } else if (has_two_selected) {
      "not_found"
    } else {
      "not_interpretable"
    }

    txt_action <- NA_character_
    txt_reason <- NA_character_

    if (included_in_txt) {
      txt_action <- "locus_included_in_txt"
      txt_reason <- NA_character_
    } else if (identical(as.character(summary_row$call_type), "no_call_missing_locus_in_report")) {
      txt_action <- "locus_excluded_from_txt_no_call"
      txt_reason <- "no_call_missing_locus_in_report"
    } else if (!has_two_selected) {
      txt_action <- "locus_excluded_from_txt_not_interpretable"
      txt_reason <- summary_row$call_type
    } else {
      txt_action <- "locus_excluded_from_txt_not_found_in_database"
      txt_reason <- "one_or_more_selected_sequences_not_found_in_database"
    }

    algorithm_note <- summary_row$review_reason
    algorithm_warning <- summary_row$warning

    allele_balance <- summary_row$allele_balance
    threshold_used <- summary_row$threshold_used
    read_1 <- summary_row$read_1
    min_homozygote_reads <- summary_row$min_homozygote_reads

    balance_lower <- threshold_used * (1 - review_balance_lower_margin)
    balance_upper <- threshold_used * (1 + review_balance_upper_margin)

    balance_near_threshold <- !is.na(allele_balance) &&
      !is.na(threshold_used) &&
      allele_balance >= balance_lower &&
      allele_balance <= balance_upper

    homozygote_near_min_reads <- zygosity == "homozygous" &&
      summary_row$call_type %in% c(
        "single_sequence_duplicated",
        "homozygous_collapsed_low_balance"
      ) &&
      !is.na(read_1) &&
      !is.na(min_homozygote_reads) &&
      read_1 < (min_homozygote_reads * review_low_homozygote_multiplier)

    review_reason <- NA_character_
    review_level <- "none"

    if (identical(as.character(summary_row$call_type), "no_call_missing_locus_in_report")) {
      review_reason <- append_reason(
        review_reason,
        "locus_missing_from_forenseq_report"
      )
      review_level <- "critical"
    } else if (!included_in_txt) {
      review_reason <- append_reason(
        review_reason,
        "locus_not_reported_in_txt"
      )
      review_level <- "critical"
    }

    if (identical(txt_action, "locus_excluded_from_txt_not_found_in_database")) {
      review_reason <- append_reason(
        review_reason,
        "selected_sequence_not_found_in_database"
      )
      review_level <- "critical"
    }

    if (isTRUE(summary_row$extra_signal_above_threshold)) {
      review_reason <- append_reason(
        review_reason,
        "third_sequence_above_extra_signal_threshold"
      )
      if (review_level == "none") {
        review_level <- "warning"
      }
    }

    if (isTRUE(balance_near_threshold)) {
      review_reason <- append_reason(
        review_reason,
        "allele_balance_near_threshold"
      )
      if (review_level == "none") {
        review_level <- "warning"
      }
    }

    if (isTRUE(homozygote_near_min_reads)) {
      review_reason <- append_reason(
        review_reason,
        "homozygous_call_near_minimum_reads"
      )
      if (review_level == "none") {
        review_level <- "warning"
      }
    }

    if (identical(as.character(summary_row$call_type), "not_interpretable_low_heterozygous_reads")) {
      review_reason <- append_reason(
        review_reason,
        "second_sequence_below_min_heterozygous_reads"
      )
      review_level <- "critical"
    }

    review_required <- review_level %in% c("warning", "critical")

    final_calls_list[[length(final_calls_list) + 1L]] <- tibble::tibble(
      TargetSNP = locus,
      `iSNP Locus` = locus,
      call_type = summary_row$call_type,
      zygosity = zygosity,
      decision = summary_row$decision,

      Haplotype_1 = final_haplotype_1,
      Haplotype_2 = final_haplotype_2,
      TargetSNP_allele_1 = final_target_allele_1,
      TargetSNP_allele_2 = final_target_allele_2,
      Sequence_1 = final_sequence_1,
      Sequence_2 = final_sequence_2,
      DB_Sequence_1 = final_db_sequence_1,
      DB_Sequence_2 = final_db_sequence_2,

      report_Haplotype_1 = report_haplotype_1,
      report_Haplotype_2 = report_haplotype_2,
      report_TargetSNP_allele_1 = report_target_allele_1,
      report_TargetSNP_allele_2 = report_target_allele_2,
      report_Sequence_1 = report_sequence_1,
      report_Sequence_2 = report_sequence_2,
      report_DB_Sequence_1 = report_db_sequence_1,
      report_DB_Sequence_2 = report_db_sequence_2,
      report_match_status_1 = report_match_status_1,
      report_match_status_2 = report_match_status_2,

      Read_1 = summary_row$read_1,
      Read_2 = summary_row$read_2,
      allele_balance = summary_row$allele_balance,
      threshold_used = summary_row$threshold_used,
      balance_lower_review_limit = balance_lower,
      balance_upper_review_limit = balance_upper,
      balance_near_threshold = balance_near_threshold,
      read_1_top2_fraction = summary_row$read_1_top2_fraction,
      read_2_top2_fraction = summary_row$read_2_top2_fraction,
      extra_signal_ratio = summary_row$extra_signal_ratio,
      extra_signal_above_threshold = summary_row$extra_signal_above_threshold,
      min_homozygote_reads = summary_row$min_homozygote_reads,
      min_heterozygous_reads = summary_row$min_heterozygous_reads,
      homozygote_near_min_reads = homozygote_near_min_reads,
      n_sequences = summary_row$n_sequences,
      n_detected_bases = summary_row$n_detected_bases,
      Polymorphic_sites = polymorphic_sites,
      Position_GRCh38 = position_grch38,
      included_in_txt = included_in_txt,
      txt_value = txt_value,
      included_in_target_txt = included_in_target_txt,
      txt_target_value = txt_target_value,
      txt_action = txt_action,
      txt_reason = txt_reason,
      review_required = review_required,
      review_level = review_level,
      review_reason = review_reason,
      algorithm_note = algorithm_note,
      warning = algorithm_warning
    )
  }

  final_calls <- dplyr::bind_rows(final_calls_list)

  included_loci <- final_calls |>
    dplyr::filter(.data$included_in_txt) |>
    dplyr::pull("iSNP Locus")

  included_target_loci <- final_calls |>
    dplyr::filter(.data$included_in_target_txt) |>
    dplyr::pull("iSNP Locus")

  familias_input <- selected_evidence |>
    dplyr::filter(.data$`iSNP Locus` %in% included_loci) |>
    dplyr::arrange(.data$`iSNP Locus`, .data$allele_copy)

  familias_target_input <- selected_evidence |>
    dplyr::filter(.data$`iSNP Locus` %in% included_target_loci) |>
    dplyr::arrange(.data$`iSNP Locus`, .data$allele_copy) |>
    dplyr::mutate(
      Haplotype = .data$TargetSNP_allele
    )

  review <- final_calls |>
    dplyr::filter(.data$review_required == TRUE) |>
    dplyr::mutate(
      final_Haplotype_1 = .data$Haplotype_1,
      final_Haplotype_2 = .data$Haplotype_2,
      final_TargetSNP_allele_1 = .data$TargetSNP_allele_1,
      final_TargetSNP_allele_2 = .data$TargetSNP_allele_2,
      final_Sequence_1 = .data$Sequence_1,
      final_Sequence_2 = .data$Sequence_2,
      final_DB_Sequence_1 = .data$DB_Sequence_1,
      final_DB_Sequence_2 = .data$DB_Sequence_2,
      ForenSeq_sequence = paste(
        stats::na.omit(c(.data$Sequence_1, .data$Sequence_2)),
        collapse = " | "
      ),
      DB_Sequence = paste(
        stats::na.omit(c(.data$DB_Sequence_1, .data$DB_Sequence_2)),
        collapse = " | "
      ),

      Haplotype_1 = .data$report_Haplotype_1,
      Haplotype_2 = .data$report_Haplotype_2,
      TargetSNP_allele_1 = .data$report_TargetSNP_allele_1,
      TargetSNP_allele_2 = .data$report_TargetSNP_allele_2,
      Sequence_1 = .data$report_Sequence_1,
      Sequence_2 = .data$report_Sequence_2
    )

  all_sequences <- all_sequences_annotated |>
    dplyr::mutate(
      TargetSNP = .data[["iSNP Locus"]],
      ForenSeq_sequence = .data$Sequence
    ) |>
    dplyr::left_join(
      final_calls |>
        dplyr::transmute(
          `iSNP Locus` = .data[["iSNP Locus"]],
          zygosity_locus = .data$zygosity,
          included_in_txt_locus = .data$included_in_txt,
          txt_value_locus = .data$txt_value,
          included_in_target_txt_locus = .data$included_in_target_txt,
          txt_target_value_locus = .data$txt_target_value
        ),
      by = "iSNP Locus"
    )

  list(
    familias_input = familias_input,
    familias_target_input = familias_target_input,
    final_calls = final_calls,
    review = review,
    all_sequences = all_sequences
  )
}

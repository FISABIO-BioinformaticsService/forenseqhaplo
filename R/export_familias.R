create_familias_row <- function(
    familias_input,
    sample_id,
    loci_order = default_loci_order(),
    include_missing_loci = TRUE,
    missing_value = "",
    dvi_info = NULL
) {
  if (nrow(familias_input) > 0) {
    check <- familias_input |>
      dplyr::count(.data$`iSNP Locus`, name = "n") |>
      dplyr::filter(.data$n != 2)

    if (nrow(check) > 0) {
      stop("Cannot create Familias row: some loci do not have exactly two rows.")
    }

    check_na <- familias_input |>
      dplyr::filter(
        is.na(.data$Haplotype) |
          .data$Haplotype == "NA" |
          .data$Haplotype == ""
      )

    if (nrow(check_na) > 0) {
      stop("Cannot create Familias row: empty or NA haplotypes found.")
    }

    haplotypes_by_locus <- familias_input |>
      dplyr::arrange(.data$`iSNP Locus`, .data$allele_copy) |>
      dplyr::group_by(.data$`iSNP Locus`) |>
      dplyr::summarise(
        familias_value = paste(.data$Haplotype, collapse = ","),
        .groups = "drop"
      )
  } else {
    haplotypes_by_locus <- tibble::tibble(
      `iSNP Locus` = character(),
      familias_value = character()
    )
  }

  loci_order <- as.character(loci_order)
  loci_order <- normalize_locus(loci_order)
  loci_order <- loci_order[!is.na(loci_order) & loci_order != ""]
  loci_order <- unique(loci_order)

  extra_loci <- setdiff(haplotypes_by_locus$`iSNP Locus`, loci_order)

  if (length(extra_loci) > 0) {
    warning(
      "Some loci have haplotypes but are not present in `loci_order`; ",
      "they will not be included in the TXT: ",
      paste(extra_loci, collapse = ", ")
    )
  }

  if (include_missing_loci) {
    final_loci <- loci_order
  } else {
    final_loci <- loci_order[loci_order %in% haplotypes_by_locus$`iSNP Locus`]
  }

  values <- rep(missing_value, length(final_loci))
  names(values) <- final_loci

  idx <- match(haplotypes_by_locus$`iSNP Locus`, names(values))
  idx_valid <- !is.na(idx)

  values[idx[idx_valid]] <- haplotypes_by_locus$familias_value[idx_valid]

  if (is.null(dvi_info)) {
    values_with_id <- c(
      sample_id = sample_id,
      values
    )
  } else {
    values_with_id <- c(
      sample_id = as.character(dvi_info$sample_id),
      relationship = as.character(dvi_info$relationship),
      family_id = as.character(dvi_info$family_id),
      values
    )
  }

  as.data.frame(
    as.list(values_with_id),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


export_familias_txt <- function(familias_rows, output_txt) {
  utils::write.table(
    familias_rows,
    file = output_txt,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "",
    fileEncoding = "UTF-8"
  )

  invisible(familias_rows)
}


add_sample_metadata <- function(
    df,
    original_sample_id,
    txt_sample_id,
    source_file
) {
  df |>
    dplyr::mutate(
      sample_id = txt_sample_id,
      original_sample_id = original_sample_id,
      source_file = basename(source_file),
      .before = 1
    )
}


excel_explanation_dictionary <- function(kind) {
  dictionaries <- list(
    call_type = c(
      no_call_missing_locus_in_report = "No call was assigned because no sequence reads were detected for this locus in the ForenSeq report.",
      not_interpretable_no_sequence = "No genotype could be interpreted because no sequence reads were detected at this locus.",
      single_sequence_duplicated = "A single sequence met the minimum read depth threshold and was called as homozygous.",
      not_interpretable_low_homozygote_reads = "No genotype could be interpreted because the read depth fell below the minimum threshold required for a homozygous call.",
      not_interpretable_low_heterozygous_reads = "No genotype could be interpreted because the second sequence met the allele balance threshold but fell below the minimum read depth required for a heterozygous call.",
      two_sequences_accepted = "The two sequences with the highest read counts were accepted because their allele balance met the threshold for this locus.",
      homozygous_collapsed_low_balance = "The second sequence fell below the allele balance threshold; therefore, the locus was called as homozygous for the primary sequence."
    ),
    decision = c(
      no_call = "No call was assigned because no sequence reads were detected for this locus in the ForenSeq report.",
      not_interpretable = "No definitive alleles could be determined for this locus.",
      duplicate_single_observed_sequence = "The single detected sequence was called as a homozygous genotype.",
      keep_two_most_abundant_sequences = "The two sequences with the highest read counts were retained as the final genotype.",
      duplicate_primary_sequence = "The locus was called as homozygous because the minor sequence failed to meet the allele balance threshold."
    ),
    warning = c(
      locus_missing_from_forenseq_report = "No sequence reads were detected for this locus in the ForenSeq report.",
      no_sequence_observed = "No sequence reads were detected at this locus.",
      single_sequence_below_min_homozygote_reads = "A single sequence was detected, but its read count fell below the minimum threshold required for a homozygous call.",
      second_sequence_below_balance_threshold = "The minor sequence failed to meet the required allele balance threshold.",
      second_sequence_below_min_heterozygous_reads = "The minor sequence met the allele balance threshold but fell below the minimum read depth required for a heterozygous call.",
      primary_sequence_below_min_homozygote_reads = "The primary sequence lacked the minimum read depth required to confirm a homozygous call.",
      third_sequence_above_extra_signal_threshold = "A third sequence exceeded the analytical threshold for extra peaks and requires manual review."
    ),
    algorithm_note = c(
      locus_missing_from_forenseq_report = "No sequence reads were detected for this locus in the ForenSeq report.",
      no_sequence_observed = "No sequence reads were detected at this locus.",
      single_sequence_below_min_homozygote_reads = "The single detected sequence was filtered out because its read depth fell below the minimum threshold.",
      second_sequence_collapsed_low_balance = "The minor sequence was filtered out because its allele balance fell below the threshold.",
      second_sequence_below_min_heterozygous_reads = "The minor sequence was filtered out because its read depth fell below the minimum threshold required for a heterozygous call.",
      primary_sequence_below_min_homozygote_reads = "The primary sequence lacked the minimum read depth required to confirm a homozygous call.",
      third_sequence_above_extra_signal_threshold = "A third sequence exceeded the extra-signal analytical threshold."
    ),
    txt_action = c(
      locus_included_in_txt = "This locus was included in the haplotype TXT file because both alleles matched the haplotype reference database.",
      locus_excluded_from_txt_no_call = "This locus was excluded from the haplotype TXT file because no sequence reads were detected in the ForenSeq report and it was flagged as 'no_call'.",
      locus_excluded_from_txt_not_interpretable = "This locus was excluded from the haplotype TXT file because the calling algorithm failed to generate an interpretable genotype.",
      locus_excluded_from_txt_not_found_in_database = "This locus was excluded from the haplotype TXT file because one or both sequences were missing from the reference database."
    ),
    txt_reason = c(
      no_call_missing_locus_in_report = "No sequence reads were detected for this locus in the ForenSeq report and it was flagged as 'no_call' during manual review.",
      one_or_more_selected_sequences_not_found_in_database = "One or both called sequences were missing from the haplotype reference database.",
      not_interpretable_no_sequence = "No genotype could be interpreted because no sequence reads were detected at this locus.",
      not_interpretable_low_homozygote_reads = "No genotype could be interpreted because the read depth fell below the minimum threshold required for a homozygous call.",
      not_interpretable_low_heterozygous_reads = "No genotype could be interpreted because the second sequence met the allele balance threshold but fell below the minimum read depth required for a heterozygous call.",
      single_sequence_duplicated = "A single sequence was called as a homozygous genotype.",
      two_sequences_accepted = "The two sequences with the highest read counts were accepted.",
      homozygous_collapsed_low_balance = "The minor sequence fell below the allele balance threshold; therefore, the locus was called as homozygous for the primary sequence."
    ),
    review_reason = c(
      locus_missing_from_forenseq_report = "This locus requires review because no sequence reads were detected in the ForenSeq report.",
      locus_not_reported_in_txt = "This locus requires manual review because it was missing from the exported haplotype TXT file.",
      selected_sequence_not_found_in_database = "This locus requires manual review because one or both called sequences were missing from the database.",
      second_sequence_below_min_heterozygous_reads = "This locus requires manual review because the second sequence met the allele balance threshold but fell below the minimum read depth required for a heterozygous call.",
      third_sequence_above_extra_signal_threshold = "This locus requires manual review because a third sequence exceeded the extra-signal analytical threshold.",
      allele_balance_near_threshold = "This locus requires manual review because the allele balance is near the threshold used to differentiate heterozygous and homozygous calls.",
      homozygous_call_near_minimum_reads = "This locus requires manual review because a homozygous call was accepted with read depth near the minimum analytical threshold."
    ),
    discard_reason = c(
      not_selected_low_homozygote_reads = "This sequence was rejected because the locus lacked sufficient read depth for a valid homozygous call.",
      not_selected_low_heterozygous_reads = "This sequence was rejected because the second sequence fell below the minimum read depth required for a valid heterozygous call.",
      no_sequence_observed = "No sequence reads were detected at this locus.",
      not_selected_for_final_call = "This sequence was not selected for the final genotype call."
    ),
    match_status = c(
      matched = "This sequence was matched in the haplotype reference database.",
      not_found = "This sequence was not found in the haplotype reference database."
    )
  )

  dictionaries[[kind]]
}


explain_code_string <- function(x, dictionary) {
  vapply(
    x,
    function(value) {
      if (is.na(value) || value == "") {
        return(NA_character_)
      }

      codes <- unlist(strsplit(as.character(value), ";", fixed = TRUE))
      codes <- trimws(codes)
      codes <- codes[codes != ""]

      if (length(codes) == 0) {
        return(NA_character_)
      }

      explanations <- dictionary[codes]
      missing <- is.na(explanations)

      if (any(missing)) {
        explanations[missing] <- paste0("Unmapped code: ", codes[missing])
      }

      paste(unname(explanations), collapse = " | ")
    },
    character(1)
  )
}


combine_explanations <- function(...) {
  parts <- list(...)

  if (length(parts) == 0) {
    return(character(0))
  }

  n <- max(vapply(parts, length, integer(1)))

  parts <- lapply(
    parts,
    function(x) {
      if (length(x) == n) {
        x
      } else if (length(x) == 1) {
        rep(x, n)
      } else {
        stop("Cannot combine explanation vectors with incompatible lengths.")
      }
    }
  )

  vapply(
    seq_len(n),
    function(i) {
      values <- vapply(parts, function(x) x[[i]], character(1))
      values <- values[!is.na(values) & values != ""]
      values <- unique(values)

      if (length(values) == 0) {
        NA_character_
      } else {
        paste(values, collapse = " | ")
      }
    },
    character(1)
  )
}


add_excel_explanations <- function(df) {
  out <- df

  if ("call_type" %in% names(out)) {
    out$call_explanation <- explain_code_string(
      out$call_type,
      excel_explanation_dictionary("call_type")
    )
  }

  if ("decision" %in% names(out)) {
    out$decision_explanation <- explain_code_string(
      out$decision,
      excel_explanation_dictionary("decision")
    )
  }

  if ("txt_action" %in% names(out)) {
    txt_action_explanation <- explain_code_string(
      out$txt_action,
      excel_explanation_dictionary("txt_action")
    )
  } else {
    txt_action_explanation <- rep(NA_character_, nrow(out))
  }

  if ("txt_reason" %in% names(out)) {
    txt_reason_explanation <- explain_code_string(
      out$txt_reason,
      excel_explanation_dictionary("txt_reason")
    )
  } else {
    txt_reason_explanation <- rep(NA_character_, nrow(out))
  }

  if ("txt_action" %in% names(out) || "txt_reason" %in% names(out)) {
    out$txt_explanation <- combine_explanations(
      txt_action_explanation,
      txt_reason_explanation
    )
  }

  if ("included_in_txt_locus" %in% names(out)) {
    out$locus_txt_explanation <- dplyr::case_when(
      is.na(out$included_in_txt_locus) ~ NA_character_,
      out$included_in_txt_locus ~ "This locus was included in the haplotype TXT.",
      !out$included_in_txt_locus ~ "This locus was not included in the haplotype TXT.",
      TRUE ~ NA_character_
    )
  }

  if ("included_in_target_txt_locus" %in% names(out)) {
    out$locus_target_txt_explanation <- dplyr::case_when(
      is.na(out$included_in_target_txt_locus) ~ NA_character_,
      out$included_in_target_txt_locus ~ "This locus was included in the target TXT.",
      !out$included_in_target_txt_locus ~ "This locus was not included in the target TXT.",
      TRUE ~ NA_character_
    )
  }

  if ("review_reason" %in% names(out)) {
    out$review_explanation <- explain_code_string(
      out$review_reason,
      excel_explanation_dictionary("review_reason")
    )
  }

  if ("algorithm_note" %in% names(out)) {
    out$algorithm_explanation <- explain_code_string(
      out$algorithm_note,
      excel_explanation_dictionary("algorithm_note")
    )
  }

  if ("warning" %in% names(out)) {
    out$warning_explanation <- explain_code_string(
      out$warning,
      excel_explanation_dictionary("warning")
    )
  }

  if ("discard_reason" %in% names(out)) {
    discard_explanation <- explain_code_string(
      out$discard_reason,
      excel_explanation_dictionary("discard_reason")
    )

    if ("selected_for_call" %in% names(out)) {
      selected_explanation <- dplyr::case_when(
        is.na(out$selected_for_call) ~ NA_character_,
        out$selected_for_call ~ "This observed sequence was selected for the final result.",
        !out$selected_for_call ~ NA_character_,
        TRUE ~ NA_character_
      )

      out$sequence_explanation <- combine_explanations(
        selected_explanation,
        discard_explanation
      )
    } else {
      out$sequence_explanation <- discard_explanation
    }
  }

  if ("match_status" %in% names(out)) {
    out$database_match_explanation <- explain_code_string(
      out$match_status,
      excel_explanation_dictionary("match_status")
    )
  }

  out
}


order_final_calls_columns <- function(df) {
  df <- add_excel_explanations(df)

  priority_columns <- c(
    "sample_id",
    "source_file",
    "TargetSNP",
    "zygosity",
    "call_explanation",
    "Haplotype_1",
    "Haplotype_2",
    "TargetSNP_allele_1",
    "TargetSNP_allele_2",
    "txt_value",
    "txt_target_value",
    "Read_1",
    "Read_2",
    "allele_balance",
    "extra_signal_ratio",
    "n_sequences",
    "Polymorphic_sites",
    "Position_GRCh38",
    "included_in_txt",
    "review_level",
    "review_explanation"
  )

  df |>
    dplyr::select(
      dplyr::any_of(priority_columns)
    )
}


order_review_columns <- function(df) {
  df <- add_excel_explanations(df)

  priority_columns <- c(
    "sample_id",
    "source_file",
    "TargetSNP",
    "review_level",
    "review_explanation",
    "zygosity",
    "call_explanation",
    "txt_explanation",
    "Haplotype_1",
    "Haplotype_2",
    "TargetSNP_allele_1",
    "TargetSNP_allele_2",
    "ForenSeq_sequence_1",
    "ForenSeq_sequence_2",
    "txt_value",
    "txt_target_value",
    "Read_1",
    "Read_2",
    "allele_balance",
    "extra_signal_ratio",
    "extra_signal_above_threshold",
    "n_sequences",
    "Polymorphic_sites",
    "Position_GRCh38",
    "final_Haplotype_1",
    "final_Haplotype_2",
    "final_TargetSNP_allele_1",
    "final_TargetSNP_allele_2"
  )

  df |>
    dplyr::select(
      dplyr::any_of(priority_columns)
    )
}


order_all_sequences_columns <- function(df) {
  df <- add_excel_explanations(df)

  priority_columns <- c(
    "sample_id",
    "source_file",
    "TargetSNP",
    "sequence_rank",
    "selected_for_call",
    "selected_copies",
    "sequence_explanation",
    "ForenSeq_sequence",
    "Read",
    "Detected Bases",
    "Haplotype",
    "TargetSNP_allele",
    "Polymorphic_sites",
    "Position_GRCh38",
    "match_status",
    "database_match_explanation"
  )

  df |>
    dplyr::select(
      dplyr::any_of(priority_columns)
    )
}


drop_excel_only_columns <- function(df) {
  df |>
    dplyr::select(
      -dplyr::any_of(c("read_3", "sequence_key"))
    )
}

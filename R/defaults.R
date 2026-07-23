#' Return the default iSNP loci order
#'
#' Returns the default loci order used to create the wide TXT file for
#' Familias. The first column in the TXT is always `sample_id`; the loci
#' returned by this function are appended after it.
#'
#' @return A character vector with iSNP locus names.
#' @export
default_loci_order <- function() {
  c(
    "rs10495407",
    "rs1294331",
    "rs1413212",
    "rs1490413",
    "rs560681",
    "rs891700",
    "rs1109037",
    "rs12997453",
    "rs876724",
    "rs907100",
    "rs993934",
    "rs1355366",
    "rs1357617",
    "rs2399332",
    "rs4364205",
    "rs6444724",
    "rs1979255",
    "rs2046361",
    "rs279844",
    "rs6811238",
    "rs13182883",
    "rs159606",
    "rs251934",
    "rs338882",
    "rs717302",
    "rs13218440",
    "rs1336071",
    "rs214955",
    "rs727811",
    "rs321198",
    "rs6955448",
    "rs737681",
    "rs917118",
    "rs10092491",
    "rs2056277",
    "rs4606077",
    "rs763869",
    "rs1015250",
    "rs10776839",
    "rs1360288",
    "rs1463729",
    "rs7041158",
    "rs3780962",
    "rs735155",
    "rs740598",
    "rs826472",
    "rs964681",
    "rs10488710",
    "rs1498553",
    "rs2076848",
    "rs901398",
    "rs10773760",
    "rs2107612",
    "rs2111980",
    "rs2269355",
    "rs2920816",
    "rs1058083",
    "rs1335873",
    "rs1886510",
    "rs354439",
    "rs1454361",
    "rs4530059",
    "rs722290",
    "rs873196",
    "rs1528460",
    "rs1821380",
    "rs8037429",
    "rs1382387",
    "rs2342747",
    "rs430046",
    "rs729172",
    "rs740910",
    "rs8078417",
    "rs938283",
    "rs9905977",
    "rs1024116",
    "rs1493232",
    "rs1736442",
    "rs9951171",
    "rs576261",
    "rs719366",
    "rs1005533",
    "rs1031825",
    "rs1523537",
    "rs445251",
    "rs221956",
    "rs2830795",
    "rs2831700",
    "rs722098",
    "rs914165",
    "rs1028528",
    "rs2040411",
    "rs733164",
    "rs987640"
  )
}

#' Return the default general allele balance threshold
#'
#' Returns the default general threshold used to decide whether the second
#' most abundant sequence should be accepted as a second haplotype.
#'
#' @return A numeric value.
#' @export
default_heterozygote_threshold <- function() {
  0.3
}

#' Return the default extra signal threshold
#'
#' Returns the threshold used to flag a relevant third sequence.
#'
#' @return A numeric value.
#' @export
default_extra_signal_threshold <- function() {
  0.25
}

#' Return the default locus-specific thresholds
#'
#' Returns the default locus-specific thresholds used when evaluating
#' heterozygous candidates. Loci not listed here use the general
#' heterozygote threshold.
#'
#' @return A named numeric vector.
#' @export
default_thresholds_by_locus <- function() {
  c(
    rs729172   = 0.1,
    rs10776839 = 0.1,
    rs1335873  = 0.1,
    rs338882   = 0.2,
    rs1493232  = 0.2,
    rs6955448  = 0.2
  )
}

#' Return the default minimum reads for homozygous duplicated calls
#'
#' Returns the minimum number of reads required to duplicate a single sequence
#' as a homozygous call.
#'
#' @return A numeric value.
#' @export
default_min_homozygote_reads <- function() {
  30
}


#' Return the default minimum reads for heterozygous calls
#'
#' Returns the minimum number of reads required for the second sequence to be
#' accepted as the second allele in a heterozygous call.
#'
#' @return A numeric value.
#' @export
default_min_heterozygous_reads <- function() {
  11
}

#' Return the default lower review balance margin
#'
#' Returns the proportional lower margin below the allele balance threshold
#' used to flag borderline calls for manual review.
#'
#' @return A numeric value.
#' @export
default_review_balance_lower_margin <- function() {
  0.15
}

#' Return the default upper review balance margin
#'
#' Returns the proportional upper margin above the allele balance threshold
#' used to flag borderline calls for manual review.
#'
#' @return A numeric value.
#' @export
default_review_balance_upper_margin <- function() {
  0.15
}

#' Return the default low homozygote review multiplier
#'
#' Returns the multiplier applied to `min_homozygote_reads` to flag homozygous
#' calls close to the minimum read threshold for manual review.
#'
#' @return A numeric value.
#' @export
default_review_low_homozygote_multiplier <- function() {
  1.5
}

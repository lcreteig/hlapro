#' Split an HLA typing string into alleles
#'
#' @description Takes in a space-separated HLA typing string and splits it into
#'   its constituent loci and alleles ("A_1", "A_2", "DRB1_1").
#'
#'   `extract_alleles_str()` takes in a single string, and returns a named
#'   character vector of alleles.
#'
#'   `extract_alleles_df()` takes in a data frame, where one column contains the
#'   typing string, and returns the same data frame along with a new column for
#'   each allele.
#'
#' @param string String, space-separated HLA typing.
#' @param df A data frame.
#' @param col_typing The name of the column in `df` that contains a
#'   space-separated HLA typing `string` for each row.
#' @param loci A string or character vector with the loci you are interested in.
#'   Only these alleles will be returned. Defaults to all. `DRB.` is used for
#'   DRB3, DRB4, and DRB5.
#'
#' @return Either a character vector or a data frame with the named alleles.
#' @export
#'
#' @examples
#' typing <- "A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53"
#' extract_alleles_str(typing, loci = "A")
#' extract_alleles_str(typing)
#'
#' df <- tidyr::tibble(typing = typing)
#' extract_alleles_df(df, "typing", loci = c("A", "B", "C"))
#'
#' # Can also handle newer nomenclature
#' extract_alleles_str("DQB1*03:01 DQB1*05:01 DRB1*04:AMR",
#'                     loci = c("DRB1", "DQB1"))

extract_alleles_str <- function(
    string,
    loci = c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1", "DRB.")) {
  loci <- rlang::arg_match(loci, multiple = TRUE)
  # TODO: document

  # get rid of any leading/trailing/double spaces
  string <- stringr::str_squish(string)

  extract_str <- function(locus, string) {
    pattern <- build_pattern(locus)
    allele_names <- stringr::str_c(locus, c("_1", "_2"))

    # return named list (locus_allele) with matches
    rlang::set_names(
      stringr::str_match(string, pattern)[, 2:3],
      allele_names
    )
  }

  purrr::map(loci, \(x) extract_str(x, string)) |>
    purrr::flatten_chr()
}

#' @rdname extract_alleles_str
#' @export
extract_alleles_df <- function(
    df,
    col_typing,
    loci = c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1", "DRB.")) {
  loci <- rlang::arg_match(loci, multiple = TRUE)
  # TODO: document

  # get rid of any leading/trailing/double spaces
  df[col_typing] <- stringr::str_squish(df[col_typing])

  extract_df <- function(locus, df, col_typing) {
    pattern <- build_pattern(locus)
    allele_names <- stringr::str_c(locus, c("_1", "_2"))

    tidyr::extract(df, tidyr::all_of(col_typing),
      into = allele_names,
      regex = pattern,
      remove = FALSE
    )
  }

  # Split typing for each locus into two columns; combine resulting data frames
  purrr::map(loci, \(x) extract_df(x, df, col_typing)) |>
    purrr::reduce(dplyr::full_join)
}

build_pattern <- function(locus) {
  regexps <- list(
    # don't match if locus is preceded by : or other capital letter. This
    # prevents matching in NMDP Multiple Allele Codes (":AABJE") or other loci
    # ("B" in "DRB1") or other prefixes (e.g. "A" in "HLA-")
    neg = "(?<![:A-Z])",
    allele = r"((?:{locus}\*?(\S+)))", # match locus and following non-spaces
    loci = list(
      A = "A",
      B = "B(?![Ww])", # B cannot be followed by "W" (pubic)
      C = "Cw?", # in serological notation, C is followed by "w"
      DPB1 = "DPB1",
      DQA1 = "DQA1", # DQA1 should always follow this format
      DQB1 = "DQ(?!A)(?:B1)?", # could be "DQB1" but also "DQ" with no A after
      DRB1 = "DR(?!A|B[2-9]|5[1-3])(?:B1)?", # either "DRB1"
      # or "DR" (excluding DRA, DR51-53)
      DRB. = r"(DR(?:(?=5[1-3])|B(?=[3-5](?!\*Neg))))" # either DR51/52/53 or
      # DRB3/4/5. Sometimes all 3 are specified, with 1-2 having suffix "*Neg".
      # These should be skipped, in order to still retrieve full typing
    ),
    # don't capture other alleles in between: i.e. any that aren't current locus
    inbetween = r"((?:\s(?!{locus})\S+)*\s?)"
  )

  # build the regular expression from its sub-parts
  stringr::str_c(
    regexps$neg,
    stringr::str_glue(regexps$allele, locus = regexps$loci[locus]),
    stringr::str_glue(regexps$inbetween, locus = regexps$loci[locus]),
    stringr::str_glue(regexps$allele, locus = regexps$loci[locus]),
    "?" # second allele is optional
  )
}

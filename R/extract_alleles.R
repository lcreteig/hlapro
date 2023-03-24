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

extract_alleles <- function(
    df,
    col_typing,
    locus = c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1", "DRB.")) {
  locus <- rlang::arg_match(locus)
  # TODO: implement for character input via stringr's str_match() or str_extract()
  # TODO: vectorize locus arg

  # get rid of any leading/trailing/double spaces
  df[col_typing] <- stringr::str_squish(df[col_typing])

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

  column_names <- stringr::str_c(locus, c("_1", "_2"))
  # build the regular expression from its sub-parts
  pattern <- stringr::str_c(
    regexps$neg,
    stringr::str_glue(regexps$allele, locus = regexps$loci[locus]),
    stringr::str_glue(regexps$inbetween, locus = regexps$loci[locus]),
    stringr::str_glue(regexps$allele, locus = regexps$loci[locus]),
    "?" # second allele is optional
  )

  tidyr::extract(df, tidyr::all_of(col_typing),
    into = column_names,
    regex = pattern,
    remove = FALSE
  )
}

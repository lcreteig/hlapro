extract_alleles <- function(df, col_typing, locus = c("A")) {
  locus <- rlang::arg_match(locus)
  #TODO: test named groups (cf. stringr::str_match() in tidyr::extract())

  regexps <- list(
    # don't match if locus is preceded by : or other capital letter. This
    # prevents matching in NMDP Multiple Allele Codes (":AABJE") or other loci
    # ("B" in "DRB1") or other prefixes (e.g. "A" in "HLA-")
    neg = "(?<![:A-Z])",
    allele = r"((?:{locus}\*?(\S+)))", # match locus and all following non-spaces
    loci = list(
      A = "A"),
    # don't capture other alleles in between: i.e. any that aren't the current locus
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

  tidyr::extract(df, {{ col_typing }},
                 into = column_names,
                 regex = pattern,
                 remove = FALSE
  )
}

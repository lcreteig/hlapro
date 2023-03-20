extract_alleles <- function(df, col_typing, locus = c("A")) {
  locus <- rlang::arg_match(locus)

  regexps <- list(
    anything = ".*", # at the start
    neg_nmdp = "(?<![:A-Z])", # don't match NMDP Multiple Allele codes
    allele = r"(\S+)", # allele
    loci = list(
      A = r"(A\*?)"
                 ),
    inbetween = list( # other HLAs inbetween
      A = r"((?:\s(?!A)\S+)*\s?)"
                     ))

  pattern <- c(regexps$anything,
    regexps$neg_nmdp, regexps$loci[[locus]],
    regexps$allele,
    regexps$inbetween[[locus]],
    regexps$neg_nmdp, regexps$loci[[locus]],
    regexps$allele
  )

  names(pattern)[4] <- paste0(locus,"_1")
  names(pattern)[8] <- paste0(locus,"_2")

  tidyr::separate_wider_regex(df, {{ col_typing }},
    patterns = pattern,
    too_few = "align_start",
    cols_remove = FALSE
  )
}

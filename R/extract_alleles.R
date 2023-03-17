extract_alleles <- function(df, col_typing) {
  # TODO: vectorize with locus argument, and named list that can be indexed with it
  anything <- ".*" # at the start
  neg_nmdp <- "(?<![:A-Z])" # don't match NMDP Multiple Allele codes
  a <- r"(A\*?)" # locus
  allele <- r"(\S+)" # allele
  inbetween <- r"((?:\s(?!A)\S+)*\s?)" # other HLAs inbetween

  pattern <- c(anything,
    neg_nmdp, a,
    A_1 = allele,
    inbetween,
    neg_nmdp, a,
    A_2 = allele
  )

  tidyr::separate_wider_regex(df, {{ col_typing }},
    patterns = pattern,
    too_few = "align_start",
    cols_remove = FALSE
  )
}

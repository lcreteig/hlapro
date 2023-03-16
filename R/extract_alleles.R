extract_alleles <- function(df, col_typing) {
  # TODO: vectorize with locus argument, and named list that can be indexed with it
  match_anything <- ".*"
  neg_nmdp <- "(?<![:A-Z])"
  match_allele <- r"(A\S+)"
  inbetween_alleles <- r"((?:\s(?!A)\S+)*\s?)"

  pattern <- c(
    match_anything,
    neg_nmdp,
    A_1 = match_allele,
    inbetween_alleles,
    neg_nmdp,
    A_2 = match_allele
  )

  tidyr::separate_wider_regex(df, {{ col_typing }},
    patterns = pattern,
    too_few = "align_start",
    cols_remove = FALSE
  )
}

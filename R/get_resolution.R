get_resolution <- function(allele) {
  # N.B. assumes all ambiguities are intermediate, even when the ambiguity is >2 field codes
  dplyr::case_when(
    # "*" followed by capital letter (but not XX), or digits and then a slash
    stringr::str_detect(
      allele,
      r"(\*\d+:?(?!XX)([A-Z]|\d+\/))"
    ) ~ "intermediate",
    # "*" followed by 4 digits with optional semicolon in between
    # TODO: return field code with high? (but some already at allelic resolution with 4 fields)
    stringr::str_detect(
      allele,
      r"(\*\d{2}:?\d{2})"
    ) ~ "high",
    is.na(allele) ~ NA_character_,
    .default = "low"
  )
}

get_resolution <- function(allele) {
  dplyr::case_when(
    # "*" followed by capital letter or digits and then a slash
    stringr::str_detect(allele, r"(\*\d+:?([A-Z]|\d+\/))") ~ "intermediate",
    # "*" followed by 4 digits with optional semicolon in between
    stringr::str_detect(allele, r"(\*\d{2}:?\d{2})") ~ "high",
    is.na(allele) ~ NA_character_,
    .default = "low"
  )
}

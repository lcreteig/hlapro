#' Determine the level of resolution of an HLA-allele
#'
#' `get_resolution()` takes in a string or character vector of HLA alleles, and
#' returns their resolution as either low, intermediate, or high.
#'
#' Resolution is defined as follows:
#'
#' - **Low**: a typing at the allele group level, such as
#'    - anything using the serological/antigen nomenclature, e.g `A2`
#'    - an XX code representing an allele group, e.g. `A*02:XX`
#'  - **Intermediate**: any typing with ambiguities, such as
#'    - a list using the forward slash notation, e.g. `HLA-A*23:26/HLA-A*23:39`
#'    - a typing with multiple allele codes, e.g. `A*01:AABJE`
#'  - **High**: any typing with more than 2 field codes, e.g.
#'    - `A*24:09` (minimum)
#'    - `HLA-A*24:02:01:02L` (maximum)
#'
#' @param allele A string or character vector with (an) HLA allele(s).
#'
#' @return A string or character vector of the same length as `allele`,
#'   with `"low"`, `"intermediate"`, or `"high"` for each element.
#' @export
#'
#' @examples
#' get_resolution("A2") # low
#' get_resolution("A*01:AABJE") # intermediate
#' get_resolution("A*24:09") # high
#'
#' # also works with character vectors, or in a data frame
#' allele_vec <- c("A2", "A*01:AABJE", "B*42:08")
#' get_resolution(allele_vec)
#'
#' tidyr::tibble(alleles = allele_vec) |>
#'   dplyr::mutate(resolution = get_resolution(allele_vec))
#'
get_resolution <- function(allele) {
  # N.B. assumes all ambiguities are intermediate, even when >2 field codes
  dplyr::case_when(
    # "*" followed by capital letter (but not XX), or digits and then a slash
    stringr::str_detect(
      allele,
      r"(\*\d+:?(?!XX)([A-Z]|\d+\/))"
    ) ~ "intermediate",
    # "*" followed by 4 digits with optional semicolon in between
    # TODO: return field code with high?
    stringr::str_detect(
      allele,
      r"(\*\d{2,3}:?\d{2,3})"
    ) ~ "high",
    is.na(allele) ~ NA_character_,
    .default = "low"
  )
}

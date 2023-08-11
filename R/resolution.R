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
#' @param extended When `TRUE`, the returned resolution also contains:
#'  - for high resolution: the number of fields ("second field", "third field",
#'    "fourth field")
#'  - for low resolution: whether the allele is a serological/molecular
#'    split/broad or
#'    [associated antigen](https://hla.alleles.org/antigens/broads_splits.html)
#'
#' @return A string or character vector of the same length as `allele`,
#'   with `"low"`, `"intermediate"`, or `"high"` for each element.
#' @export
#'
#' @examples
#' get_resolution("A2") # low
#' get_resolution("A2", extended = TRUE) # low - broad
#' get_resolution("A*01:AABJE") # intermediate
#' get_resolution("A*24:09") # high
#' get_resolution("A*24:09", extended = TRUE) # high - second field
#'
#' # also works with character vectors, or in a data frame
#' allele_vec <- c("A2", "A*01:AABJE", "B*42:08")
#' get_resolution(allele_vec)
#'
#' tidyr::tibble(alleles = allele_vec) |>
#'   dplyr::mutate(resolution = get_resolution(allele_vec))
#'
get_resolution <- function(allele, extended = FALSE) {
  check_bool(extended)
  res <- dplyr::case_when(
    # "*" followed by capital letter (but not XX), or digits and then a slash
    stringr::str_detect(
      allele,
      r"(\*\d+:?(?!XX$)([A-Z]|\d+\/))"
    ) ~ "intermediate",
    # "*" followed by 4 digits with optional semicolon in between
    stringr::str_detect(
      allele,
      r"(\*\d{2,3}:?\d{2,3})"
    ) ~ "high",
    is.na(allele) ~ NA_character_,
    .default = "low"
  )

  if (!extended) {
    return(res)
  }

  ord_seq <- c("first", "second", "third", "fourth")
  dplyr::case_when(
    res == "low" & is_associated(allele)
    ~ "serology - associated",
    res == "low" & !is_serology(allele) & !is.na(get_split(allele))
    ~ "molecular - split",
    res == "low" & is_serology(allele) & !is.na(get_split(allele))
    ~ "serology - split",
    res == "low" & !is_serology(allele) & !is.na(get_broad(allele))
    ~ "molecular - broad",
    res == "low" & is_serology(allele) & !is.na(get_broad(allele))
    ~ "serology - broad",
    res == "high" ~ stringr::str_c(
      res,
      " - ",
      ord_seq[get_n_fields(allele)],
      " field"
    ),
    .default = res
  )
}

#' Count how many fields an allele has
#'
#' Currently an allele can have 1-4 fields (see the [HLA nomenclature website]
#' (https://hla.alleles.org/nomenclature/naming.html)). These must be delimited
#' with a colon.
#'
#' @inheritParams get_resolution
#'
#' @return An integer with the amount of fields.
#' @noRd
#'
#' @examples
#' allele_vec <- c("A2", "A*01:AABJE", "B*42:08", "C*01:02:01:26")
#' get_n_fields()
#'
get_n_fields <- function(allele) {
  stringr::str_count(stringr::str_remove(allele, r"(\/.*)"), ":") + 1
}

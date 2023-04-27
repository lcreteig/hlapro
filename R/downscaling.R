#' Get serological equivalents of an HLA-allele
#'
#' `etrl_lookup()` takes in a string or character vector of HLA alleles, and
#' returns their serological equivalents as defined in the ETRL HLA tables.
#'
#' This function uses the EuroTransplant Reference Laboratory HLA ([etrl_hla])
#' tables to do the lookup. These tables define several alleles for each protein
#' (e.g. `HLA-A*01:01`, `HLA-A*01:02`), and their serological equivalents. All
#' others are grouped into an XX code (`HLA-A*01:XX`), which is also mapped to
#' a serological equivalent.
#'
#' @section Workings:
#'
#' All entered alleles will first be reduced to the two-field level (see
#' [reduce_to_nth_field()]). If this reduced allele occurs in the lookup table,
#' the corresponding rows are returned. If not, it will be converted into an XX
#' code, and that row will be returned.
#'
#' @section Exceptions:
#'
#' If the allele is already in serological notation (e.g. `A1`), the lookup
#' fails and empty rows are returned (though e.g. `A1` is of course present
#' in the other columns of the lookup table)
#'
#' If the allele has a suffix (e.g. `HLA-C*01:37N`), it has no serological
#' equivalent, and hence will also return empty rows.
#'
#' @inheritParams get_resolution
#'
#' @return A data frame with as many rows as there were elements in the input,
#'  and four columns:
#'  1. The allele that was looked up (or `NA` if it didn't exist)
#'  2. Its equivalent at the serological *split* level (or `NA` if it doesn't
#'    have one)
#'  3. Its equivalent at the serological *broad* level (or `NA` if it doesn't
#'    have one)
#'  4. Whether the allele has the public epitope (Bw4 or Bw6) (or `NA` if it
#'    doesn't)
#' @seealso [etrl_hla]: the lookup table that's used by and returned in this
#'  function
#' @export
#'
#' @examples
#' allele_vec <- c(
#'   "B15", "B*15:79N", "B*15:YETY",
#'   "B*15:01:16", "B*15:02", "B*15:85"
#' )
#' etrl_lookup(allele_vec)
etrl_lookup <- function(allele) {
  ids <- match(etrl_convert(allele), etrl_hla$Allele)
  etrl_hla[ids, ]
}

#' Retrieve broad-level equivalent of a split-level HLA-allele
#'
#' `get_broad()` takes in a string or character vector of HLA alleles. If they
#' are at the serological split-level, the corresponding broad-level allele is
#' looked up in [etrl_hla] and returned. If the input is not serology, or not a
#' split, `NA` is returned instead.
#'
#' @inheritParams get_resolution
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the corresponding broad if it exists (i.e., if `allele` was a
#'   serological split), or `NA` if none exists.
#' @export
#' @seealso
#'  - [get_public]: for looking up the public epitope of an allele
#'  - [etrl_hla]: the lookup table that's used in this function
#'  - [etrl_lookup]: for getting serological equivalents of other alleles
#'
#' @examples
#' get_broad("A24") # returns corresponding broad ("A9")
#' get_broad("A1") # is already a broad, returns `NA`
#' get_broad("A*24") # is not in serological notation, returns `NA`.
#' # use `etrl_lookup()` for these cases instead
#'
#' # Vectors also work:
#' get_broad(c("A24", "A23", "A1"))
get_broad <- function(allele) {
  unname(etrl_split_to_broad[allele])
}

#' Retrieve public epitope of serological HLA-allele
#'
#' `get_public()` takes in a string or character vector of HLA alleles. If they
#' are in serological notation (broad or split), the corresponding public
#' epitope (`Bw4` or `Bw6`) is looked up in [etrl_hla] and returned. If the
#' input is not serology, or the allele does not have the public epitope, `NA`
#' is returned instead.
#'
#' @inheritParams get_resolution
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the public epitope if it has one, or `NA` if none exists.
#' @export
#' @seealso
#'  - [get_broad]: for looking up the broad-level equivalent of a split allele
#'  - [etrl_hla]: the lookup table that's used in this function
#'  - [etrl_lookup]: for getting serological equivalents of other alleles
#'
#' @examples
#' get_public("B14") # has the epitope ("Bw6")
#' get_public("A1") # does not, returns `NA`
#' get_public("B*14") # is not in serological notation, returns `NA`.
#' # use `etrl_lookup()` for these cases instead
#'
#' # Vectors also work:
#' get_public(c("B14", "B63", "A1"))
get_public <- function(allele) {
  unname(etrl_public[allele])
}

#' Truncate an HLA-allele to a lower field
#'
#' `reduce_to_nth_field()` takes in a string or character vector of HLA alleles,
#' and reduces their resolution to the specified number of fields.
#'
#' N.B. This does not do a "proper" reduction, but simply truncates the string
#' to the last field delimiter (":"). It does *not* take into account G or P
#' groups, nor null/alternatively expressed alleles. These cannot be simply
#' be reduced to a lower field, but `reduce_to_nth_field()` simply strips off
#' these suffixes and does it anyway.
#'
#' @inheritParams get_resolution
#' @param n An integer. Can technically be anything, but the function will
#'   normally only do something with integers between 1 and 3, as these are the
#'   only possible number of fields in the current nomenclature (4 fields is the
#'   maximum and therefore cannot be further reduced).
#'
#' @return A string or character vector of the same length as `allele`, with the
#'   truncated alleles in each element. Alleles already at the desired level of
#'   resolution will be returned unchanged.
#' @export
#' @keywords internal
#'
#' @examples
#' reduce_to_nth_field("B*07:14:01", 2)
#'
#' # These won't do anything, as the resolution is equal to `n`:
#' reduce_to_nth_field("A2", 1)
#' reduce_to_nth_field("B*07:14:01", 3)
#'
#' allele_vec <- c("A*01:AABJE", "B*42:08", "C*01:02:01:26")
#' reduce_to_nth_field(allele_vec, 1)
reduce_to_nth_field <- function(allele, n) {
  allele <- remove_suffixes_groups(allele)
  # logical index of all alleles to be reduced
  res_idx <- get_n_fields(allele) > n & !is.na(allele)

  if (!any(res_idx, na.rm = TRUE)) {
    return(allele)
  }

  # Get list of field separator locations in each allele string
  field_locs <- stringr::str_locate_all(allele[res_idx], ":")
  ends <- purrr::map_int(field_locs, n) - 1 # end of nth field

  # Keep only the part of the typing up until the nth ":"
  replace(allele, res_idx, stringr::str_sub(allele[res_idx], 1, ends))
}

#' Strip broad in `split`(`broad`) notation
#'
#' @description
#' Removes parentheses and their contents, in e.g. `A24(A9)`. Because this
#' notation is not accepted by [validate_allele()], and broads can always be
#' added back with [get_broad()].
#'
#' TODO: should at some point be subsumed in a general `hla_clean()` type
#' function
#'
#' @inheritParams get_resolution
#' @param check_result If TRUE, returns input as is if the result of removing
#' the parentheses does not result in an existing split in [etrl_hla]. Sometimes
#' splits and broads are reversed (e.g. `A9(A24)`), in which case we don't want
#' to proceed.
#'
#' @return a vector of same length as allele
#' @keywords internal
#' @export
#'
#' @examples
#' strip_broad("A24(A9)")
strip_broad <- function(allele, check_result = FALSE) {
  allele_stripped <- stringr::str_remove(allele, r"(\(.*\))")

  if (check_result) {
    # return input as is if the result is not an existing split
    ifelse(is_split(allele), allele_stripped, allele)
  } else {
    allele_stripped
  }
}

etrl_convert <- function(allele) {
  allele <- remove_hla_prefix(allele)

  allele_f2 <- reduce_to_nth_field(allele, 2)

  # If allele in ETRL table, use as is, otherwise make into xx code
  ifelse(allele_f2 %in% etrl_hla$Allele, allele_f2, make_xx(allele_f2)) |>
    replace(has_suffix(allele), "") |> # suffixes cannot be reduced
    ifelse(is_serology(allele), allele, no = _) # return serology as is
}

make_xx <- function(allele) {
  stringr::str_c(reduce_to_nth_field(allele, 1), ":XX")
}

is_serology <- function(allele) {
  !stringr::str_detect(allele, r"([\*:])")
}

has_suffix <- function(allele) {
  stringr::str_detect(allele, "[NLSCAQ]$")
}

remove_suffixes_groups <- function(allele) {
  stringr::str_remove(allele, "[NLSCAQPG]$")
}

remove_hla_prefix <- function(allele) {
  stringr::str_remove(allele, "^HLA-")
}

is_split <- function(allele) {
  allele %in% names(etrl_split_to_broad)
}

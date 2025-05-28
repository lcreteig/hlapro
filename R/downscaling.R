#' Retrieve serological equivalents of an HLA-allele
#'
#' `get_serology()` takes in a string or character vector of HLA alleles. The
#' corresponding split-level (if it exists) or broad-level (if no split exists)
#' allele is looked up in [etrl_hla] and returned. If no such alleles exist,
#' `NA` is returned instead.
#'
#' @details
#'
#' This function uses the EuroTransplant Reference Laboratory HLA ([etrl_hla])
#' tables to do the lookup. These tables define several alleles for each protein
#' (e.g. `HLA-A*01:01`, `HLA-A*01:02`), and their serological equivalents. All
#' others are grouped into an XX code (`HLA-A*01:XX`), which is also mapped to
#' a serological equivalent.
#'
#' ## Workings
#'
#' All entered alleles will first be reduced to the two-field level (see
#' [reduce_to_nth_field()]). If this reduced allele occurs in the lookup table,
#' the corresponding rows are returned. If not, it will be converted into an XX
#' code, and that row will be returned.
#'
#' If the allele is already at the serological broad- or split-level, the lookup
#' will be performed using those respective columns in [etrl_hla].
#'
#' ## Exceptions
#'
#' If the allele has a suffix (e.g. `HLA-C*01:37N`), it has no serological
#' equivalent, and hence will also return nothing.
#'
#' @inheritParams get_resolution
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the corresponding serology if it exists, or `NA` if none exists.
#' @export
#' @seealso
#'  - [get_broad]: for looking up the broad-level equivalent of a split allele
#'  - [get_split]: for looking up the serological split-level equivalent of an
#'    allele
#'  - [get_public]: for looking up the public epitope of an allele
#'  - [etrl_hla]: the lookup table that's used in this function
#'
#' @examples
#' get_serology("A24") # is a serological split; returns itself ("A24")
#' get_serology("A9") # is a serological broad; returns itself ("A9")
#'
#' get_serology("A*01:01:01:50") # has a broad-level equivalent only ("A1")
#' get_serology("A*23:01:01:11") # has a split equivalent ("A23")
#' # as well as a broad ("A9"); only the former is returned
#'
#' # Vectors also work:
#' get_serology(c("A24", "A*01:XX", "B*15:15"))
get_serology <- function(allele) {
  # if not already a broad/split, convert to broad
  ifelse(is_broad(allele) | is_split(allele), allele, get_broad(allele)) |>
    # if it has a split, get that, otherwise leave unchanged
    ifelse(!is.na(get_split(allele)), get_split(allele), no = _)
}

#' Retrieve broad-level serological equivalent of an HLA-allele
#'
#' `get_broad()` takes in a string or character vector of HLA alleles. The
#' corresponding broad-level allele is looked up in [etrl_hla] and returned.
#' If no such allele exists, `NA` is returned instead.
#'
#' @inherit get_serology details
#' @inheritParams get_resolution
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the corresponding broad if it exists, or `NA` if none exists.
#' @export
#' @seealso
#'  - [get_split]: for looking up the serological split-level equivalent of an
#'    allele
#'  - [get_serology]: for looking up the split, or broad if none exists
#'  - [get_public]: for looking up the public epitope of an allele
#'  - [etrl_hla]: the lookup table that's used in this function
#'
#' @examples
#' get_broad("A24") # returns corresponding broad ("A9")
#' get_broad("A9") # is already a broad, returns itself ("A9")
#'
#' # these alleles in modern nomenclature also all return "A9"
#' get_broad("A*24")
#' get_broad("A*24:XX")
#' get_broad("A*24:02:01:102")
#'
#' # Vectors also work:
#' get_broad(c("A24", "A23", "A1"))
get_broad <- function(allele) {
  # if it's a broad return as is, else assume it's a split and lookup the broad
  ifelse(is_broad(allele), allele, unname(etrl_split_to_broad[allele])) |>
    # if it's not serology, try to lookup a two-field allele
    ifelse(is_serology(allele),
      yes = _,
      etrl_lookup(allele)$`ET MatchDeterminantBroad`
    )
}

#' Retrieve split-level serological equivalent of an HLA-allele
#'
#' `get_split()` takes in a string or character vector of HLA alleles. The
#' corresponding split-level allele is looked up in [etrl_hla] and returned.
#' If no such allele exists, `NA` is returned instead.
#'
#' @inherit get_serology details sections
#' @inheritParams get_resolution
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the corresponding split if it exists, or `NA` if none exists.
#' @export
#' @seealso
#'  - [get_broad]: for looking up the serological broad-level equivalent of an
#'    allele
#'  - [get_serology]: for looking up the split, or broad if none exists
#'  - [get_public]: for looking up the public epitope of an allele
#'  - [etrl_hla]: the lookup table that's used in this function
#'
#' @examples
#'
#' get_split("A24") # is already a split, so returns itself ("A24")
#' get_split("A9") # is a broad, so returns `NA`
#'
#' get_split("B*14:01") # returns corresponding split ("B64")
#' get_split("B*14:02") # returns corresponding split ("B65")
#' get_split("B*14:03") # no split is defined for this allele; returns `NA`
#'
#' # Vectors also work:
#' get_split(c("A24", "B*14:01", "B*14:03"))
get_split <- function(allele) {
  ifelse(is_split(allele),
    allele,
    etrl_lookup(allele)$`ET MatchDeterminantSplit`
  )
}

#' Retrieve public epitope of serological HLA-allele
#'
#' `get_public()` takes in a string or character vector of HLA alleles. The
#' corresponding public epitope (`Bw4` or `Bw6`) is looked up in [etrl_hla] and
#' returned. If the input allele does not have the public epitope, `NA`
#' is returned instead.
#'
#' @inherit get_serology details sections
#' @inheritParams get_resolution
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the public epitope if it has one, or `NA` if none exists.
#' @export
#' @seealso
#'  - [get_broad]: for looking up the broad-level equivalent of a split allele
#'  - [get_split]: for looking up the serological split-level equivalent of an
#'    allele
#'  - [get_serology]: for looking up the split, or broad if none exists
#'  - [etrl_hla]: the lookup table that's used in this function
#'
#' @examples
#' get_public("B14") # has the epitope ("Bw6")
#' get_public("B*15:12") # has "Bw6"
#' get_public("B*15:13") # but this alleles has "Bw4"
#' get_public("B*15:XX") # hence this is ambiguous (returns `NA`)
#'
#' get_public("A1") # does not have the epitope; returns `NA`
#'
#' # Vectors also work:
#' get_public(c("B14", "B63", "A1"))
get_public <- function(allele) {
  ifelse(is_serology(allele),
    unname(etrl_public[allele]),
    etrl_lookup(allele)$Public
  )
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
  # remove any null/alternative expression or group suffixes
  allele <- ifelse(has_suffix(allele) | is_group(allele),
    remove_suffixes_groups(allele),
    allele
  )
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

#' Swap a pair of HLA-alleles to match another
#'
#' `reorder_alleles()` takes in two length-2 character vectors of HLA alleles,
#'  and if they're not in the same order, reorders the 2nd to match the 1st.
#'  This can be useful when you have multiple sources of HLA typings for one
#'  individual that you're trying to match up.
#'
#'  If the vectors contain alleles in different nomenclatures / of different
#'  resolutions, the order is determined based on scaling down the resolution.
#'  First it tries to scale down to the 2-field level; if still all alleles
#'  differ from each other, it tries the serological split level, and finally
#'  the serological broad level.
#'
#' @param in_order,to_order A character vector of two HLA alleles of a certain
#'  locus. The order of `to_order` will be reversed, to match `in_order`, if
#'  necessary.
#'
#' @return A vector with the same alleles as `to_order`, possible reordered.
#' @keywords internal
#' @export
#'
#' @examples
#' reorder_alleles(in_order = c("A1", "A2"), to_order = c("A2", "A1"))
#' # can accommodate missing/homozygous typings
#' reorder_alleles(in_order = c("A1", NA), to_order = c("A2", "A1"))
#' # still works if the typings are in a different format/resolution
#' reorder_alleles(in_order = c("A1", "A2"), to_order = c("A*02:01", "A*01:01"))
reorder_alleles <- function(in_order, to_order) {
  to_order_orig <- to_order

  if (all(is.na(in_order)) || all(is.na(to_order))) {
    return(to_order_orig) # nothing to do
  }

  # scale down further and further if necessary
  downscaling_levels <- c(
    get_broad,
    get_serology,
    \(x) reduce_to_nth_field(x, 2)
  )
  ii <- 3
  # if all unique alleles, they're not in same format: scale down
  while ((length(unique(c(in_order, to_order))) > 3) && ii != 0) {
    in_order <- downscaling_levels[[ii]](in_order)
    to_order <- downscaling_levels[[ii]](to_order)
    ii <- ii - 1
  }

  # use identical() to allow for NA comparisons
  pair_1_identical <- identical(in_order[1], to_order[2])
  pair_2_identical <- identical(in_order[2], to_order[1])
  if (pair_1_identical || pair_2_identical) { # if one is in reverse order
    return(rev(to_order_orig)) # reverse the order
  }
  to_order_orig
}

#' Get rows with serological equivalents of an HLA allele from ETRL HLA table
#'
#' `etrl_lookup()` takes in a string or character vector of HLA alleles, and
#' returns their serological equivalents as defined in the ETRL HLA tables.
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
#' @keywords internal
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

is_v2 <- function(allele) {
  # 2 digits directly followed by 2 or more letters/numbers
  v2_pattern <- r"(\*\d{2}([A-Z0-9]{2,}))"
  # and v2s have no colons anywhere
  stringr::str_detect(allele, v2_pattern) & !stringr::str_detect(allele, ":")
}

has_suffix <- function(allele) {
  stringr::str_detect(allele, r"(\d[NLSCAQ]$)")
}

is_group <- function(allele) {
  stringr::str_detect(allele, r"(\d[PG]$)")
}

is_ambiguous <- function(allele) {
  stringr::str_detect(allele, c("\\/")) & stringr::str_detect(allele, c("\\*"))
}

remove_suffixes_groups <- function(allele) {
  stringr::str_remove(allele, "[NLSCAQPG]$")
}

has_hla_prefix <- function(allele) {
  stringr::str_detect(allele, "^HLA-")
}

add_hla_prefix <- function(allele) {
  ifelse(has_hla_prefix(allele), allele, stringr::str_c("HLA-", allele))
}

remove_hla_prefix <- function(allele) {
  stringr::str_remove(allele, "^HLA-")
}

is_split <- function(allele) {
  allele %in% names(etrl_split_to_broad)
}

is_broad <- function(allele) {
  allele %in% unique(etrl_hla$`ET MatchDeterminantBroad`)
}

is_associated <- function(allele) {
  is_serology(allele) & stringr::str_detect(allele, r"([A-Z]-?\d{3,4})")
}

is_public <- function(allele) {
  allele %in% c("Bw4", "Bw6")
}

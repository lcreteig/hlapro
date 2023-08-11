#' Determine mismatched antigens from typings
#'
#' `get_mismatches()` returns mismatched HLAs given a donor and recipient HLA
#' typing. That is, those antigens that occur in the *donor* but *not* the
#' recipient typing.
#'
#' @param donor,recipient A pair of strings containing the donor and recipient
#'   HLA typing, respectively. Each antigen must be separated by a space.
#'
#' @return A string with the mismatched antigens (separated by spaces).
#' @export
#'
#' @examples
#' donor <- "A1 A2"
#' recipient <- "A1 A2"
#' get_mismatches(donor, recipient)
#' #> ""
#' get_mismatches("A1 A2", "A3 A4")
#' #> "A1 A2"
#' get_mismatches("A1 A2 B5", "A3 A4 B5 B12")
#' #> "A1 A2"
get_mismatches <- function(donor, recipient) {
  stopifnot(
    "donor typing is not a single string" = length(donor) == 1,
    "recipient typing is not a single string" = length(recipient) == 1
  )
  setdiff(
    # separate alleles by spaces
    stringr::str_split_1(donor, stringr::boundary("word")),
    stringr::str_split_1(recipient, stringr::boundary("word"))
  ) |>
    stringr::str_flatten(collapse = " ") # put back into single string
}

get_mismatches <- function(donor, recipient) {
  setdiff(
    # separate alleles by spaces
    stringr::str_split_1(donor, stringr::boundary("word")),
    stringr::str_split_1(recipient, stringr::boundary("word"))
  ) |>
    stringr::str_flatten(collapse = " ") # put back into single string
}

#' Load IPD-IMGT/HLA database
#'
#' Use the Python [`py-ard` package](https://github.com/nmdp-bioinformatics/)
#' developed by the NMDP to initialize the
#' [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/).
#'
#' @param data_dir Path to the folder where the SQLite database (~600MB) is or
#'  should be stored. If no database exists on disk, it will be downloaded.
#' @param imgt_version Four-digit version number of the
#'  [database release](https://github.com/ANHIG/IMGTHLA/releases) that should
#'  be used. Defaults to the most recent release ("Latest").
#'
#' @return A Python database connection object, which should be passed to other
#'  functions that make use of the database.
#' @export
#'
#' @examples
#' \dontrun{
#' ard <- db_initialize(data_dir = "~/ipd_db/", imgt_version = "3510")
#' }
db_initialize <- function(data_dir, imgt_version = "Latest") {
  pyard <- reticulate::import("pyard")
  pyard$init(imgt_version = imgt_version, data_dir = data_dir)
}

#' Print version of initialized IPD-IMGT/HLA database
#'
#' @param ard The Python database connection object created by [db_initialize()]
#'
#' @export
#' @seealso [db_initialize()]
#'
#' @examples
#' \dontrun{
#' ard <- db_initialize(data_dir = "~/ipd_db/")
#' db_get_version(ard)
#' }
db_get_version <- function(ard) {
  ard$get_db_version()
}

reduce_to_serology <- function(ard, allele) {
  purrr::map_chr(allele, \(x) ard$redux(x, "S"))
}

reduce_to_field2 <- function(ard, allele) {
  purrr::map_chr(allele, \(x) ard$redux(x, "U2"))
}

mac_lookup <- function(ard, allele) {
  purrr::map_chr(allele, ard$lookup_mac)
}

mac_expand <- function(ard, allele) {
  purrr::map_chr(allele, ard$expand_mac)
}

is_in_ipd_db <- function(ard, allele) {
  # ard$validate() throws an error if the allele is not in the db, instead of
  # returning "FALSE", so we need to wrap it
  safe_validate <- purrr::safely(ard$validate, otherwise = FALSE)
  purrr::map(allele, safe_validate) |>
    purrr::map("result") |>
    purrr::simplify()
}

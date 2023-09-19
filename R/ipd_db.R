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
#' ard <- db_initialize(data_dir = "~/ipd_db", imgt_version = "3510")
#' }
db_initialize <- function(data_dir, imgt_version = "Latest") {
  pyard <- reticulate::import("pyard")
  pyard$init(imgt_version = imgt_version, data_dir = data_dir)
}

#' Print version of initialized IPD-IMGT/HLA database
#'
#' Thin wrapper around the `ard.db_get_version()` method of the
#' [`py-ard`](https://github.com/nmdp-bioinformatics/py-ard) Python package.
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

#' Retrieve serological equivalents of HLA alleles
#'
#' `reduce_to_serology()` is a thin wrapper around the `ard.redux()` method
#' (Reduction Type `'S'`) of the
#' [`py-ard`](https://github.com/nmdp-bioinformatics/py-ard) Python package.
#'
#' @details
#' Uses the [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/) as
#' initialized with [db_initialize()]. This function will throw an error if any
#' of the alleles in the input do not exist in the database. Use
#' [is_in_ipd_db()] to safely check if the allele(s) exist(s).
#'
#' @inheritParams db_get_version
#' @inheritParams get_resolution
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the corresponding serology if it exists, or "" if none exists (e.g.
#'   for null alleles)
#' @export
#' @seealso [reduce_to_field2()] to reduce an allele to two-field resolution
#'
#' @examples
#' \dontrun{
#' ard <- db_initialize(data_dir = "~/ipd_db/")
#' reduce_to_serology(ard, "B*13:03")
#' #> "B13"
#' reduce_to_serology(ard, "B*13:04")
#' #> "B15/B21"
#' # Also works for vectors:
#' reduce_to_serology(ard, c("B*13:03", "B*13:04"))
#' #> "B13"     "B15/B21"
#' }
reduce_to_serology <- function(ard, allele) {
  purrr::map_chr(allele, \(x) ifelse(!is.na(x), ard$redux(x, "S"), x))
}

#' Scale down HLA-alleles to two-field resolution
#'
#' `reduce_to_field2()` is a thin wrapper around the `ard.redux()` method
#' (Reduction Type `'U2'`) of the
#' [`py-ard`](https://github.com/nmdp-bioinformatics/py-ard) Python package.
#'
#' @inherit reduce_to_serology details
#' @inheritParams reduce_to_serology
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the corresponding two-field alleles.
#' @export
#' @seealso [reduce_to_serology()] to reduce an allele to
#'  serological equivalents
#'
#' @examples
#' \dontrun{
#' ard <- db_initialize(data_dir = "~/ipd_db/")
#' reduce_to_field2(ard, "A*01:04:01:01N")
#' #> "A*01:04N"
#' reduce_to_field2(ard, "B*44:270:01")
#' #> "B*44:270"
#' # N.B. serology can also be "reduced", but might lead to a long genotype:
#' reduce_to_field2(ard, "Cw10")
#' #> "C*03:02/C*03:04Q/C*03:04/C*03:06/C*03:26/C*03:28/C*03:46"
#' # Also works for vectors:
#' reduce_to_field2(ard, c("B*44:270:01", "B*44:66"))
#' #> "B*44:270" "B*44:66"
#' }
reduce_to_field2 <- function(ard, allele) {
  purrr::map_chr(allele, \(x) ifelse(!is.na(x), ard$redux(x, "lgx"), x))
}

#' Encode an ambiguous HLA typing into a Multiple Allele Code (MAC)
#'
#' @description
#' `mac_lookup()` is a thin wrapper around the `ard.lookup_mac()` method of the
#' [`py-ard`](https://github.com/nmdp-bioinformatics/py-ard) Python package.
#'
#' See the [NMDP website](https://hml.nmdp.org/MacUI/) for more on MACs.
#'
#' @inherit reduce_to_serology details
#' @inheritParams reduce_to_serology
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the corresponding MACs.
#' @export
#' @seealso [mac_expand()] for the reverse operation
#'
#' @examples
#' \dontrun{
#' ard <- db_initialize(data_dir = "~/ipd_db/")
#' mac_lookup(ard, "A*01:01/A*01:02")
#' #> "A*01:AB"
#' mac_lookup(ard, "HLA-A*25:01/HLA-A*26:01")
#' #> "HLA-A*25:BYHR"
#' # Also works for vectors:
#' mac_lookup(ard, c("A*01:01/A*01:02", "A*25:01/A*26:01"))
#' #> "A*01:AB"   "A*25:BYHR"
#' }
mac_lookup <- function(ard, allele) {
  purrr::map_chr(allele, \(x) ifelse(!is.na(x), ard$lookup_mac(x), x))
}

#' Decode a Multiple Allele Code (MAC) into an ambiguous HLA typing
#'
#' @description
#' `mac_expand()` is a thin wrapper around the `ard.expand_mac()` method of the
#' [`py-ard`](https://github.com/nmdp-bioinformatics/py-ard) Python package.
#'
#' See the [NMDP website](https://hml.nmdp.org/MacUI/) for more on MACs.
#'
#' @inherit mac_lookup details
#' @inheritParams mac_lookup
#'
#' @return A string or character vector of the same length as `allele`,
#'   with the ambiguities.
#' @export
#' @seealso [mac_lookup()] for the reverse operation
#'
#' @examples
#' \dontrun{
#' ard <- db_initialize(data_dir = "~/ipd_db/")
#' mac_expand(ard, "A*01:AB")
#' #> "A*01:01/A*01:02"
#' mac_expand(ard, "HLA-A*25:BYHR")
#' #> "HLA-A*25:01/HLA-A*26:01"
#' # Also works for vectors:
#' mac_expand(ard, c("A*01:AB", "A*25:BYHR"))
#' #> "A*01:01/A*01:02" "A*25:01/A*26:01"
#' }
mac_expand <- function(ard, allele) {
  purrr::map_chr(allele, \(x) ifelse(!is.na(x), ard$expand_mac(x), x))
}

#' Lookup whether allele exists in IPD-IMGT/HLA database
#'
#' `is_in_ipd_db()` checks whether an HLA allele has been entered into the
#' [IPD-IMGT/HLA database](https://www.ebi.ac.uk/ipd/imgt/hla/) as
#' initialized with [db_initialize()]. This function is a thin wrapper around
#' the `ard.validate()` method of the
#' [`py-ard`](https://github.com/nmdp-bioinformatics/py-ard) Python package.
#'
#' @inheritParams reduce_to_serology
#'
#' @return A Boolean or logical vector with the same length as `allele`, with
#'  `TRUE` (allele exists in database) or `FALSE` for each element.
#' @export
#'
#' @examples
#' \dontrun{
#' ard <- db_initialize(data_dir = "~/ipd_db/")
#' is_in_ipd_db(ard, "A*01:01")
#' #> TRUE
#' is_in_ipd_db(ard, "A1")
#' #> TRUE
#' is_in_ipd_db(ard, "A20")
#' #> FALSE
#' is_in_ipd_db(ard, "A*99:99")
#' #> FALSE
#' # Also works for vectors:
#' is_in_ipd_db(ard, c("A*01:AB", "XYZ*01:03"))
#' #> TRUE FALSE
#' }
is_in_ipd_db <- function(ard, allele) {
  # ard$validate() throws an error if the allele is not in the db, instead of
  # returning "FALSE", so we need to wrap it
  safe_validate <- purrr::safely(ard$validate, otherwise = FALSE)
  purrr::map(allele, safe_validate) |>
    purrr::map("result") |>
    purrr::simplify()
}

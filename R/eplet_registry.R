get_positive_eplets <- function(luminex_df, sample_col, alleles_col,
                                assignment_col, eplet_df,
                                pos_col = "eplets_pos") {
  luminex_df |>
    dplyr::mutate(eplets = lookup_eplets(eplet_df, {{ alleles_col }})) |>
    tidyr::unnest_longer("eplets") |> # one row for each eplet
    dplyr::group_by({{ sample_col }}) |> # for each sample
    dplyr::reframe({{ pos_col }} := setdiff( # discard eplets that also occur on
      .data$eplets[{{ assignment_col }}], # negative beads
      .data$eplets[!{{ assignment_col }}]
    ))
}

#' Lookup eplets corresponding to alleles in HLA Eplet Registry
#'
#' `lookup_eplets()` takes in a set of HLA alleles, and retrieves the
#' corresponding eplets from the Eplet Registry table.
#'
#' @param eplet_df Data frame containing the Eplet Registry; from output of
#'   [load_eplet_registry()].
#' @param alleles String or character vector of HLA alleles.
#'
#' @return Named list, where each element is a character vector of eplets for
#'   each allele in the input.
#' @seealso [lookup_alleles()] for the reverse operation
#' @export
#'
#' @examples
#' \dontrun{
#' df_eplets <- load_eplet_registry()
#' lookup_eplets(df_eplets, "A*01:01")
#' # Also works for vectors:
#' lookup_eplets(df_eplets, c("A*01:01", "B*08:01"))
#' }
lookup_eplets <- function(eplet_df, alleles) {
  lookup_eplet <- function(eplet_df, allele) {
    eplet_df |>
      dplyr::filter(.data$alleles == allele) |>
      dplyr::distinct(.data$name) |>
      dplyr::pull(.data$name)
  }

  purrr::map(alleles, \(x) lookup_eplet(eplet_df, x)) |>
    purrr::set_names(alleles)
}

#' Lookup alleles corresponding to eplets in HLA Eplet Registry
#'
#' `lookup_alleles()` takes in a set of eplets, and retrieves the corresponding
#' HLA alleles from the Eplet Registry table.
#'
#' @param eplet_df Data frame containing the Eplet Registry; from output of
#'   [load_eplet_registry()].
#' @param eplets String or character vector of Eplet names.
#' @param allele_set Whether to return only Luminex alleles (`"luminex"`;
#'   default) or all alleles (`"all"`).
#'
#' @return Named list, where each each element is a character vector of alleles
#'   for each eplet in the input.
#' @seealso [lookup_eplets()] for the reverse operation
#' @export
#'
#' @examples
#' \dontrun{
#' df_eplets <- load_eplet_registry()
#' lookup_alleles(df_eplets, "9F")
#' lookup_alleles(df_eplets, "3P", allele_set = "all")
#' # Also works for vectors:
#' lookup_alleles(df_eplets, c("9F", "3S"))
#' }
lookup_alleles <- function(eplet_df, eplets, allele_set = "luminex") {
  rlang::arg_match(allele_set, c("luminex", "all"))

  eplet_df <- dplyr::filter(eplet_df, .data$source == allele_set)

  lookup_allele <- function(eplet_df, eplet) {
    eplet_df |>
      dplyr::filter(.data$name == eplet) |>
      dplyr::pull(.data$alleles)
  }

  purrr::map(eplets, \(x) lookup_allele(eplet_df, x)) |>
    purrr::set_names(eplets)
}

#' Load HLA Eplet Registry table
#'
#' `load_eplet_registry()` returns a dataframe with the [HLA Eplet Registry
#' table](https://www.epregistry.com.br), which maps HLA alleles to eplets. If
#' no local copy exists, the function will prompt you whether you want to scrape
#' the online registry to download a fresh copy.
#'
#' Currently incorporates only the HLA tables (A/B/C, DRB, DQ, DP and
#' Interlocus), not MICA.
#'
#' @param folder_path Character path to the folder where the previously
#'   downloaded table is stored, or where you want a new version to be stored.
#'   Defaults to the R user cache directory.
#' @param filename Character filename for the table.
#' @param print_version Logical indicating whether to print message with version
#'   information when loading the table. Turn this off by setting to `FALSE`.
#' @param return_path Logical indicating whether to return the folder path
#'   rather than the data frame.
#' @param delete Logical indicating whether to delete the cached Eplet Registry.
#' @return Either a data frame containing the Eplet Registry (default) or a
#'   character string giving the path to the folder where it is/will be stored
#'   (if `return_path` = `TRUE`).
#' @export
#' @seealso
#'   - [lookup_alleles()] for a helper function that takes in eplets and looks
#'   up the HLA alleles they occur on in the registry
#'   - [lookup_eplets()] for a helper function that takes in alleles and looks
#'   up which eplets occur on them in the registry
#' @examples
#' \dontrun{
#' # (Down)load the eplet registry (to/)from the default path
#' df_eplets <- load_eplet_registry()
#' # Return the location of the default path
#' load_eplet_registry(return_path = TRUE)
#' # (Down)load the eplet registry (to/)from a custom path and a custom filename
#' df_eplets <- load_eplet_registry("~/eplet_registry", "eplet_table.rds")
#' # Delete the saved table
#' load_eplet_registry("~/eplet_registry", "eplet_table.rds", delete = TRUE)
#' }
#'
load_eplet_registry <- function(folder_path = NULL,
                                filename = "eplets.rds",
                                print_version = TRUE,
                                return_path = FALSE,
                                delete = FALSE) {
  if (is.null(folder_path)) {
    folder_path <- tools::R_user_dir("hlapro", "cache")
  }

  if (return_path) {
    return(folder_path)
  }

  file_path <- file.path(folder_path, filename)

  if (delete) {
    unlink(file_path)
    return(invisible())
  }

  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }

  if (file.exists(file_path)) {
    df_eplet <- readRDS(file_path)
    if (print_version) {
      message(
        stringr::str_glue(
          "Loaded Eplet Registry table ({attr(df_eplet, 'notes')}),\n",
          "released {attr(df_eplet, 'date')}, ",
          "downloaded from {attr(df_eplet, 'url')}"
        )
      )
    }

    return(invisible(df_eplet))
  }

  if (rlang::is_interactive() && scrape_permission() == 2) {
    return(invisible())
  }

  invisible(scrape_eplet_registry(file_path))
}

fetch_registry_version <- function() {
  registry_url <- "https://www.epregistry.com.br"
  version_text <- rvest::read_html(registry_url) |>
    rvest::html_elements("body > footer > div.container > p:nth-child(2)") |>
    rvest::html_text2()

  version_date <- stringr::str_extract(version_text, r"(\d{4}-\d{2}-\d{2})")

  invisible(c(
    date = version_date,
    # extract text following version date
    notes = stringr::str_extract(version_text,
      stringr::str_glue(r"((?:{version_date}\.\s)(.*))"),
      group = 1
    ),
    url = registry_url
  ))
}

scrape_permission <- function() {
  q_title <- paste(
    "Do you want to download the HLA Eplet Registry tables",
    "(to lookup which eplets occur on which alleles, and vice versa?)"
  )

  utils::menu(choices = c("Yes", "No"), title = q_title)
}

scrape_eplet_registry <- function(file_path) {
  base_url <- "https://www.epregistry.com.br/index/databases/database/"
  databases <- c("ABC", "DRB", "DQ", "DP", "DRDQDP")
  # CSS selector paths to the individual columns
  # (scraping entire table with rvest::html_table resulted in misaligned rows/
  # columns)
  base_path <- "#table-result > div > table > tbody > tr > td:nth-child"
  col_paths <- c(
    id = "(1)",
    name = "(2)",
    description = "(3)",
    exposition = "(4)",
    confirmation = "(6)",
    alleles_luminex = "(9)",
    alleles_all = "(10) > div > div:nth-of-type(2) > p:nth-of-type(2)"
  )
  col_paths[] <- paste0(base_path, col_paths)

  scrape_column <- function(page_html, col_path) {
    rvest::html_elements(page_html, col_path) |>
      rvest::html_text2()
  }

  # scrape all columns, add each to a list, and store the database used
  scrape_table <- function(base_url, database, col_paths) {
    Sys.sleep(0.5) # wait a little between scrapes
    page_html <- rvest::read_html(paste0(base_url, database))
    purrr::map(col_paths, \(x) scrape_column(page_html, x)) |>
      purrr::list_assign(database = database)
  }

  # for each database (i.e. page), scrape all columns, store in another list
  df <- purrr::map(databases,
    \(x) scrape_table(base_url, x, col_paths),
    .progress = "Collecting tables from HLA Eplet Registry website"
  ) |>
    purrr::map(tidyr::as_tibble) |> # make a dataframe out of each scraped db
    purrr::list_rbind() |> # combine into one dataframe
    # clean up the column: text always starts with "Yes" if eplet confirmed
    dplyr::mutate(confirmation = stringr::str_starts(
      .data$confirmation, "Yes"
    )) |>
    # one column for the alleles, and another for if they're luminex or not
    tidyr::pivot_longer(c(.data$alleles_luminex, .data$alleles_all),
      names_to = "source",
      names_prefix = "alleles_",
      values_to = "alleles"
    ) |>
    dplyr::group_by(.data$id) |> # one row per allele
    tidyr::separate_longer_delim(.data$alleles, delim = ",") |>
    dplyr::filter(.data$alleles != "") |> # get rid of trailing comma artefact
    # clean up whitespace at start/end
    dplyr::ungroup() |>
    dplyr::mutate(dplyr::across(
      dplyr::where(is.character),
      ~ stringr::str_trim(.x)
    ))

  registry_info <- fetch_registry_version()
  attr(df, "date") <- registry_info[["date"]]
  attr(df, "notes") <- registry_info[["notes"]]
  attr(df, "url") <- registry_info[["url"]]

  saveRDS(df, file_path)
  df
}

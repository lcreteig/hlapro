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

scrape_eplet_registry <- function(base_url, databases) {
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
    page_html <- rvest::read_html(paste0(base_url, database))
    purrr::map(col_paths, \(x) scrape_column(page_html, x)) |>
      purrr::list_assign(database = database)
  }

  # for each database (i.e. page), scrape all columns, store in another list
  df <- purrr::map(databases, \(x) scrape_table(base_url, x, col_paths)) |>
    purrr::map(tidyr::as_tibble) |> # make a dataframe out of each scraped db
    purrr::list_rbind() |> # combine into one dataframe
    # clean up the column: text always starts with "Yes" if eplet confirmed
    dplyr::mutate(confirmation = stringr::str_starts(confirmation, "Yes")) |>
    # one column for the alleles, and another for if they're luminex or not
    tidyr::pivot_longer(c(alleles_luminex, alleles_all),
      names_to = "source",
      names_prefix = "alleles_",
      values_to = "alleles"
    ) |>
    dplyr::group_by(id) |>
    tidyr::separate_longer_delim(alleles, delim = ",") |> # one row per allele
    # clean up whitespace at start/end
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ stringr::str_trim(.x)))

  registry_info <- fetch_registry_version()
  attr(df, "date") <- registry_info[["date"]]
  attr(df, "notes") <- registry_info[["notes"]]
  attr(df, "url") <- registry_info[["url"]]

  df
}

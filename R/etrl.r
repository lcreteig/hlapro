load_etrl_tables <- function(print_version = FALSE,
                             return_path = FALSE,
                             delete = FALSE) {
  folder_path <- tools::R_user_dir("hlapro", "cache")
  file_name <- "etrl.rds"
  file_path <- file.path(folder_path, file_name)

  if (return_path) {
    return(folder_path)
  }

  if (delete) {
    unlink(file_path)
    return(invisible())
  }

  if (!dir.exists(folder_path)) {
    dir.create(folder_path)
  }

  if (file.exists(file_path)) {
    df_etrl <- readRDS(file_path)
    if (print_version) {
      message(
        stringr::str_glue(
          "Loaded v{attr(df_etrl, 'version')} of ETRL tables ",
          "(released {attr(df_etrl, 'date')}), ",
          "downloaded from {attr(df_etrl, 'url')} "
        )
      )
    }

    return(invisible(df_etrl))
  }

  if (dl_permission() == 2) {
    return(invisible())
  }

  invisible(download_etrl(file_path))
}

fetch_etrl_version <- function() {
  etrl_url <- "https://etrl.eurotransplant.org/resources/new-hla-tables/"
  tbls_page <- rvest::read_html(etrl_url)

  version_text <- tbls_page |>
    rvest::html_elements(".entry-content") |>
    rvest::html_elements("p") |>
    rvest::html_text2() |>
    stringr::str_flatten()

  invisible(c(
    version = stringr::str_extract(version_text, r"((?<=v)\d[\.\d]*\b)"),
    date = stringr::str_extract(version_text, r"(\d{1,2}-\d{1,2}-\d{4})"),
    url = etrl_url
  ))
}

download_etrl <- function(file_path) {
  base_url <- "https://etrl.eurotransplant.org/resources/"
  tbl_suffixes <- c(
    "a", "b", "c",
    "drb1", "drb3-4-5", "dqb1", "dqa1", "dpb1", "dpa1"
  )
  tbl_urls <- stringr::str_c(base_url, "hla-", tbl_suffixes, "/")

  scrape_etrl <- function(tbl_url) {
    Sys.sleep(0.1) # wait a little between scrapes
    rvest::read_html(tbl_url) |>
      rvest::html_element("table") |>
      rvest::html_table(na.strings = "")
  }

  df_etrl <- purrr::map(tbl_urls, scrape_etrl,
    .progress = "Collecting HLA tables from ETRL website"
  ) |>
    purrr::list_rbind()

  etrl_info <- fetch_etrl_version()
  attr(df_etrl, "version") <- etrl_info[["version"]]
  attr(df_etrl, "date") <- etrl_info[["date"]]
  attr(df_etrl, "url") <- etrl_info[["url"]]

  saveRDS(df_etrl, file_path)
  df_etrl
}

dl_permission <- function() {
  q_title <- paste(
    "Do you want to download the EuroTransplant Reference",
    "Laboratory (ETRL) HLA tables (to convert allele-level",
    "HLA typings to serological equivalents)?"
  )

  utils::menu(choices = c("Yes", "No"), title = q_title)
}

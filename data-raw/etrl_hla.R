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

download_etrl <- function() {
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

  df_etrl
}

q_title <- paste(
  "Do you want to download the EuroTransplant Reference",
  "Laboratory (ETRL) HLA tables (to convert allele-level",
  "HLA typings to serological equivalents)?"
)

if (dl_permission(q_title) == 1) {
  rlang::check_installed("rvest", reason = "to scrape the ETRL HLA tables")
  etrl_hla <- download_etrl()

  # Fix some errors in the table
  etrl_hla <- etrl_hla |>
    dplyr::mutate(
      Allele =
        dplyr::case_match(Allele,
          "A*24:xx" ~ "A*24:XX",
          "A32:XX" ~ "A*32:XX",
          "A34:XX" ~ "A*34:XX",
          "A36:XX" ~ "A*36:XX",
          "A43:XX" ~ "A*43:XX",
          "A66:XX" ~ "A*66:XX",
          "B37:XX" ~ "B*37:XX",
          "B51:XX" ~ "B*51:XX",
          "DRTB1*03:XX" ~ "DRB1*03:XX",
          "DQB1:03:09" ~ "DQB1*03:09",
          "DQB1:02:10" ~ "DQB1*02:10",
          "DQB1:04:XX" ~ "DQB1*04:XX",
          "DQB1:05:XX" ~ "DQB1*05:XX",
          "DQB1:06:XX" ~ "DQB1*06:XX",
          .default = Allele
        )
    )

  rlang::check_installed("usethis", reason = "to save the data")
  usethis::use_data(etrl_hla, overwrite = TRUE)
} else {
  stop("Download aborted")
}

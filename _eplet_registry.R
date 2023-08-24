library(tidyverse)
library(rvest)

base_url <- "https://www.epregistry.com.br/index/databases/database/"
databases <- c("ABC", "DRB", "DQ", "DP", "DRDQDP")

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
  html_elements(page_html, col_path) |>
    html_text2()
}

scrape_table <- function(base_url, database, col_paths) {
  page_html <- read_html(paste0(base_url, database))
  map(col_paths, \(x) scrape_column(page_html, x)) |>
    list_assign(database = database)
}

tbl <- map(databases, \(x) scrape_table(base_url, x, col_paths)) |>
  map(as_tibble) |>
  list_rbind() |>
  mutate(confirmation = str_starts(confirmation, "Yes")) |>
  pivot_longer(c(alleles_luminex, alleles_all),
    names_to = "source",
    names_prefix = "alleles_",
    values_to = "alleles"
  ) |>
  group_by(id) |>
  separate_longer_delim(alleles, delim = ",") |>
  mutate(across(everything(), ~ str_trim(.x)))

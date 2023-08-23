library(tidyverse)
library(rvest)

page_url <- "https://www.epregistry.com.br/index/databases/database/ABC/"
tbls_page <- read_html(page_url)

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


df <- map(col_paths, \(x) html_text2(html_elements(tbls_page, x))) |>
  as_tibble() |>
  mutate(confirmation = str_starts(confirmation, "Yes")) |>
  pivot_longer(c(alleles_luminex, alleles_all),
    names_to = "source",
    names_prefix = "alleles_",
    values_to = "alleles"
  ) |>
  group_by(id) |>
  separate_longer_delim(alleles, delim = ",") |>
  mutate(across(everything(), ~ str_trim(.x)))

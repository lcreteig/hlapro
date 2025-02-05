q_title <- paste(
  "Do you want to download a file from GitHub that contains",
  "deleted and changed HLA alleles?"
)

if (dl_permission(q_title) == 1) {
  rlang::check_installed("readr", reason = "to read the data")

  url_ipd <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/"

  deleted_changed <- readr::read_csv(paste0(url_ipd, "Deleted_alleles.txt"),
    skip = 7
  ) |>
    dplyr::filter(stringr::str_starts(AlleleID, "HLA")) |>
    dplyr::rename(allele_old = Allele) |>
    # extract allele name from description (first 2+field allele)
    dplyr::mutate(allele_new = stringr::str_extract(
      Description, r"(\w+\*\d{2,}:[\d:]+[NLSCAQ]?)"
    )) |>
    dplyr::mutate(date_changed = as.Date(
      paste("1", stringr::str_extract( # assume first of the month
        Description, r"(\w+\s\d+(?=\)))" # date is "Month Year)"
      )),
      format = "%d %B %Y"
    )) |>
    dplyr::select(allele_old, allele_new, date_changed)

  rlang::check_installed("usethis", reason = "to save the data")
  usethis::use_data(deleted_changed, overwrite = TRUE) # save dataframe
} else {
  stop("Download aborted")
}

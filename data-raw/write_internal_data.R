load(file = "data/etrl_hla.rda")
load(file = "data/v2_to_v3.rda")
load(file = "data/deleted_changed.rda")

make_public_lookup <- function(etrl_hla) {
  etrl_hla |>
    tidyr::pivot_longer(
      c(
        `ET MatchDeterminantSplit`,
        `ET MatchDeterminantBroad`,
      ),
      names_to = NULL,
      values_to = "broad_split"
    ) |>
    # Filter out the broads/splits that map to more than one public epitope
    dplyr::distinct(broad_split, Public) |>
    dplyr::add_count(broad_split) |>
    dplyr::filter(n == 1 & !(is.na(broad_split) | is.na(Public))) |>
    dplyr::select(-n) |>
    tibble::deframe()
}

make_broad_split_lookup <- function(etrl_hla) {
  etrl_hla |>
    dplyr::filter(!is.na(`ET MatchDeterminantSplit`)) |>
    dplyr::distinct(`ET MatchDeterminantSplit`, `ET MatchDeterminantBroad`) |>
    tibble::deframe()
}

etrl_split_to_broad <- make_broad_split_lookup(etrl_hla)
etrl_public <- make_public_lookup(etrl_hla)

lookup_v3 <- tibble::deframe(v2_to_v3)
lookup_del_chg <- tibble::deframe(dplyr::select(deleted_changed, -date_changed))

usethis::use_data(etrl_hla, etrl_split_to_broad, etrl_public,
  lookup_v3, lookup_del_chg,
  overwrite = TRUE, internal = TRUE
)

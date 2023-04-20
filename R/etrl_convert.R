etrl_lookup <- function(allele) {
  df_etrl <- load_etrl_tables()
  ids <- match(etrl_convert(allele), df_etrl$Allele)
  df_etrl[ids, ]
}

etrl_convert <- function(allele) {
  allele <- remove_hla_prefix(allele)

  df_etrl <- load_etrl_tables()
  allele_f2 <- reduce_to_nth_field(allele, 2)

  # If allele in ETRL table, us as is, otherwise make into xx code
  ifelse(allele_f2 %in% df_etrl$Allele, allele_f2, make_xx(allele_f2)) |>
    replace(has_suffix(allele), "") |> # suffixes cannot be reduced
    ifelse(is_serology(allele), allele, no = _) # return serology as is
}

make_xx <- function(allele) {
  stringr::str_c(reduce_to_nth_field(allele, 1), ":XX")
}

is_serology <- function(allele) {
  !stringr::str_detect(allele, r"([\*:])")
}

has_suffix <- function(allele) {
  stringr::str_detect(allele, "[NLSCAQ]$")
}

remove_suffixes_groups <- function(allele) {
  stringr::str_remove(allele, "[NLSCAQPG]$")
}

remove_hla_prefix <- function(allele) {
  stringr::str_remove(allele, "^HLA-")
}

is_mac <- function(allele) {
  stringr::str_detect(allele, r"(\*\d+:(?!XX)[A-Z])")
}

is_ambiguous <- function(allele) {
  stringr::str_detect(allele, r"(\*\d+:\d+\/)")
}

reduce_to_nth_field <- function(allele, n) {
  allele <- remove_suffixes_groups(allele)
  # logical index of all alleles to be reduced
  res_idx <- get_n_fields(allele) > n

  if (!any(res_idx)) {
    return(allele)
  }

  # Get list of field separator locations in each allele string
  field_locs <- stringr::str_locate_all(allele[res_idx], ":")
  ends <- purrr::map_int(field_locs, n) - 1 # end of nth field

  # Keep only the part of the typing up until the nth ":"
  replace(allele, res_idx, stringr::str_sub(allele[res_idx], 1, ends))
}

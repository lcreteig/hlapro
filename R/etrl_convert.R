etrl_convert <- function(allele) {
  allele <- remove_group_suffix(remove_hla_prefix(allele))

  if (is_serology(allele)) {
    return(allele)
  }

  if (has_suffix(allele)) {
    return("")
  }

  allele_to_f2 <- reduce_to_nth_field(allele, 2)
  df_etrl <- load_etrl_tables()
  if (allele_to_f2 %in% df_etrl$Allele) { # if in ETRL table
    allele_to_f2
  } else {
    make_xx(allele)
  }
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

remove_group_suffix <- function(allele) {
  stringr::str_remove(allele, "[PG]$")
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
  if (get_n_fields(allele) <= n) {
    return(allele)
  }
  field_locs <- stringr::str_locate_all(allele, ":")[[1]]
  stringr::str_sub(allele, 1, field_locs[n] - 1)
}

upscale_typing <- function(filepath,
                           typing,
                           loci = c("A", "B", "DRB1", "DRB.", "DQB1"),
                           population = "EURCAU",
                           n_haplos = NULL,
                           n_genos = 1) {
  typing_alleles <- extract_alleles_str(typing,
    strip_locus = FALSE, loci = loci
  )
  typing_alleles <- typing_alleles[!is.na(typing_alleles)]

  translate_top_haplos(filepath,
    loci = loci,
    population = population,
    n_haplos = n_haplos
  ) |>
    select_compatible_haplos(typing_alleles) |>
    make_phased_genotypes(typing_alleles) |>
    make_unphased_genotypes() |>
    dplyr::slice_max(unphased_prob, n = n_genos)
}

translate_top_haplos <- function(filepath, loci, population, n_haplos) {
  haplo_df <- readxl::read_xlsx(filepath, na = "NA") |>
    dplyr::rename(`DRB.` = `DRB3-4-5`) |> # `DRB.` is what rest of package used
    # make generic name for the "freq" and "rank" cols, for ease of use later on
    dplyr::rename_with(\(x) stringr::str_replace_all(x, population, "haplo")) |>
    dplyr::select(c(A, B, C, `DRB.`, DRB1, DQB1, haplo_freq, haplo_rank)) |>
    dplyr::mutate(dplyr::across(all_of(loci), ~ stringr::str_remove(.x, "g")))

  if (is.null(n_haplos)) {
    n_haplos <- max(haplo_df$haplo_rank, na.rm = TRUE)
  }

  haplo_df |>
    dplyr::slice_min(haplo_rank, n = n_haplos) |>
    # translate alleles to serological equivalents
    dplyr::mutate(
      dplyr::across(all_of(loci), ~ get_broad(.x), .names = "{.col}_broad"),
      dplyr::across(all_of(loci), ~ get_split(.x), .names = "{.col}_split")
    )
}

select_compatible_haplos <- function(df, typing_alleles) {
  loci_present <- unique(stringr::str_extract(
    names(typing_alleles),
    r"(^.*(?=_))"
  ))

  df |>
    tidyr::pivot_longer(
      cols = tidyr::ends_with(c("_broad", "_split")),
      names_to = c("locus", ".value"),
      names_sep = "_"
    ) |>
    dplyr::filter(locus %in% loci_present) |>
    dplyr::group_by(haplo_rank) |>
    dplyr::filter(all(broad %in% typing_alleles | split %in% typing_alleles)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(
      names_from = locus,
      values_from = c(broad, split),
      names_glue = "{locus}_{.value}"
    )
}

make_phased_genotypes <- function(df, typing_alleles) {
  df |>
    dplyr::cross_join(df, suffix = c("_1", "_2")) |> # combine all permutations
    # make heterozygous phased genotypes: keep only unique haplo combinations
    dplyr::filter(haplo_rank_1 < haplo_rank_2) |>
    tibble::rowid_to_column(var = "id_phased_geno") |>
    # make genotype out of serological equivalents  in phased genotypes
    tidyr::pivot_longer(
      cols = tidyr::contains(c("_broad_", "_split_")),
      names_to = "locus_res_allele",
      values_to = "ser_typing"
    ) |>
    dplyr::group_by(id_phased_geno) |>
    dplyr::filter(all(typing_alleles %in% ser_typing)) |>
    dplyr::ungroup() |>
    dplyr::select(!c(locus_res_allele, ser_typing)) |>
    dplyr::distinct(id_phased_geno, .keep_all = TRUE) |>
    dplyr::mutate(
      phased_freq = 2 * haplo_freq_1 * haplo_freq_2,
      phased_prob = phased_freq / sum(phased_freq)
    )
}

make_unphased_genotypes <- function(df) {
  df |>
    tidyr::pivot_longer(
      cols = tidyr::ends_with(c("_1", "_2")),
      names_to = c(".value", "haplo"),
      names_pattern = "(.*)_(\\d)"
    ) |>
    tidyr::pivot_longer(
      cols = c(A, B, C, `DRB.`, DRB1, DQB1),
      names_to = "locus",
      values_to = "typing"
    ) |>
    dplyr::group_by(id_phased_geno) |>
    dplyr::mutate(unphased_geno = stringr::str_flatten(
      stringr::str_sort(typing), " "
    )) |>
    dplyr::select(!c(locus, typing)) |>
    dplyr::distinct(id_phased_geno, haplo_rank, .keep_all = TRUE) |>
    tidyr::pivot_wider(
      names_from = haplo,
      values_from = c(haplo_freq, haplo_rank)
    ) |>
    dplyr::group_by(unphased_geno) |>
    dplyr::mutate(
      id_unphased_geno = dplyr::cur_group_id(), .before = dplyr::everything(),
      unphased_freq = sum(phased_freq),
      unphased_prob = sum(phased_prob)
    ) |>
    dplyr::ungroup() |>
    dplyr::relocate(unphased_geno, .after = id_unphased_geno) |>
    dplyr::select(-id_phased_geno)
}

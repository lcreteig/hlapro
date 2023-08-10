upscale_typings <- function(filepath,
                            typings,
                            loci = c("A", "B", "DRB1", "DRB.", "DQB1"),
                            population = "EURCAU",
                            n_haplos = NULL,
                            n_genos = 1,
                            as_list = FALSE) {
  # Load data and select haplotypes for given population
  haplo_df <- translate_top_haplos(filepath,
    loci = loci,
    population = population,
    n_haplos = n_haplos
  )

  upscale_typing <- function(haplo_df, typing, loci) {
    alleles <- extract_alleles_str(typing, strip_locus = FALSE, loci = loci)
    alleles <- alleles[!is.na(alleles)]

    haplo_df |>
      select_compatible_haplos(alleles) |>
      make_phased_genotypes(alleles) |>
      make_unphased_genotypes() |>
      dplyr::slice_max(.data$unphased_prob, n = n_genos)
  }

  res <- purrr::map(typings, \(x) upscale_typing(haplo_df, x, loci))
  if (as_list) { return(res) }
  purrr::list_rbind(res, names_to = "id_input_typing")
}

translate_top_haplos <- function(filepath, loci, population, n_haplos) {
  haplo_df <- readxl::read_xlsx(filepath, na = "NA") |>
    dplyr::rename(`DRB.` = .data$`DRB3-4-5`) |> # `DRB.` cf. rest of package
    # make generic name for the "freq" and "rank" cols, for ease of use later on
    dplyr::rename_with(\(x) stringr::str_replace_all(x, population, "haplo")) |>
    dplyr::select(dplyr::all_of(c(
      "A", "B", "C",
      "DRB.", "DRB1", "DQB1",
      "haplo_freq", "haplo_rank"
    ))) |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(loci), ~ stringr::str_remove(.x, "g")
    ))

  if (is.null(n_haplos)) {
    n_haplos <- max(haplo_df$haplo_rank, na.rm = TRUE)
  }

  haplo_df |>
    dplyr::slice_min(.data$haplo_rank, n = n_haplos) |>
    # translate alleles to serological equivalents
    dplyr::mutate(
      dplyr::across(dplyr::all_of(loci),
        ~ get_broad(.x),
        .names = "{.col}_broad"
      ),
      dplyr::across(dplyr::all_of(loci),
        ~ get_split(.x),
        .names = "{.col}_split"
      )
    )
}

select_compatible_haplos <- function(df, alleles) {
  loci_present <- unique(stringr::str_extract(
    names(alleles),
    r"(^.*(?=_))"
  ))

  df |>
    tidyr::pivot_longer(
      cols = tidyr::ends_with(c("_broad", "_split")),
      names_to = c("locus", ".value"),
      names_sep = "_"
    ) |>
    dplyr::filter(.data$locus %in% loci_present) |>
    dplyr::group_by(.data$haplo_rank) |>
    dplyr::filter(all(.data$broad %in% alleles | .data$split %in% alleles)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(
      names_from = .data$locus,
      values_from = c(.data$broad, .data$split),
      names_glue = "{locus}_{.value}"
    )
}

make_phased_genotypes <- function(df, alleles) {
  df |>
    dplyr::cross_join(df, suffix = c("_1", "_2")) |> # combine all permutations
    # make heterozygous phased genotypes: keep only unique haplo combinations
    dplyr::filter(.data$haplo_rank_1 < .data$haplo_rank_2) |>
    tibble::rowid_to_column(var = "id_phased_geno") |>
    # make genotype out of serological equivalents  in phased genotypes
    tidyr::pivot_longer(
      cols = tidyr::contains(c("_broad_", "_split_")),
      names_to = "locus_res_allele",
      values_to = "ser_typing"
    ) |>
    dplyr::group_by(.data$id_phased_geno) |>
    dplyr::filter(all(alleles %in% .data$ser_typing)) |>
    dplyr::ungroup() |>
    dplyr::select(!c(.data$locus_res_allele, .data$ser_typing)) |>
    dplyr::distinct(.data$id_phased_geno, .keep_all = TRUE) |>
    dplyr::mutate(
      phased_freq = 2 * .data$haplo_freq_1 * .data$haplo_freq_2,
      phased_prob = .data$phased_freq / sum(.data$phased_freq)
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
      cols = dplyr::all_of(c("A", "B", "C", "DRB.", "DRB1", "DQB1")),
      names_to = "locus",
      values_to = "typing"
    ) |>
    dplyr::group_by(.data$id_phased_geno) |>
    dplyr::mutate(unphased_geno = stringr::str_flatten(
      stringr::str_sort(.data$typing), " "
    )) |>
    dplyr::select(!c(.data$locus, .data$typing)) |>
    dplyr::distinct(.data$id_phased_geno, .data$haplo_rank, .keep_all = TRUE) |>
    tidyr::pivot_wider(
      names_from = .data$haplo,
      values_from = c(.data$haplo_freq, .data$haplo_rank)
    ) |>
    dplyr::group_by(.data$unphased_geno) |>
    dplyr::mutate(
      id_unphased_geno = dplyr::cur_group_id(), .before = dplyr::everything(),
      unphased_freq = sum(.data$phased_freq),
      unphased_prob = sum(.data$phased_prob)
    ) |>
    dplyr::ungroup() |>
    dplyr::relocate(.data$unphased_geno, .after = .data$id_unphased_geno) |>
    dplyr::select(-.data$id_phased_geno)
}

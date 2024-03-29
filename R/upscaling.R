#' Impute the resolution of an HLA typing from low to high (second field)
#'
#' `upscale_typings()` takes in a (vector of) low-resolution HLA typing
#' string(s), and "upscales" it to(two-field) high resolution based on
#' [haplotype frequencies from the NMDP
#' registry](http://dx.doi.org/10.1016/j.humimm.2013.06.025).
#'
#' This function uses haplotype frequencies published by the NMDP at
#' <https://frequency.nmdp.org>. You'll need to login and accept the license to
#' obtain the data, which is why we cannot distribute it with this package.
#'
#' ## Imputation algorithm
#'
#' Roughly, the function performs the following steps:
#'
#'  1. Translate all haplotype alleles from high-resolution to their serological
#' equivalents
#'  2. Select compatible haplotypes, i.e. those haplotypes with alleles that
#'  _all_ occur in the input genotype
#'  3. Combine all compatible haplotypes into phased genotypes (i.e. unique
#'  haplotype pairs). Fully homozygous haplotypes are not considered
#'  4. Calculate the frequency and probability of the phased genotypes
#'  5. Combine phased genotypes into unique unphased genotypes, and calculate
#'  their frequency and probability
#'
#' For a detailed explanation of the terminology and imputation algorithm,
#' see the following references:
#'
#' - Geffard et. al., Easy-HLA: a validated web application suite to reveal the
#'   full details of HLA typing, _Bioinformatics, Volume 36, Issue 7_, April
#'   2020, Pages 2157-2164, <https://doi.org/10.1093/bioinformatics/btz875>
#' - Madbouly, A., Gragert, L., Freeman, J., Leahy, N., Gourraud, P.-.-A.,
#'   Hollenbach, J.A., Kamoun, M., Fernandez-Vina, M. and Maiers, M. (2014),
#'   Validation of statistical imputation of allele-level multilocus phased
#'   genotypes from ambiguous HLA assignments. _Tissue Antigens, 84_: 285-292.
#'   <https://doi.org/10.1111/tan.12390>
#'
#' @param filepath String with the path to the `HLA-A~C~B~DRB1~DQB1.xlsx` file
#'   as downloaded from <https://frequency.nmdp.org>.
#' @param typings String or character vector of HLA typings, with
#'   space-separated alleles in serological notation.
#' @param loci_input Character vector of loci in the input genotype to be used
#'   for upscaling. Must be a subset of `c("A", "B", "C", "DRB1", "DRB.",
#'   "DQB1")`, where `DRB.` is DRB3/4/5. If one of these loci does not occur in
#'   the input typing, it will be ignored.
#' @param loci_output Character vector of loci that the upscaled, output
#'   genotype should contain. Can be different from `loci_input`, for example if
#'   there's no typing at all for a certain locus and you want to infer it from
#'   the haplotypes. Must be a subset of `c("A", "B", "C", "DRB1", "DRB.",
#'   "DQB1")`, where `DRB.` is DRB3/4/5.
#' @param population String specifying which population group to use the
#'   haplotype frequencies from. Must correspond to one of the abbreviations
#'   that can be found on <https://haplostats.org>
#'   Defaults to `"EURCAU"` (the largest population in the dataset).
#' @param n_haplos Number of most frequent haplotypes to use for the upscaling,
#'   e.g. `5000` to consider only the 5000 haplotypes with the highest frequency
#'   (the rest is discarded). Defaults to all haplotypes with a non-zero
#'   frequency in the selected population.
#' @param n_genos Number of output genotypes to return for each input genotype,
#'   sorted by probability (frequency) of the output genotypes. Defaults to the
#'   most likely genotype only (i.e., `1`).
#' @param as_list Boolean (`TRUE` or `FALSE`) that determines whether to return
#'   the result as a single dataframe, or a list of dataframes: one for each
#'   input typing (the latter can be useful when the input typings also live in
#'   a data frame; see examples).
#'
#' @return A (list of) dataframe(s) with the upscaling results, containing the
#' following columns:
#'    1. `id_input_typing`: Sequential identifier for the input typing
#'    2. `id_unphased_geno`: Unique identifier for each output unphased genotype
#'    3. `unphased_geno`: Upscaled genotype
#'    4. `unphased_freq`: Frequency of upscaled genotype
#'    5. `unphased_prob`: Probability of upscaled genotype
#'    6. `phased_freq`: Frequency of phased genotype (= pair of haplotypes) that
#'        make up the unphased genotype (can be many-to-one)
#'    7. `phased_prob`: Probability of the phased genotype
#'    8. `haplo_freq_1`: Frequency of the first haplotype in the phased genotype
#'    9. `haplo_rank_1`: Rank (descending) of the frequency of this haplotype
#'    8. `haplo_freq_2`: Frequency of the 2nd haplotype in the phased genotype
#'    9. `haplo_rank_2`: Rank (descending) of the frequency of this haplotype
#' @export
#' @examples
#' \dontrun{
#' upscale_typings(
#'   filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
#'   typing = "A24 A28 B35 B61 DR4 DR11"
#' )
#'
#' # If your GL Strings are in a data frame with some ID'ing columns that you
# want to keep attached, call `gl_to_df()` on the GL String column in your
# data frame:
#' # If you've got more than one typing to upscale, perhaps along with some
#' # ID'ing columns (e.g. patient ID), it's probably best to put them in a
#' # data frame and call `upscale_typings()` on your data frame:
#' library(tidyverse)
#'
#' typing_df <- tibble(
#'   id = c("001", "002"),
#'   input_typings = c(
#'     "A24 A28 B35 B61 DR4 DR11",
#'     "A2 A3 B52 B35 Cw4 DR11 DR52 DQ3"
#'   )
#' )
#' typing_df |>
#'   mutate(geno_df = upscale_typings(
#'     "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
#'     input_typings,
#'     as_list = TRUE
#'   )) |>
#'   unnest(geno_df)
#' }
upscale_typings <- function(filepath,
                            typings,
                            loci_input = c(
                              "A", "B",
                              "DRB1", "DRB.", "DQB1"
                            ),
                            loci_output = c(
                              "A", "B", "C",
                              "DRB1", "DRB.", "DQB1"
                            ),
                            population = "EURCAU",
                            n_haplos = NULL,
                            n_genos = 1,
                            as_list = FALSE) {
  loci_input <- rlang::arg_match(loci_input, c(
    "A", "B", "C",
    "DRB1", "DRB.", "DQB1"
  ), multiple = TRUE)
  loci_output <- rlang::arg_match(loci_output, multiple = TRUE)
  check_number_whole(n_genos)

  # Load data and select haplotypes for given population
  haplo_df <- translate_top_haplos(filepath,
    loci_input = loci_input,
    loci_output = loci_output,
    population = population,
    n_haplos = n_haplos
  )

  upscale_typing <- function(haplo_df, typing, loci_input) {
    # get alleles in input genotyping into a vector; only loci we want
    # consider for upscaling
    alleles <- extract_alleles_str(typing,
      strip_locus = FALSE,
      loci = loci_input
    )
    alleles <- alleles[!is.na(alleles)]

    haplo_df |>
      select_compatible_haplos(alleles) |>
      make_phased_genotypes(alleles) |>
      make_unphased_genotypes(loci_output) |>
      dplyr::slice_max(.data$unphased_prob, n = n_genos) # select top n_genos
  }

  # Perform upscaling for each typing
  res <- purrr::map(typings, \(x) upscale_typing(haplo_df, x, loci_input))
  if (as_list) {
    return(res)
  }
  purrr::list_rbind(res, names_to = "id_input_typing")
}

translate_top_haplos <- function(filepath, loci_input, loci_output,
                                 population, n_haplos) {
  haplo_df <- readxl::read_xlsx(filepath, na = "NA")

  stopifnot(
    "`population` must occur in NMDP dataset" =
      any(stringr::str_detect(colnames(haplo_df), population))
  )

  haplo_df <- haplo_df |>
    dplyr::rename(`DRB.` = dplyr::matches("DRB3-4-5")) |> # `DRB3-4-5` to `DRB.`
    # make generic name for the "freq" and "rank" cols, for ease of use later on
    dplyr::rename_with(\(x) stringr::str_replace_all(x, population, "haplo")) |>
    dplyr::select(dplyr::all_of(c(loci_output, "haplo_freq", "haplo_rank")))

  # determine how many haplos to include if not user specified
  if (is.null(n_haplos)) {
    n_haplos <- max(haplo_df$haplo_rank, na.rm = TRUE)
  }
  check_number_whole(n_haplos)

  haplo_df |>
    dplyr::slice_min(.data$haplo_rank, n = n_haplos) |>
    # translate alleles to serological equivalents
    dplyr::mutate(
      dplyr::across(dplyr::all_of(loci_input), # new columns for broads
        ~ get_broad(stringr::str_remove(.x, "g")),
        .names = "{.col}_broad"
      ),
      dplyr::across(dplyr::all_of(loci_input), # new columns for splits
        ~ get_split(stringr::str_remove(.x, "g")),
        .names = "{.col}_split"
      )
    )
}

select_compatible_haplos <- function(df, alleles) {
  loci_present <- unique(stringr::str_extract(names(alleles), r"(^.*(?=_))"))

  df |>
    tidyr::pivot_longer(
      cols = tidyr::ends_with(c("_broad", "_split")),
      names_to = c("locus", ".value"),
      names_sep = "_"
    ) |>
    dplyr::filter(.data$locus %in% loci_present) |>
    dplyr::group_by(.data$haplo_rank) |>
    # keep only downscaled haplos where *all* alleles occur in input genotype,
    # to prevent having to evaluate all haplos
    dplyr::filter(all(.data$broad %in% alleles | .data$split %in% alleles)) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider( # one row per haplotype pair ()
      names_from = "locus",
      values_from = tidyr::all_of(c("broad", "split")),
      names_glue = "{locus}_{.value}"
    )
}

make_phased_genotypes <- function(df, alleles) {
  df |>
    dplyr::cross_join(df, suffix = c("_1", "_2")) |> # combine all permutations
    # make phased genotypes: keep only unique haplo combinations
    dplyr::filter(.data$haplo_rank_1 <= .data$haplo_rank_2) |>
    tibble::rowid_to_column(var = "id_phased_geno") |>
    # make genotype out of serological equivalents  in phased genotypes
    tidyr::pivot_longer(
      cols = tidyr::contains(c("_broad_", "_split_")),
      names_to = "locus_res_allele",
      values_to = "ser_typing"
    ) |>
    dplyr::group_by(.data$id_phased_geno) |>
    # keep only phased genotypes where *all* alleles occur in input genotype
    dplyr::filter(all(alleles %in% .data$ser_typing)) |>
    dplyr::ungroup() |>
    dplyr::select(!dplyr::all_of(c("locus_res_allele", "ser_typing"))) |>
    dplyr::distinct(.data$id_phased_geno, .keep_all = TRUE) |>
    dplyr::mutate( # calculate phased genotype frequency / probability
      # for phased frequency calculation, change multiplier if homo/hetero
      multiplier = ifelse(.data$haplo_rank_1 == .data$haplo_rank_2, 1, 2),
      phased_freq = .data$multiplier * .data$haplo_freq_1 * .data$haplo_freq_2,
      phased_prob = .data$phased_freq / sum(.data$phased_freq)
    )
}

make_unphased_genotypes <- function(df, loci_output) {
  df |>
    tidyr::pivot_longer(
      cols = tidyr::ends_with(c("_1", "_2")),
      names_to = c(".value", "haplo"),
      names_pattern = "(.*)_(\\d)"
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(loci_output),
      names_to = "locus",
      values_to = "typing"
    ) |>
    dplyr::group_by(.data$id_phased_geno) |>
    # for each phased genotype: sort alleles and combine into unphased genotype
    dplyr::mutate(unphased_geno = stringr::str_flatten(
      stringr::str_sort(.data$typing), " "
    )) |>
    dplyr::select(!dplyr::all_of(c("locus", "typing"))) |>
    dplyr::distinct(.data$id_phased_geno, .data$haplo_rank, .keep_all = TRUE) |>
    tidyr::pivot_wider( # one row with both haplotypes
      names_from = "haplo",
      values_from = tidyr::all_of(c("haplo_freq", "haplo_rank"))
    ) |>
    dplyr::group_by(.data$unphased_geno) |> # for each distinct unphased geno
    dplyr::mutate( # calculate frequency and probability
      id_unphased_geno = dplyr::cur_group_id(), .before = dplyr::everything(),
      unphased_freq = sum(.data$phased_freq),
      unphased_prob = sum(.data$phased_prob)
    ) |>
    dplyr::ungroup() |>
    dplyr::relocate("unphased_geno", .after = "id_unphased_geno") |>
    dplyr::select(!("id_phased_geno"))
}

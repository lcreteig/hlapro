library(tidyverse)
library(readxl)
library(hlapro)

haplo_df <- read_xlsx("~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx")

top_haplos_ser <- haplo_df |>
  rename(`DRB.` = `DRB3-4-5`) |>
  mutate(across(c(A, B, `DRB.`, DRB1, DQB1), ~ str_remove(.x, "g"))) |> # remove "g" suffix
  slice_max(EURCAU_freq, n = 5000) |> # cf. BW, pick top 5000 haplos
  select(A, B, C, `DRB.`, DRB1, DQB1, EURCAU_freq, EURCAU_rank) |>
  # translate alleles to serological equivalents
  mutate(
    across(c(A, B, `DRB.`, DRB1, DQB1), ~ get_broad(.x), .names = "{.col}_broad"),
    across(c(A, B, `DRB.`, DRB1, DQB1), ~ get_split(.x), .names = "{.col}_split")
  )

typing <- "A2 A3 B52 B35 Cw4 DR11 DR52 DQ3"
typing <- "A24 A28 B35 B61 DR4 DR11"

typing_alleles <- extract_alleles_str(typing,
  strip_locus = FALSE,
  loci = c("A", "B", "DRB1", "DRB.", "DQB1")
)
typing_alleles <- typing_alleles[!is.na(typing_alleles)]

loci_present <- unique(str_extract(names(typing_alleles), r"(^.*(?=_))"))

# keep only haplos with alleles compatible with typing
comp_haplos <- top_haplos_ser |>
  pivot_longer(
    cols = ends_with(c("_broad", "_split")),
    names_to = c("locus", ".value"),
    names_sep = "_"
  ) |>
  filter(locus %in% loci_present) |>
  group_by(EURCAU_rank) |>
  filter(all(broad %in% typing_alleles | split %in% typing_alleles)) |>
  ungroup() |>
  pivot_wider(
    names_from = locus,
    values_from = c(broad, split),
    names_glue = "{locus}_{.value}"
  )

phased_genotypes <- comp_haplos |>
  cross_join(comp_haplos, suffix = c("_1", "_2")) |> # combine all permutations of haplos
  filter(EURCAU_rank_1 < EURCAU_rank_2) |> # make heterozygous phased genotypes: keep only unique combinations of haplos
  rowid_to_column(var = "id_phased_geno") |>
  # make genotype out of serological equivalents of alleles in phased genotypes
  pivot_longer(
    cols = contains(c("_broad_", "_split_")),
    names_to = "locus_res_allele",
    values_to = "ser_typing"
  ) |>
  group_by(id_phased_geno) |>
  filter(all(typing_alleles %in% ser_typing)) |>
  ungroup() |>
  select(!c(locus_res_allele, ser_typing)) |>
  distinct(id_phased_geno, .keep_all = TRUE) |>
  mutate(
    freq_phased = 2 * EURCAU_freq_1 * EURCAU_freq_2,
    prob_phased = freq_phased / sum(freq_phased)
  )

unphased_genotypes <- phased_genotypes |>
  pivot_longer(
    cols = ends_with(c("_1", "_2")),
    names_to = c(".value", "haplo"),
    names_pattern = "(.*)_(\\d)"
  ) |>
  pivot_longer(
    cols = c(A, B, C, `DRB.`, DRB1, DQB1),
    names_to = "locus",
    values_to = "typing"
  ) |>
  group_by(id_phased_geno) |>
  mutate(unphased_geno = str_flatten(str_sort(typing), " ")) |>
  select(!c(locus, typing)) |>
  distinct(id_phased_geno, EURCAU_rank, .keep_all = TRUE) |>
  pivot_wider(
    names_from = haplo,
    names_prefix = "haplo_",
    values_from = c(EURCAU_freq, EURCAU_rank)
  ) |>
  group_by(unphased_geno) |>
  mutate(
    id_unphased_geno = cur_group_rows(), .before = everything(),
    freq_unphased = sum(freq_phased),
    prob_unphased = sum(prob_phased)
  ) |>
  ungroup()

unphased_genotypes |>
  relocate(unphased_geno, .after = id_unphased_geno) |>
  slice_max(prob_unphased)

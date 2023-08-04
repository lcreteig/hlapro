library(tidyverse)
library(readxl)
library(hlapro)

haplo_df <- read_xlsx("~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx")

top_haplos_ser <- haplo_df |>
  mutate(across(c(A, B, C, `DRB3-4-5`, DRB1, DQB1), ~ str_remove(.x, "g"))) |> # remove "g" suffix
  slice_max(EURCAU_freq, n = 5000) |> # cf. BW, pick top 5000 haplos
  select(A, B, C, `DRB3-4-5`, DRB1, DQB1, EURCAU_freq, EURCAU_rank) |>
  # translate alleles to serological equivalents
  mutate(across(c(A, B, DRB1), ~ get_broad(.x), .names = "{.col}_broad"),
         across(c(A, B, DRB1), ~ get_split(.x), .names = "{.col}_split"))

typing <- "A24 A28 B35 B61 DR4 DR11"

typing_alleles <- extract_alleles_str(typing,
                                      strip_locus = FALSE,
                                      loci = c("A", "B", "DRB1"))

comp_haplos <- top_haplos_ser |> # keep only haplos with alleles compatible with typing
  filter(A_broad %in% typing_alleles | A_split %in% typing_alleles,
         B_broad %in% typing_alleles | B_split %in% typing_alleles,
         DRB1_broad %in% typing_alleles | DRB1_split %in% typing_alleles)

# combine all permutations of haplos (i.e., make phased genotypes)
cross_join(comp_haplos, comp_haplos, suffix = c("_1", "_2")) |>
  # make genotype out of serological equivalents of alleles in phased genotype
  mutate(genotype_ser = pmap(
    .l = list(A_broad_1, A_broad_2, A_split_1, A_split_2,
              B_broad_1, B_broad_2, B_split_1, B_split_2,
              DRB1_broad_1, DRB1_broad_2, DRB1_split_1, DRB1_split_2),
    .f = c
  )) |> # keep only phased genotypes that are same as input typing
  filter(map_lgl(genotype_ser, \(x) all(typing_alleles %in% x))) |>
  mutate(
    phased_freq = 2 * EURCAU_freq_1 * EURCAU_freq_2,
    likelihood = phased_freq / sum(phased_freq)
  ) |>
  slice_max(likelihood)


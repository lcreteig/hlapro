library(tidyverse)
library(readxl)
library(hlapro)

typing <- "A1 A2 B7 B8 DR1 DR7"

typing_alleles <- extract_alleles_str(typing, strip_locus = FALSE, loci = c("A", "B", "DRB1"))

haplo_df <- read_xlsx("~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx")

# compatible_haplos <- haplo_df |>
#   select(A, C, B, `DRB3-4-5`, DRB1, DQB1, EURCAU_freq, EURCAU_rank) |>
#   filter(get_serology(A) %in% typing_alleles[c("A_1", "A_2")],
#          get_serology(B) %in% typing_alleles[c("B_1", "B_2")],
#          get_serology(DRB1) %in% typing_alleles[c("DRB1_1", "DRB1_2")]) |>
#   mutate(haplo = str_c(A, B, C, `DRB3-4-5`, DRB1, DQB1, sep = "~"), .before = EURCAU_freq)
#
# expand(compatible_haplos, haplo, haplo)

top_haplos_ser <- haplo_df |>
  slice_max(EURCAU_freq, n = 5000) |>
  mutate(haplo = str_c(A, B, C, `DRB3-4-5`, DRB1, DQB1, sep = "~"), .before = EURCAU_freq) |>
  select(A, B, DRB1, haplo, EURCAU_freq, EURCAU_rank) |>
  mutate(across(c(A, B, DRB1), ~ get_serology(.x)))

comp_haplos <- top_haplos_ser |>
  filter(
    A %in% typing_alleles[c("A_1", "A_2")],
    B %in% typing_alleles[c("B_1", "B_2")],
    DRB1 %in% typing_alleles[c("DRB1_1", "DRB1_2")]
  )

cross_join(comp_haplos, comp_haplos, suffix = c("_1", "_2")) |>
  mutate(genotype_ser = pmap(
    .l = list(A_1, A_2, B_1, B_2, DRB1_1, DRB1_2),
    .f = c
  )) |>
  filter(map_lgl(genotype_ser, \(x) all(typing_alleles %in% x))) |>
  mutate(
    phased_freq = 2 * EURCAU_freq_1 * EURCAU_freq_2,
    likelihood = phased_freq / sum(phased_freq)
  ) |>
  slice_max(likelihood)

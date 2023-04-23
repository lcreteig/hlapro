library(tidyverse)

test_df <- tribble(
  ~ETnumber, ~typing,
  "ET1",   "A*01",
  "ET1",   "A2",
  "ET2",   "A*29:01:01",
  "ET2",   "C*233:01:01",
  "ET3",   "A*01:08",
  "ET3",   "B*07:02",
)

test_df_notr <- tribble(
  ~ETnumber, ~notr_typing,
  "ET1",   "A1 A2 B7 B8",
  "ET2",   "A29 A19 B7",
  "ET3",   "A1 B7 Bw6"
)

test_df |>
  mutate(etrl = etrl_lookup(typing)) |>
  unnest(cols = c(etrl)) |>
  rename(allele = `Allele`,
         split = `ET MatchDeterminantSplit`,
         broad = `ET MatchDeterminantBroad`,
         public = `Public`) |>
  mutate(other_alleles = if_else(is.na(allele), typing, NA)) |>
  unite("etrl_typing", c(split, broad, public, other_alleles), sep = " ", remove = FALSE, na.rm = TRUE) |>
  group_by(ETnumber) |>
  reframe(etrl_typing = str_flatten(etrl_typing, collapse = " ")) |>
  left_join(test_df_notr) |>
  mutate(etrl_typing = str_split(etrl_typing, " "),
         notr_typing = str_split(notr_typing, " ")) |>
  rowwise() |>
  mutate(typings_identical = setequal(etrl_typing, notr_typing),
         etrl_only = list(setdiff(notr_typing, etrl_typing)),
         notr_only = list(setdiff(etrl_typing, notr_typing))
         ) |>
  mutate(across(c(etrl_typing, notr_typing, etrl_only, notr_only), ~ str_flatten(.x, collapse = " ")))


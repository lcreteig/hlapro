# TODO: add test for homozygous haplotype frequency
# TODO: add test for phased genotype frequencies (ground truth from HaploStats)
# TODO: add test for outcome genotype (ground truth from HaploStats)
test_that("serological typings return expected EURCAU genotype frequencies", {
  # NMDP data cannot be made available; users have to ID and accept license
  skip_on_cran()
  skip_on_ci()

  input_typings <- c(
    "A24 A28 B35 B61 DR4 DR11",
    "A2 A3 B52 B35 Cw4 DR11 DR52 DQ3"
  )

  geno_freqs <- upscale_typings(
    filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
    typing = input_typings,
    n_haplos = 5000
  )

  freqs <- geno_freqs |>
    dplyr::group_by(id_input_typing) |>
    dplyr::slice_max(unphased_prob, n = 1, with_ties = FALSE) |>
    dplyr::pull(unphased_freq)

  expect_equal(freqs, c(2.542e-07, 5.13531e-08), tolerance = 1e-03)
})

test_that("serological typings return expected EURCAU genotype frequencies", {
  # NMDP data cannot be made available; users have to ID and accept license
  skip_on_cran()
  skip_on_ci()

  input_typing_1 <- "A24 A28 B35 B61 DR4 DR11"
  geno_freq <- upscale_typing(
    filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
    typing = input_typing_1,
    n_haplos = 5000
  )
  expect_equal(geno_freq$unphased_freq[1], 2.542e-07, tolerance = 1e-03)

  input_typing_2 <- "A2 A3 B52 B35 Cw4 DR11 DR52 DQ3"
  geno_freq <- upscale_typing(
    filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
    typing = input_typing_2,
    n_haplos = 5000
  )
  expect_equal(geno_freq$unphased_freq[1], 5.13531e-08, tolerance = 1e-05)
})

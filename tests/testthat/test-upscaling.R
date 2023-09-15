# TODO: add test for phased genotype frequencies (ground truth from HaploStats)
test_that("serological typings return expected high-res EURCAU genotypes", {
  # Compared to PROCARE 1.0 data

  # NMDP data cannot be made available; users have to ID and accept license
  skip_on_cran()
  skip_on_ci()

  input_typings <- c(
    "A24 A28 B35 B61 DR4 DR11",
    "A2 A3 B52 B35 Cw4 DR11 DR52 DQ3"
  )

  res <- upscale_typings(
    filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
    typing = input_typings,
    n_haplos = 5000
  ) |>
    dplyr::group_by(id_input_typing) |>
    dplyr::slice_max(unphased_prob, n = 1, with_ties = FALSE)

  # frequency
  expect_equal(dplyr::pull(res, unphased_freq),
    c(2.542e-07, 5.13531e-08),
    tolerance = 1e-03
  )
  # genotype
  expect_equal(
    dplyr::pull(res, unphased_geno),
    c(
      paste(
        "A*24:02g A*68:01g B*35:03g B*40:02g C*02:02g C*04:01g",
        "DQB1*03:01g DQB1*03:02g DRB1*04:03 DRB1*11:01g DRB3*02:02g DRB4*01:01g"
      ),
      paste(
        "A*02:01g A*03:01g B*35:01g B*52:01g C*04:01g C*12:02",
        "DQB1*03:01g DQB1*03:01g DRB1*11:01g DRB1*11:04 DRB3*02:02g DRB3*02:02g"
      )
    )
  )
})

test_that("homozygous low-res genotype returns homozygous high-res genotypes", {
  # Compared to HaploStats output

  # NMDP data cannot be made available; users have to ID and accept license
  skip_on_cran()
  skip_on_ci()

  input_typings <- c(
    "A1 B8 Cw7 DR17 DR52 DQ2",
    "A2 B7 DR15 DQ1"
  )

  res <- upscale_typings(
    filepath = "~/Downloads/A~C~B~DRB3-4-5~DRB1~DQB1.xlsx",
    typing = input_typings
  ) |>
    dplyr::group_by(id_input_typing) |>
    dplyr::slice_max(unphased_prob, n = 1, with_ties = FALSE)

  # frequency
  expect_equal(dplyr::pull(res, unphased_freq),
    c(4.259e-03, 3.896e-04),
    tolerance = 1e-03
  )
  # genotype
  expect_equal(
    dplyr::pull(res, unphased_geno),
    c(
      paste(
        "A*01:01g A*01:01g B*08:01g B*08:01g C*07:01g C*07:01g",
        "DQB1*02:01g DQB1*02:01g DRB1*03:01 DRB1*03:01 DRB3*01:01 DRB3*01:01"
      ),
      paste(
        "A*02:01g A*02:01g B*07:02g B*07:02g C*07:02g C*07:02g",
        "DQB1*06:02 DQB1*06:02 DRB1*15:01 DRB1*15:01 DRB5*01:01 DRB5*01:01"
      )
    )
  )
})

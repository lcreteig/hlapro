test_that("both HLA-A alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(A_1 = "1", A_2 = "2", df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("both HLA-A alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQ5 A1 A2")
  df_out <- tidyr::tibble(A_1 = "1", A_2 = "2", df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("both HLA-A alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "B7 B8 A1 A2 DR1 DQ5")
  df_out <- tidyr::tibble(A_1 = "1", A_2 = "2", df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("both HLA-A alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "A1 B7 B8 A2 DR1 DQ5")
  df_out <- tidyr::tibble(A_1 = "1", A_2 = "2", df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("second allele is NA when only 1 HLA-A allele is present", {
  df_in <- tidyr::tibble(typing = "B7 A1 B8")
  df_out <- tidyr::tibble(A_1 = "1", A_2 = NA_character_, df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A1 B7 B8")
  df_out <- tidyr::tibble(A_1 = "1", A_2 = NA_character_, df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
  df_in <- tidyr::tibble(typing = "B7 B8 A1")
  df_out <- tidyr::tibble(A_1 = "1", A_2 = NA_character_, df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("both alleles are NA when no HLA-A allele is present", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQ5")
  # TODO: consider turning this into a true NA in the function itself
  df_out <- tidyr::tibble(A_1 = "", A_2 = NA_character_, df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("double-digit HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A10 A19 B7")
  df_out <- tidyr::tibble(A_1 = "10", A_2 = "19", df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("high resolution HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*01:01 A*01:02 B07*01:02 B07*01:02")
  df_out <- tidyr::tibble(A_1 = "01:01", A_2 = "01:02", df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("intermediate resolution HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC B*08:NMTJ")
  df_out <- tidyr::tibble(A_1 = "02:HBMC", A_2 = NA_character_, df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("'A's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC DRB1*04:AMR B*08:NMTJ")
  df_out <- tidyr::tibble(A_1 = "02:HBMC", A_2 = NA_character_, df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A*02:HBMC B*08:BAD")
  df_out <- tidyr::tibble(A_1 = "02:HBMC", A_2 = NA_character_, df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A*01:AABJE A*02:HBMC")
  df_out <- tidyr::tibble(A_1 = "01:AABJE", A_2 = "02:HBMC", df_in)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

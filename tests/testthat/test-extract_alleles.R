# HLA-A -------------------------------------------------------------------

test_that("both HLA-A alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("both HLA-A alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQ5 A1 A2")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("both HLA-A alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "B7 B8 A1 A2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("both HLA-A alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "A1 B7 B8 A2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("second allele is NA when only 1 HLA-A allele is present", {
  df_in <- tidyr::tibble(typing = "B7 A1 B8")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A1 B7 B8")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
  df_in <- tidyr::tibble(typing = "B7 B8 A1")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("both alleles are NA when no HLA-A allele is present", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQ5")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, A_1 = NA_character_, A_2 = NA_character_)
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("double-digit HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A10 A19 B7")
  df_out <- tidyr::tibble(df_in, A_1 = "10", A_2 = "19")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("high resolution HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*01:01 A*01:02 B*07:02 B*07:02")
  df_out <- tidyr::tibble(df_in, A_1 = "01:01", A_2 = "01:02")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("A's in HLA- prefix are not extracted", {
  df_in <- tidyr::tibble(typing = "HLA-DRB1*04:01 A*01:01 A*01:02 B*07:02")
  df_out <- tidyr::tibble(df_in, A_1 = "01:01", A_2 = "01:02")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("intermediate resolution HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, A_1 = "02:HBMC", A_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

test_that("'A's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC DRB1*04:AMR B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, A_1 = "02:HBMC", A_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A*02:HBMC B*08:BAD")
  df_out <- tidyr::tibble(df_in, A_1 = "02:HBMC", A_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A*01:AABJE A*02:HBMC")
  df_out <- tidyr::tibble(df_in, A_1 = "01:AABJE", A_2 = "02:HBMC")
  expect_equal(extract_alleles(df_in, typing, locus = "A"), df_out)
})

# HLA-B -------------------------------------------------------------------

test_that("public is ignored", {
  df_in <- tidyr::tibble(typing = "A1 A2 Bw4 BW6 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("both HLA-B alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "B7 B8 A1 A2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("both HLA-B alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "A1 A2 DR1 DQ5 B7 B8")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("both HLA-B alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("both HLA-B alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "B7 A1 B8 A2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("second allele is NA when only 1 HLA-B allele is present", {
  df_in <- tidyr::tibble(typing = "A1 B7 A2")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
  df_in <- tidyr::tibble(typing = "B7 A1 A2")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
  df_in <- tidyr::tibble(typing = "A1 A2 B7")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("both alleles are NA when no HLA-B allele is present", {
  df_in <- tidyr::tibble(typing = "A1 A2 DR1 DQ5")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, B_1 = NA_character_, B_2 = NA_character_)
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("double-digit HLA-B alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A10 B15 B40")
  df_out <- tidyr::tibble(df_in, B_1 = "15", B_2 = "40")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("high resolution HLA-B alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*01:01 A*01:02 B*07:02 B*07:03")
  df_out <- tidyr::tibble(df_in, B_1 = "07:02", B_2 = "07:03")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("B's in other loci are not extracted", {
  df_in <- tidyr::tibble(typing = "DRB1*04:01 A*01:01 A*01:02 B*07:02 B*07:03")
  df_out <- tidyr::tibble(df_in, B_1 = "07:02", B_2 = "07:03")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("intermediate resolution HLA-B alleles are extracted", {
  df_in <- tidyr::tibble(typing = "B*08:NMTJ A*02:HBMC ")
  df_out <- tidyr::tibble(df_in, B_1 = "08:NMTJ", B_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

test_that("'B's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC DRB1*04:AMR B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, B_1 = "08:NMTJ", B_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
  df_in <- tidyr::tibble(typing = "C*01:BJZ B*08:BAD")
  df_out <- tidyr::tibble(df_in, B_1 = "08:BAD", B_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
  df_in <- tidyr::tibble(typing = "B*08:BAD B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, B_1 = "08:BAD", B_2 = "08:NMTJ")
  expect_equal(extract_alleles(df_in, typing, locus = "B"), df_out)
})

# HLA-C -------------------------------------------------------------------

test_that("both HLA-C alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "Cw1 Cw2 A1 A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("both HLA-C alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DQ5 Cw1 Cw2")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("both HLA-C alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw1 Cw2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("both HLA-C alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "A1 A2 Cw1 B7 B8 Cw2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("second allele is NA when only 1 HLA-C allele is present", {
  df_in <- tidyr::tibble(typing = "A1 Cw1 B7 B8")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
  df_in <- tidyr::tibble(typing = "Cw1 A1 B7 B8")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
  df_in <- tidyr::tibble(typing = "A1 B7 B8 Cw1")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("both alleles are NA when no HLA-C allele is present", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQ5")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, C_1 = NA_character_, C_2 = NA_character_)
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("double-digit HLA-C alleles are extracted", {
  df_in <- tidyr::tibble(typing = "Cw10 Cw10 A10 A19 B7")
  df_out <- tidyr::tibble(df_in, C_1 = "10", C_2 = "10")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("high resolution HLA-C alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*01:01 A*01:02 B*07:02C*01:02 C*02:08")
  df_out <- tidyr::tibble(df_in, C_1 = "01:02", C_2 = "02:08")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("intermediate resolution HLA-C alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC B*08:NMTJ C*01:BJZ")
  df_out <- tidyr::tibble(df_in, C_1 = "01:BJZ", C_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

test_that("'C's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "C*01:BJZ A*02:HBMC DRB1*07:CYMD B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, C_1 = "01:BJZ", C_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
  df_in <- tidyr::tibble(typing = "C*01:BJZ DRB1*07:GCM")
  df_out <- tidyr::tibble(df_in, C_1 = "01:BJZ", C_2 = "")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
  df_in <- tidyr::tibble(typing = "C*04:CYMD C*01:BJZ")
  df_out <- tidyr::tibble(df_in, C_1 = "04:CYMD", C_2 = "01:BJZ")
  expect_equal(extract_alleles(df_in, typing, locus = "C"), df_out)
})

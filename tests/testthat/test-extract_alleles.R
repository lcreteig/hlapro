# Vectorization -----------------------------------------------------------

test_that("multiple loci can be extracted", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2", B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = c("A", "B")), df_out)
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2",
                          B_1 = "7", B_2 = "8",
                          C_1 = "3", C_2 = "",
                          DPB1_1 = NA_character_,
                          DPB1_2 = NA_character_,
                          DQA1_1 = NA_character_,
                          DQA1_2 = NA_character_,
                          DQB1_1 = "5", DQB1_2 = "8",
                          DRB1_1 = "4", DRB1_2 = "11",
                          DRB._1 = "52", DRB._2 = "53")
  expect_equal(extract_alleles_df(df_in, "typing"), df_out)
})

# String issues -----------------------------------------------------------

test_that("extraction is robust to double spaces", {
  df_in <- tidyr::tibble(typing = "A1  A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  df_out["typing"] <- stringr::str_squish(df_out["typing"])
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("extraction is robust to tabs", {
  df_in <- tidyr::tibble(typing = "A1   A2 B7 B8 Cw1 Cw2 ")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  df_out["typing"] <- stringr::str_squish(df_out["typing"])
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})


test_that("extraction is robust to trailing spaces", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw1 Cw2 ")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  df_out["typing"] <- stringr::str_squish(df_out["typing"])
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("extraction is robust to leading spaces", {
  df_in <- tidyr::tibble(typing = " A1 A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  df_out["typing"] <- stringr::str_squish(df_out["typing"])
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

# HLA-A -------------------------------------------------------------------

test_that("both HLA-A alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("both HLA-A alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQ5 A1 A2")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("both HLA-A alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "B7 B8 A1 A2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("both HLA-A alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "A1 B7 B8 A2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "2")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("second allele is NA when only 1 HLA-A allele is present", {
  df_in <- tidyr::tibble(typing = "B7 A1 B8")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A1 B7 B8")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
  df_in <- tidyr::tibble(typing = "B7 B8 A1")
  df_out <- tidyr::tibble(df_in, A_1 = "1", A_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("both alleles are NA when no HLA-A allele is present", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQ5")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, A_1 = NA_character_, A_2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("double-digit HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A10 A19 B7")
  df_out <- tidyr::tibble(df_in, A_1 = "10", A_2 = "19")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("high resolution HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*01:01 A*02:101:01:02N B*07:02 B*07:02")
  df_out <- tidyr::tibble(df_in, A_1 = "01:01", A_2 = "02:101:01:02N")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("A's in HLA- prefix are not extracted", {
  df_in <- tidyr::tibble(typing = "HLA-DRB1*04:01 A*01:01 A*01:02 B*07:02")
  df_out <- tidyr::tibble(df_in, A_1 = "01:01", A_2 = "01:02")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("intermediate resolution HLA-A alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, A_1 = "02:HBMC", A_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

test_that("'A's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC DRB1*04:AMR B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, A_1 = "02:HBMC", A_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A*02:HBMC B*08:BAD")
  df_out <- tidyr::tibble(df_in, A_1 = "02:HBMC", A_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
  df_in <- tidyr::tibble(typing = "A*01:AABJE A*02:HBMC")
  df_out <- tidyr::tibble(df_in, A_1 = "01:AABJE", A_2 = "02:HBMC")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "A"), df_out)
})

# HLA-B -------------------------------------------------------------------

test_that("public is ignored", {
  df_in <- tidyr::tibble(typing = "A1 A2 Bw4 BW6 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("both HLA-B alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "B7 B8 A1 A2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("both HLA-B alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "A1 A2 DR1 DQ5 B7 B8")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("both HLA-B alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("both HLA-B alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "B7 A1 B8 A2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("second allele is NA when only 1 HLA-B allele is present", {
  df_in <- tidyr::tibble(typing = "A1 B7 A2")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
  df_in <- tidyr::tibble(typing = "B7 A1 A2")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
  df_in <- tidyr::tibble(typing = "A1 A2 B7")
  df_out <- tidyr::tibble(df_in, B_1 = "7", B_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("both alleles are NA when no HLA-B allele is present", {
  df_in <- tidyr::tibble(typing = "A1 A2 DR1 DQ5")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, B_1 = NA_character_, B_2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("double-digit HLA-B alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A10 B15 B40")
  df_out <- tidyr::tibble(df_in, B_1 = "15", B_2 = "40")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("high resolution HLA-B alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*01:01 A*01:02 B*07:02 B*07:03:01:02N")
  df_out <- tidyr::tibble(df_in, B_1 = "07:02", B_2 = "07:03:01:02N")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("B's in other loci are not extracted", {
  df_in <- tidyr::tibble(typing = "DRB1*04:01 A*01:01 A*01:02 B*07:02 B*07:03")
  df_out <- tidyr::tibble(df_in, B_1 = "07:02", B_2 = "07:03")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("intermediate resolution HLA-B alleles are extracted", {
  df_in <- tidyr::tibble(typing = "B*08:NMTJ A*02:HBMC")
  df_out <- tidyr::tibble(df_in, B_1 = "08:NMTJ", B_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

test_that("'B's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC DRB1*04:AMR B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, B_1 = "08:NMTJ", B_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
  df_in <- tidyr::tibble(typing = "C*01:BJZ B*08:BAD")
  df_out <- tidyr::tibble(df_in, B_1 = "08:BAD", B_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
  df_in <- tidyr::tibble(typing = "B*08:BAD B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, B_1 = "08:BAD", B_2 = "08:NMTJ")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "B"), df_out)
})

# HLA-C -------------------------------------------------------------------

test_that("both HLA-C alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "Cw1 Cw2 A1 A2 B7 B8 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("both HLA-C alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DQ5 Cw1 Cw2")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("both HLA-C alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw1 Cw2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("both HLA-C alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "A1 A2 Cw1 B7 B8 Cw2 DR1 DQ5")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "2")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("second allele is NA when only 1 HLA-C allele is present", {
  df_in <- tidyr::tibble(typing = "A1 Cw1 B7 B8")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
  df_in <- tidyr::tibble(typing = "Cw1 A1 B7 B8")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
  df_in <- tidyr::tibble(typing = "A1 B7 B8 Cw1")
  df_out <- tidyr::tibble(df_in, C_1 = "1", C_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("both alleles are NA when no HLA-C allele is present", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQ5")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, C_1 = NA_character_, C_2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("double-digit HLA-C alleles are extracted", {
  df_in <- tidyr::tibble(typing = "Cw10 Cw10 A10 A19 B7")
  df_out <- tidyr::tibble(df_in, C_1 = "10", C_2 = "10")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("high resolution HLA-C alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*01:01 B*07:02 C*01:02 C*02:03:01:02N")
  df_out <- tidyr::tibble(df_in, C_1 = "01:02", C_2 = "02:03:01:02N")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("intermediate resolution HLA-C alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*02:HBMC B*08:NMTJ C*01:BJZ")
  df_out <- tidyr::tibble(df_in, C_1 = "01:BJZ", C_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

test_that("'C's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "C*01:BJZ A*02:HBMC DRB1*07:CYMD B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, C_1 = "01:BJZ", C_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
  df_in <- tidyr::tibble(typing = "C*01:BJZ DRB1*07:GCM")
  df_out <- tidyr::tibble(df_in, C_1 = "01:BJZ", C_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
  df_in <- tidyr::tibble(typing = "C*04:CYMD C*01:BJZ")
  df_out <- tidyr::tibble(df_in, C_1 = "04:CYMD", C_2 = "01:BJZ")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "C"), df_out)
})

# HLA-DPB1 ----------------------------------------------------------------

test_that("both HLA-DPB1 alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "DPB1*04:02 DPB1*05:01 B*39:09 DRB1*04:11")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "04:02", DPB1_2 = "05:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
})

test_that("both HLA-DPB1 alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "B*39:09 DRB1*04:11 DPB1*04:02 DPB1*05:01")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "04:02", DPB1_2 = "05:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
})

test_that("both HLA-DPB1 alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "B*39:09 DPB1*04:02 DPB1*05:01 DRB1*04:11")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "04:02", DPB1_2 = "05:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
})

test_that("both HLA-DPB1 alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "DPB1*04:02 DPA1*02:02 DPB1*05:01 B*39:09")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "04:02", DPB1_2 = "05:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
})

test_that("second allele is NA when only 1 HLA-DPB1 allele is present", {
  df_in <- tidyr::tibble(typing = "DPB1*01:01 DPA1*02:02 B*39:09")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "01:01", DPB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
  df_in <- tidyr::tibble(typing = "DPA1*02:02 DPB1*01:01 B*39:09")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "01:01", DPB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
  df_in <- tidyr::tibble(typing = "DPA1*02:02 B*39:09 DPB1*01:01")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "01:01", DPB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
})

test_that("both alleles are NA when no HLA-DPB1 allele is present", {
  df_in <- tidyr::tibble(typing = "DPA1*02:02 DQB1*03:24 DRB1*04:11")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, DPB1_1 = NA_character_, DPB1_2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
})

test_that("more than single-digit HLA-DPB1 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "DPB1*129:01 DPB1*15:01:01:01")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "129:01", DPB1_2 = "15:01:01:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
})

test_that("intermediate resolution HLA-DPB1 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*02:BMC DPB1*15:FNWN B*08:DPB DPB1*04:BDVU")
  df_out <- tidyr::tibble(df_in, DPB1_1 = "15:FNWN", DPB1_2 = "04:BDVU")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DPB1"), df_out)
})

# HLA-DQA1 ----------------------------------------------------------------

test_that("both HLA-DQA1 alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "DQA1*01:03 DQA1*05:05 DPB1*05:01 DQB1*03:01")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "01:03", DQA1_2 = "05:05")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
})

test_that("both HLA-DQA1 alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "DPB1*05:01 DQB1*03:01 DQA1*01:03 DQA1*05:05")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "01:03", DQA1_2 = "05:05")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
})

test_that("both HLA-DQA1 alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "DPB1*05:01 DQA1*01:03 DQA1*05:05 DQB1*03:01")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "01:03", DQA1_2 = "05:05")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
})

test_that("both HLA-DQA1 alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "DPB1*05:01 DQA1*01:03 DQB1*03:01 DQA1*05:05")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "01:03", DQA1_2 = "05:05")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
})

test_that("second allele is NA when only 1 HLA-DQA1 allele is present", {
  df_in <- tidyr::tibble(typing = "DQA1*01:03 DPA1*02:02 B*39:09 DQB1*03:01")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "01:03", DQA1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
  df_in <- tidyr::tibble(typing = "DPA1*02:02 B*39:09 DQA1*01:03 DQB1*03:01")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "01:03", DQA1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
  df_in <- tidyr::tibble(typing = "DPA1*02:02 B*39:09 DQB1*03:01 DQA1*01:03")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "01:03", DQA1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
})

# TODO: add serology
test_that("both alleles are NA when no HLA-DQA1 allele is present", {
  df_in <- tidyr::tibble(typing = "DPA1*02:02 DQB1*03:24 DQB1*03:01 DRB1*04:11")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, DQA1_1 = NA_character_, DQA1_2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DQ3 DQ5")
  df_out <- tidyr::tibble(df_in, DQA1_1 = NA_character_, DQA1_2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
})

test_that("more than single-digit HLA-DQA1 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "DQA1*01:01:01 DQA1*01:02:01")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "01:01:01", DQA1_2 = "01:02:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
})

test_that("intermediate resolution HLA-DQA1 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*02:BMC DQA1*05:BCHD B*08:DQA DQA1*01:UKJA")
  df_out <- tidyr::tibble(df_in, DQA1_1 = "05:BCHD", DQA1_2 = "01:UKJA")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQA1"), df_out)
})

# HLA-DQB1 ----------------------------------------------------------------

test_that("both HLA-DQB1 alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "DQ5 DQ8 A1 A2 B7 B8 Cw1 Cw2 DR1")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "5", DQB1_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

test_that("both HLA-DQB1 alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw1 Cw2 DR1 DQ5 DQ8")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "5", DQB1_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

test_that("both HLA-DQB1 alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw1 Cw2 DQ5 DQ8 DR1")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "5", DQB1_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

test_that("both HLA-DQB1 alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DQ5 Cw1 Cw2 DQ8 DR1")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "5", DQB1_2 = "8")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

test_that("second allele is NA when only 1 HLA-DQB1 allele is present", {
  df_in <- tidyr::tibble(typing = "A1 Cw1 B7 DQ5")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "5", DQB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
  df_in <- tidyr::tibble(typing = "DQ5 A1 B7 B8")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "5", DQB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
  df_in <- tidyr::tibble(typing = "A1 B7 DQ5 Cw1")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "5", DQB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

test_that("both alleles are NA when no HLA-DQB1 allele is present", {
  df_in <- tidyr::tibble(typing = "B7 B8 DR1 DQA1*01:03")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, DQB1_1 = NA_character_, DQB1_2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

test_that("high resolution HLA-DQB1 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "A*01:01 DQA1*01:03 DQB1*03:01 DQB1*05:01")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "03:01", DQB1_2 = "05:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

test_that("intermediate resolution HLA-DQB1 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "DQA1*05:BC B*08:NM DQB1*03:AG DQB1*03:AFYYJ")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "03:AG", DQB1_2 = "03:AFYYJ")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

test_that("'DQ's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "DQB1*03:AG A*02:DQMC DRB1*07:CYMD B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "03:AG", DQB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
  df_in <- tidyr::tibble(typing = "DQB1*03:AFYYJ DRB1*07:GDQM")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "03:AFYYJ", DQB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
  df_in <- tidyr::tibble(typing = "C*04:CYMDQ DQB1*03:AG C*01:BZ DQB1*03:AFYYJ")
  df_out <- tidyr::tibble(df_in, DQB1_1 = "03:AG", DQB1_2 = "03:AFYYJ")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DQB1"), df_out)
})

# HLA-DRB1 ----------------------------------------------------------------

test_that("both HLA-DRB1 alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "DR1 DR3 A1 A2 B7 B8 Cw1 Cw2 DQ5 DQ8 DR52")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "1", DRB1_2 = "3")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("both HLA-DRB1 alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 Cw1 Cw2 DQ5 DQ8 DR52 DR1 DR3")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "1", DRB1_2 = "3")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("both HLA-DRB1 alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DR3 Cw1 Cw2 DQ5 DQ8 DR52")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "1", DRB1_2 = "3")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("both HLA-DRB1 alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "A1 A2 DR1 DR52 B7 B8 DR3 Cw1 Cw2 DQ5 DQ8")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "1", DRB1_2 = "3")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("second allele is NA when only 1 HLA-DRB1 allele is present", {
  df_in <- tidyr::tibble(typing = "A1 Cw1 B7 DQ5 DR1")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "1", DRB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
  df_in <- tidyr::tibble(typing = "DR1 A1 B7 B8 DQ5")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "1", DRB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
  df_in <- tidyr::tibble(typing = "A1 B7 DR1 DQ5 Cw1")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "1", DRB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("both alleles are NA when no HLA-DRB1 allele is present", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DQA1*01:03 DQ3")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, DRB1_1 = NA_character_, DRB1_2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("high resolution HLA-DRB1 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "DRA*01:01 DRB1*03:01 DQA1*01:03 DRB1*07:01")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "03:01", DRB1_2 = "07:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("intermediate resolution HLA-DRB1 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "DQA1*05:BC B*08:NM DRB1*04:AMR DRB1*07:GC")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "04:AMR", DRB1_2 = "07:GC")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("'DR's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "DQB1*03:AG A*02:DRMC DRB1*07:CYMD B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "07:CYMD", DRB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
  df_in <- tidyr::tibble(typing = "DRB1*07:GC DQB1*07:GDRM")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "07:GC", DRB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
  df_in <- tidyr::tibble(typing = "C*04:CYMDR DRB1*07:GC C*01:BZ DRB1*04:AMR")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "07:GC", DRB1_2 = "04:AMR")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

test_that("DRB3/4/5 is not extracted", {
  df_in <- tidyr::tibble(typing = "DR51 DR52 DR53 DR1 DR5")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "1", DRB1_2 = "5")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
  df_in <- tidyr::tibble(typing = "DRB3*04:01 DRB4*03:01 DRB1*02:01 DRB5*01:01")
  df_out <- tidyr::tibble(df_in, DRB1_1 = "02:01", DRB1_2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB1"), df_out)
})

# HLA-DRB3/4/5 ------------------------------------------------------------

test_that("all variants are extracted", {
  df_in <- tidyr::tibble(typing = "DR51 DR52")
  df_out <- tidyr::tibble(df_in, DRB._1 = "51", DRB._2 = "52")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
  df_in <- tidyr::tibble(typing = "DR53 DRB3*01:32")
  df_out <- tidyr::tibble(df_in, DRB._1 = "53", DRB._2 = "3*01:32")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
  df_in <- tidyr::tibble(typing = "DRB4*01:05 DRB5*02:04")
  df_out <- tidyr::tibble(df_in, DRB._1 = "4*01:05", DRB._2 = "5*02:04")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("DRB1 is not extracted", {
  df_in <- tidyr::tibble(typing = "DR5 DR52 DRB1*01:01 DRB4*04:01")
  df_out <- tidyr::tibble(df_in, DRB._1 = "52", DRB._2 = "4*04:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("negatives are not extracted", {
  df_in <- tidyr::tibble(typing = "DRB3*Neg DRB4*Neg DRB5*Neg DRB4*04:01")
  df_out <- tidyr::tibble(df_in, DRB._1 = "4*04:01", DRB._2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("both HLA-DRB3/4/5 alleles are extracted from beginning", {
  df_in <- tidyr::tibble(typing = "DR52 DR53 A1 A2 B7 B8 DQ5 DQ8 DR4 DR11")
  df_out <- tidyr::tibble(df_in, DRB._1 = "52", DRB._2 = "53")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("both HLA-DRB3/4/5 alleles are extracted from end", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DQ5 DQ8 DR4 DR11 DR52 DR53")
  df_out <- tidyr::tibble(df_in, DRB._1 = "52", DRB._2 = "53")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("both HLA-DRB3/4/5 alleles are extracted from middle", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR52 DR53 DQ5 DQ8 DR4 DR11")
  df_out <- tidyr::tibble(df_in, DRB._1 = "52", DRB._2 = "53")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("both HLA-DRB3/4/5 alleles are extracted when separated by others", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR52 DQ5 DQ8 DR4 DR53 DR11")
  df_out <- tidyr::tibble(df_in, DRB._1 = "52", DRB._2 = "53")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("second allele is NA when only 1 HLA-DRB3/4/5 allele is present", {
  df_in <- tidyr::tibble(typing = "A1 Cw1 B7 DQ5 DR5 DR52")
  df_out <- tidyr::tibble(df_in, DRB._1 = "52", DRB._2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
  df_in <- tidyr::tibble(typing = "A1 Cw1 B7 DR52 DQ5 DR5")
  df_out <- tidyr::tibble(df_in, DRB._1 = "52", DRB._2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
  df_in <- tidyr::tibble(typing = "DR52 A1 Cw1 B7 DQ5 DR5")
  df_out <- tidyr::tibble(df_in, DRB._1 = "52", DRB._2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("both alleles are NA when no HLA-DRB3/4/5 allele is present", {
  df_in <- tidyr::tibble(typing = "A1 A2 B7 B8 DR1 DR5 DQA1*01:03 DQ3")
  # TODO: empty string if only no 2nd allele. Useful that it's NA here?
  df_out <- tidyr::tibble(df_in, DRB._1 = NA_character_, DRB._2 = NA_character_)
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("high resolution HLA-DRB3/4/5 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "DRB5*02:03 DRB1*05:01 DRB1*03:01 DRB3*01:01")
  df_out <- tidyr::tibble(df_in, DRB._1 = "5*02:03", DRB._2 = "3*01:01")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("intermediate resolution HLA-DRB3/4/5 alleles are extracted", {
  df_in <- tidyr::tibble(typing = "DRB5*03:MAC DRB1*04:AMR DRB3*05:MAC")
  df_out <- tidyr::tibble(df_in, DRB._1 = "5*03:MAC", DRB._2 = "3*05:MAC")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

test_that("'DRB's in MACs are not extracted", {
  df_in <- tidyr::tibble(typing = "DQB1*03:AG A*02:DRBC DRB5*03:MAC B*08:NMTJ")
  df_out <- tidyr::tibble(df_in, DRB._1 = "5*03:MAC", DRB._2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
  df_in <- tidyr::tibble(typing = "DRB5*03:MAC DQB1*07:DRBA")
  df_out <- tidyr::tibble(df_in, DRB._1 = "5*03:MAC", DRB._2 = "")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
  df_in <- tidyr::tibble(typing = "C*04:YDRB DRB5*03:MAC DRB1*07:GC DRB5*01:AM")
  df_out <- tidyr::tibble(df_in, DRB._1 = "5*03:MAC", DRB._2 = "5*01:AM")
  expect_equal(extract_alleles_df(df_in, "typing", loci = "DRB."), df_out)
})

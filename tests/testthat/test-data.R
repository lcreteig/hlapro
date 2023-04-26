test_that("etrl_version matches site", {
  expect_equal(attr(etrl_hla, "version"), "2.1")
})

test_that("etrl release date matches site", {
  expect_equal(attr(etrl_hla, "date"), "24-04-2023")
})

test_that("table has 4 columns", {
  expect_equal(length(etrl_hla), 4)
})

test_that("column names and types are correct", {
  etrl_hla_info <- c(
    Allele = "character",
    `ET MatchDeterminantSplit` = "character",
    `ET MatchDeterminantBroad` = "character",
    Public = "character"
  )
  expect_equal(purrr::map_chr(etrl_hla, class), etrl_hla_info)
})

test_that("table has no empty rows", {
  expect_true(all(rowSums(is.na(etrl_hla)) != ncol(etrl_hla)))
})

test_that("alleles span all loci", {
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(A\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(B\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(C\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(DRB1\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(DRB3\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(DRB4\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(DRB5\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(DQB1*\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(DQA1*\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(DPB1*\\*)")))
  expect_true(any(stringr::str_detect(etrl_hla$Allele, r"(DPA1*\\*)")))
})

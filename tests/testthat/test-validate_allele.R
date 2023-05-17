test_that("serological/antigen notation is valid", {
  expect_true(validate_allele("A2"))
  expect_true(validate_allele("A80"))
  expect_true(validate_allele("A203"))
  expect_false(validate_allele("X1")) # no valid locus
  expect_true(validate_allele("A*6603")) # fine (old nomenclature)
  expect_false(validate_allele("A*66035")) # too many digits
  expect_false(validate_allele("A&6603")) # no valid separator
  expect_true(validate_allele("B1"))
  expect_true(validate_allele("Bw6"))
  expect_true(validate_allele("B50"))
  expect_true(validate_allele("Cw3"))
  expect_true(validate_allele("DR8"))
  expect_false(validate_allele("DRB11")) # if there's a B, needs "*" notation
  expect_false(validate_allele("DRA1")) # needs "*" notation
  expect_true(validate_allele("DR10"))
  expect_true(validate_allele("DR52"))
  expect_true(validate_allele("DQ1"))
  expect_true(validate_allele("DQA-01")) # ETRL nomenclature
  expect_false(validate_allele("DQA101")) # if "A1", needs "*" after
  expect_false(validate_allele("DQB101")) # if "B1", needs "*" after
  expect_true(validate_allele("DP-01"))
  expect_false(validate_allele("DJ-01")) # no valid locus
  expect_true(validate_allele("DP-0201"))
  expect_false(validate_allele("DPA101")) # if "A1", needs "*" after
  expect_false(validate_allele("DPB101")) # if "B1", needs "*" after
})

test_that("XX codes are valid", {
  expect_true(validate_allele("A*02:XX"))
  expect_true(validate_allele("B*78:XX"))
  expect_true(validate_allele("C*17:XX"))
  expect_true(validate_allele("DRB1*09:XX"))
  expect_true(validate_allele("DQB1*04:XX"))
  expect_true(validate_allele("DQA1*05:XX"))
})

test_that("MACs are valid", {
  expect_true(validate_allele("A*01:AABJE"))
  expect_true(validate_allele("DRB1*07:GC"))
  expect_true(validate_allele("B*15CFRG")) # fine (old nomenclature)
  expect_true(validate_allele("A*01KG"))
  expect_true(validate_allele("DPB1*04BDVU"))
  expect_false(validate_allele("DPB1*04:BDVUABC")) # too long
  expect_false(validate_allele("DPB1*04:B")) # too short
})

test_that("prefixing is allowed", {
  expect_true(validate_allele("HLA-A*01:01"))
  expect_false(validate_allele("HLAA*01:01")) # no separator
  expect_true(validate_allele("HLA-DRB1*10:03"))
})

test_that("protein-level typings are valid", {
  expect_true(validate_allele("A*74:06"))
  expect_false(validate_allele("A*10001:01")) # too many digits in protein
  expect_true(validate_allele("B*81:01"))
  expect_true(validate_allele("C*01:02"))
  expect_true(validate_allele("DRA*01:01"))
  expect_true(validate_allele("DRB1*01:01"))
  expect_false(validate_allele("DRB*01:01")) # DR not disambiguated
  expect_true(validate_allele("DRB1*12:02"))
  expect_true(validate_allele("DQB1*02:01"))
  expect_true(validate_allele("DQA1*01:05"))
  expect_true(validate_allele("DPB1*01:01"))
  expect_true(validate_allele("DPA1*01:03"))
})

test_that(">2 field codes are valid", {
  expect_true(validate_allele("C*01:02:31"))
  expect_true(validate_allele("C*01:02:30:01"))
  expect_false(validate_allele("C*01:02:30:01:01")) # one field too many
  expect_false(validate_allele("C*01:02:3000:01:01")) # too many digits
  expect_false(validate_allele("C*01:02:30:1111")) # too many digits
})

test_that("suffixes are valid", {
  expect_true(validate_allele("A*01:27N"))
  expect_true(validate_allele("A*24:02:01:02L"))
  expect_true(validate_allele("B*44:02:01:02S"))
  expect_true(validate_allele("C*01:121Q"))
})

test_that("P groups are valid", {
  expect_true(validate_allele("A*01:01P"))
  expect_true(validate_allele("DPB1*184:01P"))
})

test_that("G groups are valid", {
  expect_true(validate_allele("A*01:01:01G"))
  expect_true(validate_allele("DPB1*879:01:01G"))
})

test_that("ambiguities are allowed", {
  expect_true(validate_allele("A*74:06/A*22:03N"))
  expect_false(validate_allele("A*74:06/B*22:03N")) # must be same locus
  expect_true(validate_allele("A*30:01/30:14/30:15"))
  expect_true(validate_allele("A*30:01/30:14L/30:15")) # with suffix
  expect_true(validate_allele("A*30:01/14/15L")) # shorthand
  expect_false(validate_allele("A*01:01/A*01:02N:01")) # suffix in middle
  expect_true(validate_allele("C*01:02/C*01:03/C*01:04/C*01:05/C*01:06"))
  expect_true(validate_allele("DRB1*11/13"))
})

test_that("no illegal (combinations of) characters", {
  expect_false(validate_allele("B42:02")) # missing "*"
})

test_that("vectors work", {
  expect_equal(
    validate_allele(c("A2", "A*01:AABJE", "A*24:02:01:02L", "not")),
    c(TRUE, TRUE, TRUE, FALSE)
  )
})

test_that("mutate in dataframe works", {
  df <- tidyr::tibble(alleles = c("A2", "A*01:AABJE", "A*24:02:01:02L", "not"))
  df <- dplyr::mutate(df, alleles_check = validate_allele(alleles))
  expect_equal(df$alleles_check, c(TRUE, TRUE, TRUE, FALSE))
})

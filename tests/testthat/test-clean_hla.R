# strip_redundant() -----------------------------------------------------------

test_that("split(broad) works", {
  expect_equal(strip_redundant("A24(9)"), "A24")
})

test_that("broad(split) works", {
  expect_equal(strip_redundant("A9(24)"), "A24")
})

test_that("broads without split are not removed", {
  expect_equal(strip_redundant("A10"), "A10")
})

test_that("homozygotes are retained", {
  expect_equal(strip_redundant("A29 A29"), "A29 A29")
})

test_that("splits are removed if there's something higher", {
  expect_equal(strip_redundant("DR5 DR11 DRB1*11:XX"), "DRB1*11:XX")
  expect_equal(
    strip_redundant("DR1 DRB1*01:02 DR7 DRB1*07:XX"),
    c("DRB1*01:02 DRB1*07:XX")
  )
})

test_that("typing string works", {
  expect_equal(strip_redundant("A24(9) A10 A25 B70"), "A24 A25 B70")
  expect_equal(strip_redundant("A24(9) B64(14)"), "A24 B64")
})

test_that("NA returns NA", {
  expect_equal(strip_redundant(NA), NA_character_)
})

test_that("stripping is vectorized", {
  expect_equal(
    strip_redundant(c("A24(9)", "A24(9) A10 A25 B70", NA)),
    c("A24", "A24 A25 B70", NA)
  )
})

# add_xx_suffix -----------------------------------------------------------

test_that("serology is untouched", {
  expect_equal(add_xx_suffix("A1"), "A1")
  expect_equal(add_xx_suffix("A203"), "A203")
})

test_that("intermediate/high is untouched", {
  expect_equal(add_xx_suffix("A*01:AABJE"), "A*01:AABJE")
  expect_equal(add_xx_suffix("A*02:07"), "A*02:07")
})

test_that("XX codes are untouched", {
  expect_equal(add_xx_suffix("A*01:XX"), "A*01:XX")
})

test_that("low-res moleculars become XX", {
  expect_equal(add_xx_suffix("A*01:XX"), "A*01:XX")
})

test_that("NAs are handled", {
  expect_equal(add_xx_suffix(NA), NA)
  expect_equal(add_xx_suffix(""), "")
})

test_that("vectors work", {
  expect_equal(add_xx_suffix(c("A*01", NA, "A1")), c("A*01:XX", NA, "A1"))
})

# prefix_ambiguity() ------------------------------------------------

test_that("alleles with no ambiguities are untouched", {
  expect_equal(prefix_ambiguity("Cw3"), "Cw3")
  expect_equal(prefix_ambiguity("DQB1*04:XX"), "DQB1*04:XX")
  expect_equal(prefix_ambiguity("B*44:02:01:02S"), "B*44:02:01:02S")
})

test_that("complete ambiguities are untouched", {
  expect_equal(prefix_ambiguity("A*01:01/A*01:02"), "A*01:01/A*01:02")
  expect_equal(
    prefix_ambiguity("A*01:01/23/A*36:04"),
    "A*01:01/A*01:23/A*36:04"
  )
})

# TODO: what to do when ambiguity is in allele group (e.g. DRB1*11/13)

test_that("missing loci are filled in", {
  expect_equal(
    prefix_ambiguity("DRB1*01:01/01:02"),
    "DRB1*01:01/DRB1*01:02"
  )
})

test_that("missing allele groups and loci are filled in", {
  expect_equal(
    prefix_ambiguity("DRB1*01:01/02"),
    "DRB1*01:01/DRB1*01:02"
  )
})

test_that("NA returns NA", {
  expect_equal(prefix_ambiguity(NA), NA_character_)
})

test_that("vectors work", {
  expect_equal(
    prefix_ambiguity(c("A1", NA, "DRB1*01:01/02")),
    c("A1", NA, "DRB1*01:01/DRB1*01:02")
  )
})

# add_leading_zero() ------------------------------------------------------

test_that("single digits in serology are untouched", {
  expect_equal(add_leading_zero("A1"), "A1")
  expect_equal(add_leading_zero("Cw3"), "Cw3")
})

test_that("single digits in loci are untouched", {
  expect_equal(add_leading_zero("DPB1*01:01"), "DPB1*01:01")
  expect_equal(add_leading_zero("DQB1*04:XX"), "DQB1*04:XX")
})

test_that("leading zeros are added", {
  expect_equal(add_leading_zero("DRB1*1:1"), "DRB1*01:01")
  expect_equal(add_leading_zero("DQB1*8:XX"), "DQB1*08:XX")
  expect_equal(add_leading_zero("C*03:01/2/3"), "C*03:01/02/03")
  expect_equal(add_leading_zero("DP-1"), "DP-01")
})

test_that("NA returns NA", {
  expect_equal(add_leading_zero(NA), NA_character_)
})

test_that("vectors work", {
  expect_equal(
    add_leading_zero(c("A1", NA, "DRB1*1:1")),
    c("A1", NA, "DRB1*01:01")
  )
})


# remove_punctuation() ----------------------------------------------------

test_that("locus- and field separators are untouched", {
  expect_equal(remove_punctuation("C*03:01"), "C*03:01")
})

test_that("ambiguities are untouched", {
  expect_equal(remove_punctuation("A24/A25"), "A24/A25")
})

test_that("dashes are untouched", {
  expect_equal(remove_punctuation("DP-0201"), "DP-0201")
  expect_equal(remove_punctuation("HLA-A1"), "HLA-A1")
})

test_that("common punctuation is removed", {
  expect_equal(remove_punctuation("()[]{}A,;'1"), "A1")
})

test_that("common symbols are removed", {
  expect_equal(remove_punctuation("@!$%&A,;'1"), "A1")
})

test_that("GL string syntax is preserved", {
  expect_equal(
    remove_punctuation("hla#3.25.0#HLA-A*01:02+HLA-A*24:02:01:01"),
    "hla#3.25.0#HLA-A*01:02+HLA-A*24:02:01:01"
  )
  expect_equal(
    remove_punctuation("hla#3.29.0#HLA-DRB1*03:01:02~HLA-DRB5*01:01:01"),
    "hla#3.29.0#HLA-DRB1*03:01:02~HLA-DRB5*01:01:01"
  )
  expect_equal(
    remove_punctuation("HLA-A*02:69+HLA-A*23:30|HLA-A*02:302+HLA-A*23:26"),
    "HLA-A*02:69+HLA-A*23:30|HLA-A*02:302+HLA-A*23:26"
  )
  expect_equal(
    remove_punctuation("kir#2.12#KIR2DL5A*0010101+KIR2DL5A*0010201?"),
    "kir#2.12#KIR2DL5A*0010101+KIR2DL5A*0010201?"
  )
  expect_equal(
    remove_punctuation("HLA-DPB1*04:ANKZX+HLA-DPB1*04:FNVS^HLA-DRB3*XXXX"),
    "HLA-DPB1*04:ANKZX+HLA-DPB1*04:FNVS^HLA-DRB3*XXXX"
  )
})

test_that("NA returns NA", {
  expect_equal(remove_punctuation(NA), NA_character_)
})

test_that("vectors work", {
  expect_equal(
    remove_punctuation(c("A1", NA, "DPA1*,01$:01;'")),
    c("A1", NA, "DPA1*01:01")
  )
})

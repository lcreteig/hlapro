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

test_that("v2 alleles are untouched", {
  expect_equal(add_xx_suffix("Cw*0203"), "Cw*0203")
  expect_equal(add_xx_suffix("A*01MTN"), "A*01MTN")
})

test_that("low-res moleculars become XX", {
  expect_equal(add_xx_suffix("A*01"), "A*01:XX")
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

test_that("prefixing works for suffixes", {
  expect_equal(
    prefix_ambiguity("A*02:01/01L/04/07/09/15N"),
    "A*02:01/A*02:01L/A*02:04/A*02:07/A*02:09/A*02:15N"
  )
})

test_that("prefixing works for v2 alleles as well", {
  expect_equal(
    prefix_ambiguity("DRB4*0101/03/06"),
    "DRB4*0101/DRB4*0103/DRB4*0106"
  )
  expect_equal(
    prefix_ambiguity("A*0101/0103"),
    "A*0101/A*0103"
  )
  # alleles with suffixes
  expect_equal(
    prefix_ambiguity("A*0201/01L/04/07/09/15N"),
    "A*0201/A*0201L/A*0204/A*0207/A*0209/A*0215N"
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


# convert_v2_to_v3() ------------------------------------------------------

test_that("non-v2 alleles are untouched", {
  expect_equal(convert_v2_to_v3("A1"), "A1")
  expect_equal(convert_v2_to_v3("A2403"), "A2403")
  expect_equal(convert_v2_to_v3("A*01:XX"), "A*01:XX")
  expect_equal(convert_v2_to_v3("A*01:01"), "A*01:01")
  expect_equal(convert_v2_to_v3("DQB1*05:01:16"), "DQB1*05:01:16")
  expect_equal(convert_v2_to_v3("DPB1*1082:01"), "DPB1*1082:01")
})

test_that("v2 exceptions are handled", {
  expect_equal(convert_v2_to_v3("A*0105N"), "A*01:04:01:01")
  expect_equal(convert_v2_to_v3("A*1150"), "A*11:50Q")
  expect_equal(convert_v2_to_v3("B*9526"), "B*15:126")
  expect_equal(convert_v2_to_v3("Cw*010201"), "C*01:02:01")
  expect_equal(convert_v2_to_v3("DPB1*0802"), "DPB1*106:01")
})

test_that("v2 heuristic conversion works", {
  expect_equal(convert_v2_to_v3("A*0101"), "A*01:01")
  expect_equal(convert_v2_to_v3("A*010102"), "A*01:01:02")
  expect_equal(convert_v2_to_v3("DRB1*14125"), "DRB1*14:125") # odd # digits
  expect_equal(convert_v2_to_v3("A*01010203"), "A*01:01:02:03")
  # suffixes work
  expect_equal(convert_v2_to_v3("B*39010102L"), "B*39:01:01:02L")
  expect_equal(convert_v2_to_v3("A*02113N"), "A*02:113N")
  # Cw is C
  expect_equal(convert_v2_to_v3("Cw*0202"), "C*02:02")
  expect_equal(convert_v2_to_v3("Cw*02:02"), "C*02:02")
})

test_that("MAC and XX codes work", {
  # exceptions
  expect_equal(convert_v2_to_v3("A*24CECH"), "A*24:CWEN")
  expect_equal(convert_v2_to_v3("DPB1*03WFR"), "DPB1*03:FNYE")
  # heuristics
  expect_equal(convert_v2_to_v3("A*01XX"), "A*01:XX")
  expect_equal(convert_v2_to_v3("A*01AB"), "A*01:AB")
  expect_equal(convert_v2_to_v3("Cw*01BJZ"), "C*01:BJZ")
  expect_equal(convert_v2_to_v3("B*08YETY"), "B*08:YETY")
})

test_that("ambiguities work", {
  expect_equal(
    convert_v2_to_v3("DRB4*0101/DRB4*0103/DRB4*0106"),
    "DRB4*01:01/DRB4*01:03/DRB4*01:06"
  )
  expect_equal(
    convert_v2_to_v3("B*3919/B*3920/B*3921"),
    "B*39:19/B*39:20/B*39:24:01:01"
  )
})

test_that("NAs work", {
  expect_equal(convert_v2_to_v3(NA), NA)
})

test_that("vectors work", {
  expect_equal(
    convert_v2_to_v3(c("A*9209", "A*01:01", "A*011201", "A*0101/A*0102", NA)),
    c("A*02:109", "A*01:01", "A*01:12:01", "A*01:01/A*01:02", NA_character_)
  )
})


# convert_deleted() -------------------------------------------------------

test_that("conversion works with known gotchas", {
  expect_equal(convert_deleted("A*0105N"), "A*01:04:01:01N") # suffix
  expect_equal(convert_deleted("A*02:100"), NA_character_) # unassigned
  expect_equal(convert_deleted("C*03:99:01"), "C*01:169:01") # 2nd allele
  expect_equal(convert_deleted("C*03:12"), "C*03:19") # 1st allele
  expect_equal(convert_deleted("DPB1*35:01:02"), "DPB1*621:01") # not with date
})

test_that("conversion works with ambiguities", {
  expect_equal(convert_deleted("A*33:37/A*33:38"), "A*33:37/A*33:44")
})

test_that("not deleted is not converted", {
  expect_equal(convert_deleted("A1"), "A1")
  expect_equal(convert_deleted("B*08:01"), "B*08:01")
  expect_equal(convert_deleted("C*01:BJZ"), "C*01:BJZ")
  expect_equal(convert_deleted("A*0101"), "A*0101")
})

test_that("NAs work", {
  expect_equal(convert_deleted(NA), NA)
})

test_that("vectors work", {
  expect_equal(
    convert_deleted(c("A1", "C*12:139", "B*44:246N", "DQB1*02:01:14", NA)),
    c("A1", "C*12:139Q", NA_character_, "DQB1*02:02:17", NA_character_)
  )
  expect_equal(
    convert_deleted(c("A*24:22/A*24:329/A*24:378", NA, "DPB1*0701")),
    c("A*24:22/A*24:329Q/A*24:378Q", NA_character_, NA_character_)
  )
})

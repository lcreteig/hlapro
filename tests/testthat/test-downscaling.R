# reduce_to_nth_field -----------------------------------------------------

# TODO: handle old nomenclature like A*6603

test_that("serology is returned as is", {
  expect_equal(reduce_to_nth_field("A1", 1), "A1")
})

test_that("allele is returned as is if > n fields", {
  expect_equal(reduce_to_nth_field("A*01", 2), "A*01")
  expect_equal(reduce_to_nth_field("A*01", 3), "A*01")
  expect_equal(reduce_to_nth_field("A*01", 4), "A*01")
})

test_that("1st field reduction works", {
  expect_equal(reduce_to_nth_field("A*01:01:01", 1), "A*01")
  expect_equal(reduce_to_nth_field("A*01:01:01:02", 1), "A*01")
  expect_equal(reduce_to_nth_field("A*01:01:01G", 1), "A*01")
  expect_equal(reduce_to_nth_field("A*01:01P", 1), "A*01")
})

test_that("2nd field reduction works", {
  expect_equal(reduce_to_nth_field("A*01:01:01", 2), "A*01:01")
  expect_equal(reduce_to_nth_field("A*01:01:01:02", 2), "A*01:01")
  expect_equal(reduce_to_nth_field("A*01:01:01G", 2), "A*01:01")
  expect_equal(reduce_to_nth_field("A*01:01P", 2), "A*01:01")
})

test_that("3rd field reduction works", {
  expect_equal(reduce_to_nth_field("A*01:01:01", 3), "A*01:01:01")
  expect_equal(reduce_to_nth_field("A*01:01:01:02", 3), "A*01:01:01")
  expect_equal(reduce_to_nth_field("A*01:01:01G", 3), "A*01:01:01")
})

test_that("reduction is vectorized", {
  allele_in <- c("A1", "A*01", "A*01:01:01", "A*01:01:01:02")
  allele_out <- c("A1", "A*01", "A*01", "A*01")
  expect_equal(reduce_to_nth_field(allele_in, 1), allele_out)
})

test_that("NAs are handled", {
  expect_equal(reduce_to_nth_field(NA, 1), NA)
  expect_equal(
    reduce_to_nth_field(c("A*01:01", NA), 1),
    c("A*01", NA_character_)
  )
})


# etrl_convert() ----------------------------------------------------------

test_that("serology is not changed", {
  expect_equal(etrl_convert("A9"), "A9")
  expect_equal(etrl_convert("Cw3"), "Cw3")
})

test_that("first field gets XX", {
  expect_equal(etrl_convert("A*01"), "A*01:XX")
  expect_equal(etrl_convert("C*03"), "C*03:XX")
})

test_that("2nd field alleles in ETRL are returned as is", {
  expect_equal(etrl_convert("A*01:01"), "A*01:01")
})

test_that("2nd field alleles not in ETRL become XX", {
  expect_equal(etrl_convert("A*01:84"), "A*01:XX")
})

# TODO: consider expanding MAC and throwing warning if it cannot be
# converted unambiguously
test_that("MACs are XX codes", {
  expect_equal(etrl_convert("A*01:AABJE"), "A*01:XX")
  expect_equal(etrl_convert("DRB1*03:DWJAA"), "DRB1*03:XX")
})

# TODO: consider throwing warning if cannot be converted unambiguously
test_that("ambiguities are XX codes", {
  expect_equal(etrl_convert("C*01:02/C*01:03/C*01:04"), "C*01:XX")
})

test_that("null/alternative expression is not reduced", {
  expect_equal(etrl_convert("A*01:27N"), "")
  expect_equal(etrl_convert("A*24:02:01:02L"), "")
  expect_equal(etrl_convert("B*44:02:01:02S"), "")
  expect_equal(etrl_convert("C*01:121Q"), "")
})

test_that("HLA-prefix is removed", {
  expect_equal(etrl_convert("HLA-A*01:01:01)"), "A*01:01")
})

test_that("conversion is vectorized", {
  allele_in <- c(
    "A1", "A*01", "A*01:01N", "A*01:01", "A*01:AABJE",
    "A*01:01:01", "A*01:08"
  )
  allele_out <- c(
    "A1", "A*01:XX", "", "A*01:01", "A*01:XX",
    "A*01:01", "A*01:XX"
  )
  expect_equal(etrl_convert(allele_in), allele_out)
})

test_that("NAs are handled", {
  expect_equal(etrl_convert(NA), NA)
  expect_equal(etrl_convert(c("A*01", NA)), c("A*01:XX", NA))
})

# etrl_lookup() -----------------------------------------------------------

test_that("dataframe is returned", {
  etrl_out <- tidyr::tibble(
    Allele = NA_character_,
    `ET MatchDeterminantSplit` = NA_character_,
    `ET MatchDeterminantBroad` = NA_character_,
    `Public` = NA_character_
  )
  expect_equal(etrl_lookup(""), etrl_out, ignore_attr = TRUE)
  expect_equal(etrl_lookup(NA), etrl_out, ignore_attr = TRUE)
})

test_that("output dimensions correct for vector input", {
  expect_equal(dim(etrl_lookup(c("A1", "A*01:XX", "B00", "A*01", ""))), c(5, 4))
})


# get_broad() -------------------------------------------------------------

test_that("splits return broad", {
  expect_equal(get_broad("A24"), "A9")
  expect_equal(get_broad("A25"), "A10")
})

test_that("broads return NA", {
  expect_equal(get_broad("A1"), NA_character_)
  expect_equal(get_broad("A2"), NA_character_)
})

test_that("other inputs return NA", {
  expect_equal(get_broad("A*01:01"), NA_character_)
})

test_that("broad lookup is vectorized", {
  expect_equal(
    get_broad(c("A24", "A25", "A*01:01", NA)),
    c("A9", "A10", NA, NA)
  )
})

# get_public --------------------------------------------------------------

test_that("splits return public", {
  expect_equal(get_public("B64"), "Bw6")
  expect_equal(get_public("B77"), "Bw4")
  expect_equal(get_public("A23"), NA_character_)
  expect_equal(get_public("A24"), NA_character_)
})

test_that("broads return public", {
  expect_equal(get_public("B7"), "Bw6")
  expect_equal(get_public("B13"), "Bw4")
  expect_equal(get_public("A1"), NA_character_)
  expect_equal(get_public("A2"), NA_character_)
})

test_that("other inputs return NA", {
  expect_equal(get_public("A*01:01"), NA_character_)
  expect_equal(get_public("B*14:XX"), NA_character_)
})

test_that("public lookup is vectorized", {
  expect_equal(
    get_public(c("B64", "B13", "A*01:01", NA)),
    c("Bw6", "Bw4", NA, NA)
  )
})

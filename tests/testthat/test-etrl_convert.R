
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

# TODO: vectorize and test?

test_that("serology is not changed", {
  expect_equal(etrl_convert("A9"), "A9")
  expect_equal(etrl_convert("Cw3"), "Cw3")
})

test_that("first field gets XX", {
  # TODO: handle old nomenclature like A*6603
  expect_equal(etrl_convert("A*01"), "A*01:XX")
  expect_equal(etrl_convert("A*01"), "A*01:XX")
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

test_that("reduction works", {
  expect_equal(etrl_convert("A*01:01:01"), "A*01:01")
  expect_equal(etrl_convert("A*01:01:01:02"), "A*01:01")
  expect_equal(etrl_convert("A*01:28"), "A*01:XX")
  expect_equal(etrl_convert("A*01:230"), "A*01:XX")
})

test_that("null/alternative expression is not reduced", {
  expect_equal(etrl_convert("A*01:27N"), "")
  expect_equal(etrl_convert("A*24:02:01:02L"), "")
  expect_equal(etrl_convert("B*44:02:01:02S"), "")
  expect_equal(etrl_convert("C*01:121Q"), "")
})

test_that("groups work", {
  expect_equal(etrl_convert("A*01:01:01G"), "A*01:01")
  expect_equal(etrl_convert("A*01:01P"), "A*01:01")
  expect_equal(etrl_convert("A*01:09P"), "A*01:XX")
})

test_that("HLA-prefix is removed", {
  expect_equal(etrl_convert("HLA-A*01:01:01)"), "A*01:01")
})

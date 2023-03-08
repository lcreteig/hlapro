test_that("identical typings return empty string", {
  expect_equal(get_mismatches("A1", "A1"), "")
  expect_equal(get_mismatches("A1 A2", "A1 A2"), "")
  expect_equal(get_mismatches("A1 B8 C7 DR3 DQ1", "A1 B8 C7 DR3 DQ1"), "")
})

test_that("recipient HLAs not in donor typing are not mismatches", {
  expect_equal(get_mismatches("A1", "A1 B7"), "")
  expect_equal(get_mismatches("A1 B7", "A1 B7 B8"), "")
})

test_that("donor HLAs not in recipient typing are mismatches", {
  expect_equal(get_mismatches("DR4 DR7", "DR2 DR6"), "DR4 DR7")
  expect_equal(get_mismatches("A1 B7", "A1"), "B7")
  expect_equal(get_mismatches("A1 B7 B8", "A1 B7"), "B8")
  expect_equal(
    get_mismatches("A2 B7 C4 DR4 DQ4", "A1 B8 C7 DR3 DQ1"),
    "A2 B7 C4 DR4 DQ4"
  )
})

test_that("a single string is returned", {
  expect_length(get_mismatches("A2 B10 C7 DR3 DQ1", "A1 B8 C7 DR3 DQ1"), 1)
})

test_that("empty (donor) typing(s) equal empty string", {
  expect_equal(get_mismatches("", ""), "")
  expect_equal(get_mismatches("", "A1"), "")
})

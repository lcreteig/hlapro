test_that("basic positivity implementation works", {
  specificities <- c(
    "A*01:01", "A*02:01",
    "B*07:01", "B*08:02",
    "C*03:01", "C*04:01",
    "DQB1*01:01", "DQB1*02:01"
  )
  mfi <- c(15093, 394, 1, 499, 1000, 2000, 27, 1200)

  expect_equal(
    determine_positive_beads(specificities, mfi),
    c(
      "positive", "negative",
      "negative", "negative",
      "negative", "negative",
      "negative", "positive"
    )
  )
})

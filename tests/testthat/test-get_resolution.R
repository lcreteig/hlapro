test_that("serologicals are low", {
  expect_equal(get_resolution("A2"), "low")
  expect_equal(get_resolution("B70"), "low")
  expect_equal(get_resolution("Cw2"), "low")
  expect_equal(get_resolution("DQA-01"), "low")
  expect_equal(get_resolution("DP-0201"), "low")
  expect_equal(get_resolution("DP0201"), "low")
  expect_equal(get_resolution("HLA-DRB1*13"), "low")
})

test_that("v3 MACs are intermediate", {
  expect_equal(get_resolution("A*01:AABJE"), "intermediate")
  expect_equal(get_resolution("DRB1*07:GC"), "intermediate")
  # TODO: what about XX codes?
})

test_that("v2 MACs are intermediate", {
  expect_equal(get_resolution("B*15CFRG"), "intermediate")
  expect_equal(get_resolution("A*01KG"), "intermediate")
  expect_equal(get_resolution("DPB1*04BDVU"), "intermediate")
})

test_that("Ambiguous alleles are intermediate", {
  expect_equal(get_resolution("C*01:02/C*01:03/C*01:04"), "intermediate")
  expect_equal(get_resolution("HLA-A*23:26/HLA-A*23:39"), "intermediate")
  expect_equal(get_resolution("DRB1*13:02/13:36/13:67/13:96"), "intermediate")
  # TODO: what about ambiguities at the sub-protein level?
})

test_that(">2 field codes is high", {
  expect_equal(get_resolution("B*42:08"), "high")
  expect_equal(get_resolution("A*02:101:01"), "high")
  expect_equal(get_resolution("A*01:101:01:02"), "high")
  expect_equal(get_resolution("A*01:101:01:02N"), "high")
  expect_equal(get_resolution("HLA-DRB1*13:01:02"), "high")
  expect_equal(get_resolution("HLA-DRB1*13:01:01:02"), "high")
  expect_equal(get_resolution("HLA-A*24:09N"), "high")
  expect_equal(get_resolution("HLA-A*30:14L"), "high")
  expect_equal(get_resolution("HLA-A*24:02:01:02L"), "high")
  expect_equal(get_resolution("HLA-B*44:02:01:02S"), "high")
  expect_equal(get_resolution("HLA-A*32:11Q"), "high")
  # TODO: What about G groups, P groups?
  # TODO: return field code with high?
})

# get_resolution() --------------------------------------------------------

test_that("serologicals are low", {
  expect_equal(get_resolution("A2"), "low")
  expect_equal(get_resolution("A32(19)"), "low")
  expect_equal(get_resolution("B70"), "low")
  expect_equal(get_resolution("Cw2"), "low")
  expect_equal(get_resolution("DQA-01"), "low")
  expect_equal(get_resolution("DP-0201"), "low")
  expect_equal(get_resolution("DP0201"), "low")
  expect_equal(get_resolution("HLA-DRB1*13"), "low")
  expect_equal(get_resolution("A*01:XX"), "low")
})

test_that("MACs are intermediate", {
  # v3
  expect_equal(get_resolution("A*01:AABJE"), "intermediate")
  expect_equal(get_resolution("DRB1*07:GC"), "intermediate")
  expect_equal(get_resolution("C*17:XXHN"), "intermediate") # starts with XX
  # v2
  expect_equal(get_resolution("B*15CFRG"), "intermediate")
  expect_equal(get_resolution("A*01KG"), "intermediate")
  expect_equal(get_resolution("DPB1*04BDVU"), "intermediate")
})

test_that("Ambiguous alleles are handled", {
  # v3
  expect_equal(get_resolution("C*01:02/C*01:03/C*01:04"), "intermediate")
  expect_equal(get_resolution("C*01:02:01/C*01:02:02"), "high")
  expect_equal(get_resolution("HLA-A*23:26/HLA-A*23:39"), "intermediate")
  expect_equal(get_resolution("DRB1*13:02/13:36/13:67/13:96"), "intermediate")
  # v2
  expect_equal(get_resolution("DRB4*0101/DRB4*0103/DRB4*0106"), "intermediate")
})

test_that(">2 field codes is high", {
  expect_equal(get_resolution("B*42:08"), "high")
  expect_equal(get_resolution("DPB1*296:01"), "high")
  expect_equal(get_resolution("A*02:101:01"), "high")
  expect_equal(get_resolution("A*01:101:01:02"), "high")
  expect_equal(get_resolution("A*01:101:01:02N"), "high")
  expect_equal(get_resolution("HLA-DRB1*13:01:02"), "high")
  expect_equal(get_resolution("HLA-DRB1*13:01:01:02"), "high")
  # Suffixes
  expect_equal(get_resolution("HLA-A*24:09N"), "high")
  expect_equal(get_resolution("HLA-A*30:14L"), "high")
  expect_equal(get_resolution("HLA-A*24:02:01:02L"), "high")
  expect_equal(get_resolution("HLA-B*44:02:01:02S"), "high")
  expect_equal(get_resolution("HLA-A*32:11Q"), "high")
  # G/P groups
  expect_equal(get_resolution("A*01:01P"), "high")
  expect_equal(get_resolution("DQA1*02:01:01G"), "high")
  # v2
  expect_equal(get_resolution("A*0101"), "high")
})

test_that("extended mode works", {
  expect_equal(
    get_resolution("B*42:08", extended = TRUE),
    "high - second field"
  )
  expect_equal(
    get_resolution("A*9289", extended = TRUE),
    "high - second field"
  )
  expect_equal(
    get_resolution("A*02:101:01", extended = TRUE),
    "high - third field"
  )
  expect_equal(
    get_resolution("B*150303", extended = TRUE),
    "high - third field"
  )
  expect_equal(
    get_resolution("A*01:101:01:02N", extended = TRUE),
    "high - fourth field"
  )
  expect_equal(
    get_resolution("DRB1*15030101", extended = TRUE),
    "high - fourth field"
  )
  expect_equal(
    get_resolution("C*01:02/C*01:03/C*01:04", extended = TRUE),
    "intermediate - CAL"
  )
  expect_equal(
    get_resolution("A*01:AABJE", extended = TRUE),
    "intermediate - MAC"
  )
  expect_equal(get_resolution("A2", extended = TRUE), "serology - broad")
  expect_equal(get_resolution("A32", extended = TRUE), "serology - split")
  expect_equal(get_resolution("A*02", extended = TRUE), "molecular - broad")
  expect_equal(get_resolution("A*32:XX", extended = TRUE), "molecular - split")
  expect_equal(get_resolution("A203", extended = TRUE), "serology - associated")
  expect_equal(
    get_resolution("DR1403", extended = TRUE),
    "serology - associated"
  )
  expect_equal(get_resolution(NA, extended = TRUE), NA_character_)
})

test_that("vectors work", {
  expect_equal(
    get_resolution(c("A2", "A*01:AABJE", "B*42:08")),
    c("low", "intermediate", "high")
  )
})

test_that("mutate in dataframe works", {
  df <- tidyr::tibble(alleles = c("A2", "A*01:AABJE", "B*42:08"))
  df <- dplyr::mutate(df, allele_res = get_resolution(alleles))
  expect_equal(df$allele_res, c("low", "intermediate", "high"))
})

# get_n_fields() ----------------------------------------------------------

test_that("fields are counted accurately", {
  expect_equal(get_n_fields("A1"), 1)
  expect_equal(get_n_fields("A*01"), 1)
  expect_equal(get_n_fields("A*01:01"), 2)
  expect_equal(get_n_fields("A*01:01:01"), 3)
  expect_equal(get_n_fields("A*01:01:01:01"), 4)
})

test_that("count stops at ambiguous alleles", {
  expect_equal(get_n_fields("A*01:01/A*01:02"), 2)
})

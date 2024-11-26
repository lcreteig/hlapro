# reduce_to_nth_field -----------------------------------------------------

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

# get_serology() ----------------------------------------------------------

test_that("broads are returned as is", {
  expect_equal(get_serology("A1"), "A1")
  expect_equal(get_serology("B15"), "B15")
})

test_that("splits are returned as is", {
  expect_equal(get_serology("A24"), "A24")
  expect_equal(get_serology("B63"), "B63")
})

test_that("modern nomenclature with splits return splits", {
  expect_equal(get_serology("A*23:01"), "A23")
  expect_equal(get_serology("A*23:XX"), "A23")
  expect_equal(get_serology("A*23:01:01:11"), "A23")
  expect_equal(get_serology("B*14:01"), "B64")
})

test_that("modern nomenclature with only broad returns broad", {
  expect_equal(get_serology("B*07:02"), "B7")
  expect_equal(get_serology("B*14:XX"), "B14")
  expect_equal(get_serology("C*03:16"), "Cw16")
  expect_equal(get_serology("A*01:01:01:50"), "A1")
})

test_that("non-existing inputs return NA", {
  expect_equal(get_serology("A20"), NA_character_)
  expect_equal(get_serology("A*20"), NA_character_)
})

test_that("NA returns NA", {
  expect_equal(get_serology(NA), NA)
})

test_that("serology lookup is vectorized", {
  expect_equal(
    get_serology(c("A24", "B*14:01", NA)),
    c("A24", "B64", NA)
  )
})


# get_broad() -------------------------------------------------------------

test_that("splits return broads", {
  expect_equal(get_broad("A24"), "A9")
  expect_equal(get_broad("A25"), "A10")
})

test_that("broads return broads", {
  expect_equal(get_broad("A1"), "A1")
  expect_equal(get_broad("A2"), "A2")
})

test_that("modern nomenclature returns broads", {
  expect_equal(get_broad("A*01"), "A1")
  expect_equal(get_broad("A*01:XX"), "A1")
  expect_equal(get_broad("A*01:01:01"), "A1")
  expect_equal(get_broad("A*24"), "A9")
  expect_equal(get_broad("A*24:XX"), "A9")
  expect_equal(get_broad("A*24:02:01:102"), "A9")
})

test_that("non-existing inputs return NA", {
  expect_equal(get_broad("A20"), NA_character_)
  expect_equal(get_broad("A*20"), NA_character_)
})

test_that("NA returns NA", {
  expect_equal(get_broad(NA), NA)
})

test_that("broad lookup is vectorized", {
  expect_equal(
    get_broad(c("A24", "A25", "A*01:01", NA)),
    c("A9", "A10", "A1", NA)
  )
})

# get_split() -------------------------------------------------------------

test_that("splits return splits", {
  expect_equal(get_split("A24"), "A24")
  expect_equal(get_split("A*24"), "A24")
  expect_equal(get_split("A*24:02:01:102"), "A24")
  expect_equal(get_split("A25"), "A25")
})

test_that("broads return NA", {
  expect_equal(get_split("A1"), NA_character_)
  expect_equal(get_split("A9"), NA_character_)
})

test_that("allele nomenclature returns splits", {
  expect_equal(get_split("B*14:01"), "B64")
  expect_equal(get_split("B*14:02"), "B65")
  expect_equal(get_split("B*14:03"), NA_character_)
  expect_equal(get_split("B*14"), NA_character_)
})

test_that("non-existing inputs return NA", {
  expect_equal(get_split("A20"), NA_character_)
  expect_equal(get_split("A*20"), NA_character_)
})

test_that("NA returns NA", {
  expect_equal(get_split(NA), NA_character_)
})

test_that("split lookup is vectorized", {
  expect_equal(
    get_split(c("A24", "B*14:01", NA)),
    c("A24", "B64", NA)
  )
})

# get_public --------------------------------------------------------------

test_that("splits return public", {
  expect_equal(get_public("B64"), "Bw6")
  expect_equal(get_public("B77"), "Bw4")
  # exceptions!
  expect_equal(get_public("B62"), NA_character_)
  expect_equal(get_public("B71"), NA_character_)
})

test_that("broads return public", {
  expect_equal(get_public("B7"), "Bw6")
  expect_equal(get_public("B13"), "Bw4")
  expect_equal(get_public("A1"), NA_character_)
  expect_equal(get_public("A2"), NA_character_)
  # exceptions!
  expect_equal(get_public("B15"), NA_character_)
  expect_equal(get_public("B27"), NA_character_)
  expect_equal(get_public("B47"), NA_character_)
})

test_that("allele nomenclature works", {
  expect_equal(get_public("B*07:02"), "Bw6")
  expect_equal(get_public("B*13:01"), "Bw4")
  expect_equal(get_public("B*14:XX"), "Bw6")
  expect_equal(get_public("B*15:XX"), NA_character_)
  expect_equal(get_public("B*40:02:01:01"), "Bw6")
})

test_that("other loci return NA", {
  expect_equal(get_public("A23"), NA_character_)
  expect_equal(get_public("A24"), NA_character_)
  expect_equal(get_public("A*24:XX"), NA_character_)
  expect_equal(get_public("C*07:626"), NA_character_)
})

test_that("NA returns NA", {
  expect_equal(get_public(NA), NA)
})

test_that("public lookup is vectorized", {
  expect_equal(
    get_public(c("B64", "B13", "A*01:01", "B*07:09", NA)),
    c("Bw6", "Bw4", NA, "Bw6", NA)
  )
})

# reorder_alleles() -------------------------------------------------------

# TODO: add test for single alleles that don't match?

test_that("nothing is done when already in order", {
  in_order <- c("A1", "A2")
  to_order <- c("A1", "A2")
  expect_equal(reorder_alleles(in_order, to_order), to_order)
})

test_that("reordering with two unique alleles works", {
  in_order <- c("A2", "A1")
  to_order <- c("A1", "A2")
  expect_equal(reorder_alleles(in_order, to_order), c("A2", "A1"))
})

test_that("missings in to-be-ordered vector are respected", {
  in_order <- c("A1", "A2")
  to_order <- c("A1", NA)
  expect_equal(reorder_alleles(in_order, to_order), c("A1", NA))
  in_order <- c("A2", "A1")
  to_order <- c("A1", NA)
  expect_equal(reorder_alleles(in_order, to_order), c(NA, "A1"))
  in_order <- c("A2", "A1")
  to_order <- c(NA, NA)
  expect_equal(reorder_alleles(in_order, to_order), c(NA, NA))
})

test_that("missings in in-order vector are respected", {
  in_order <- c(NA, NA)
  to_order <- c("A1", "A2")
  expect_equal(reorder_alleles(in_order, to_order), c("A1", "A2"))
  in_order <- c("A1", NA)
  to_order <- c("A1", "A2")
  expect_equal(reorder_alleles(in_order, to_order), c("A1", "A2"))
  in_order <- c("A2", NA)
  to_order <- c("A1", "A2")
  expect_equal(reorder_alleles(in_order, to_order), c("A2", "A1"))
  in_order <- c("DQ3", NA)
  to_order <- c("DQ2", "DQ7") # one split
  expect_equal(reorder_alleles(in_order, to_order), c("DQ7", "DQ2"))
})

test_that("alleles with no serological equivalent are handled", {
  in_order <- c("C*04:09N", "C*01:02")
  to_order <- c("Cw1", NA)
  expect_equal(reorder_alleles(in_order, to_order), c(NA, "Cw1"))
  in_order <- c("Cw1", NA)
  to_order <- c("C*01:02", "C*04:09N")
  expect_equal(reorder_alleles(in_order, to_order), c("C*01:02", "C*04:09N"))
  in_order <- c("Cw1", NA)
  to_order <- c("C*04:09N", "C*01:02")
  expect_equal(reorder_alleles(in_order, to_order), c("C*01:02", "C*04:09N"))
  in_order <- c("C*04:09N", "C*01:02")
  to_order <- c("C*01:02", "C*04:09N")
  expect_equal(reorder_alleles(in_order, to_order), c("C*04:09N", "C*01:02"))
})

test_that("reordering works when resolutions differ", {
  in_order <- c("A1", "A2")
  to_order <- c("A*01:01", "A*02:01:01")
  expect_equal(reorder_alleles(in_order, to_order), c("A*01:01", "A*02:01:01"))
  in_order <- c("A*01:01", "A*02:01:01")
  to_order <- c("A1", "A2")
  expect_equal(reorder_alleles(in_order, to_order), c("A1", "A2"))
  in_order <- c("A1", "A2")
  to_order <- c("A*02:01:01", "A*01:01")
  expect_equal(reorder_alleles(in_order, to_order), c("A*01:01", "A*02:01:01"))
  in_order <- c("A*02:01:01", "A*01:01")
  to_order <- c("A1", "A2")
  expect_equal(reorder_alleles(in_order, to_order), c("A2", "A1"))
  # one split one broad
  in_order <- c("A24", "A2") # has a split
  to_order <- c("A9", "A2") # has a broad
  expect_equal(reorder_alleles(in_order, to_order), c("A9", "A2"))
  in_order <- c("A24", "A2") # has a split
  to_order <- c("A2", "A9") # has a broad
  expect_equal(reorder_alleles(in_order, to_order), c("A9", "A2"))
  # two splits with different order: reverse
  in_order <- c("A23", "A24")
  to_order <- c("A*24:02", "A*23:01")
  expect_equal(reorder_alleles(in_order, to_order), c("A*23:01", "A*24:02"))
  # two splits with same order: don't reverse
  in_order <- c("A23", "A24")
  to_order <- c("A*23:01", "A*24:02")
  expect_equal(reorder_alleles(in_order, to_order), c("A*23:01", "A*24:02"))
  # already in order; don't reverse
  in_order <- c("A*32:01:01:01", "A*26:01:01:01")
  to_order <- c("A25", "A26")
  expect_equal(reorder_alleles(in_order, to_order), c("A25", "A26"))
  # out of order; reverse
  in_order <- c("A*32:01:01:01", "A*26:01:01:01")
  to_order <- c("A26", "A25")
  expect_equal(reorder_alleles(in_order, to_order), c("A25", "A26"))
})

test_that("data frames with groupings work", {
  df_in <- tibble::tibble(
    id = c("001", "001", "002", "002"),
    locus = c("A", "A", "A", "A"),
    allele = c("1", "2", "1", "2"),
    in_order = c("A2", "A1", "A3", NA),
    to_order = c("A*01:01", "A*02:01:01", "A3", NA)
  )
  df_out <- df_in |>
    dplyr::mutate(to_order = c("A*02:01:01", "A*01:01", "A3", NA))
  expect_equal(
    df_in |>
      dplyr::group_by(id, locus) |>
      dplyr::mutate(
        to_order =
          reorder_alleles(in_order, to_order)
      ) |>
      dplyr::ungroup(),
    df_out
  )
})

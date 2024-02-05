# get_positive_eplets -----------------------------------------------------

mock_eplet_df <- dplyr::tribble(
  ~name, ~alleles,
  "D111", "D*01:01",
  "D112", "D*01:01",
  "D113", "D*01:01",
  "D113", "D*01:02"
)

test_that("negative eplets are 'subtracted' from positive eplets", {
  luminex_df <- dplyr::tribble(
    ~sampleID, ~allele, ~positive,
    "001", "D*01:01", TRUE,
    "001", "D*01:02", FALSE
  )

  expect_equal(
    get_positive_eplets(luminex_df, sampleID, allele, positive, mock_eplet_df),
    dplyr::tibble(
      sampleID = c("001", "001"),
      eplets_pos = c("D111", "D112")
    )
  )
})

test_that("output differs per ID", {
  luminex_df <- dplyr::tribble(
    ~sampleID, ~allele, ~positive,
    "001", "D*01:01", TRUE,
    "001", "D*01:02", FALSE,
    "002", "D*01:01", TRUE,
    "002", "D*01:02", TRUE,
  )

  expect_equal(
    get_positive_eplets(luminex_df, sampleID, allele, positive, mock_eplet_df),
    dplyr::tibble(
      sampleID = c("001", "001", "002", "002", "002"),
      eplets_pos = c("D111", "D112", "D111", "D112", "D113")
    )
  )
})

test_that("all rows are returned if no overlapping negative eplets", {
  luminex_df <- dplyr::tribble(
    ~sampleID, ~allele, ~positive,
    "001", "D*01:01", TRUE,
    "001", "D*01:02", TRUE,
    "002", "D*01:02", TRUE
  )

  expect_equal(
    get_positive_eplets(luminex_df, sampleID, allele, positive, mock_eplet_df),
    dplyr::tibble(
      sampleID = c("001", "001", "001", "002"),
      eplets_pos = c("D111", "D112", "D113", "D113")
    )
  )
})

test_that("empty rows are returned if no eplets are positive", {
  luminex_df <- dplyr::tribble(
    ~sampleID, ~allele, ~positive,
    "001", "D*01:01", FALSE,
    "001", "D*01:02", FALSE
  )

  expect_equal(
    get_positive_eplets(luminex_df, sampleID, allele, positive, mock_eplet_df),
    dplyr::tibble(
      sampleID = character(),
      eplets_pos = character()
    )
  )

  luminex_df <- luminex_df |>
    dplyr::add_row(sampleID = "002", allele = "D*01:02", positive = TRUE)

  expect_equal(
    get_positive_eplets(luminex_df, sampleID, allele, positive, mock_eplet_df),
    dplyr::tibble(sampleID = c("002"), eplets_pos = c("D113"))
  )
})

# eplet_registry ----------------------------------------------------------

df_eplets <- load_eplet_registry()

# lookup_alleles() --------------------------------------------------------

test_that("right alleles are returned", {
  expect_equal(lookup_alleles(df_eplets, "23L"), list(`23L` = "DQB1*04:01"))
  expect_equal(
    lookup_alleles(df_eplets, "4Q", allele_set = "all"),
    list(`4Q` = c(
      "DRB1*01:144", "DRB1*07:01", "DRB1*07:03", "DRB1*07:04", "DRB1*07:07",
      "DRB1*07:09", "DRB1*07:121", "DRB1*07:136", "DRB1*07:139", "DRB1*07:151",
      "DRB1*07:152", "DRB1*07:153", "DRB1*09:01", "DRB1*09:20", "DRB1*09:21",
      "DRB1*09:31", "DRB1*09:32", "DRB4*01:01", "DRB4*01:02", "DRB4*01:03",
      "DRB4*01:07", "DRB4*01:151", "DRB4*01:152", "DRB4*01:155", "DRB4*01:156",
      "DRB4*01:168", "DRB4*01:173", "DRB4*01:176"
    ))
  )
})

test_that("allele lookup is vectorized", {
  expect_equal(
    lookup_alleles(df_eplets, c("9D", "9T")),
    list(
      `9D` = c(
        "B*08:01", "C*06:02", "C*07:01", "C*07:02",
        "C*07:04", "C*18:01", "C*18:02"
      ),
      `9T` = c(
        "A*29:01", "A*29:02", "A*31:01", "A*33:01", "A*33:03"
      )
    )
  )
})


test_that("NAs are dealt with", {
  expect_equal(
    lookup_alleles(df_eplets, c("23L", NA)),
    set_names(list("DQB1*04:01", character(0)), c("23L", NA))
  )
})

# lookup_eplets() ---------------------------------------------------------

test_that("right eplets are returned", {
  expect_equal(
    lookup_eplets(df_eplets, "DPA1*03:01"),
    list(`DPA1*03:01` = c(
      "11M", "28E", "31M", "50Q", "56A", "65I", "66S", "127L", "160F", "190T"
    ))
  )
})

test_that("allele lookup is vectorized", {
  expect_equal(
    lookup_eplets(df_eplets, c("DPA1*03:01", "DQA1*01:04")),
    list(
      `DPA1*03:01` = c(
        "11M", "28E", "31M", "50Q", "56A", "65I", "66S", "127L", "160F", "190T"
      ),
      `DQA1*01:04` = c(
        "2G", "25YT", "40E", "52SK", "75I", "129QS", "160A", "160AD", "185I"
      )
    )
  )
})

test_that("allele lookup is vectorized", {
  expect_equal(
    lookup_eplets(df_eplets, c("DPA1*03:01", NA_character_)),
    set_names(
      list(
        c(
          "11M", "28E", "31M", "50Q", "56A", "65I",
          "66S", "127L", "160F", "190T"
        ),
        NA_character_
      ),
      c("DPA1*03:01", NA)
    )
  )
})

# load_eplet_registry() ---------------------------------------------------

test_that("returns path when return_path = TRUE", {
  path <- load_eplet_registry(return_path = TRUE)
  expect_type(path, "character")
})

test_that("load_eplet_registry prints message when print_version = TRUE", {
  expect_message(load_eplet_registry(print_version = TRUE))
})

test_that("table has 9 columns", {
  expect_equal(length(df_eplets), 9)
})

test_that("column names and types are correct", {
  eplet_registry_info <- c(
    id = "character",
    name = "character",
    residue_type = "character",
    description = "character",
    exposition = "character",
    confirmation = "character",
    locus_group = "character",
    source = "character",
    alleles = "character"
  )
  expect_equal(purrr::map_chr(df_eplets, class), eplet_registry_info)
})

test_that("low cardinality character columns contain expected values", {
  # exposition
  values_exposition <- df_eplets |>
    dplyr::filter(locus_group != "DRDQDP") |> # these are undefined
    dplyr::pull(exposition) |>
    unique()
  expect_setequal(
    values_exposition,
    c("Very Low", "Low", "Intermediate", "High", NA)
  )

  # residue type
  expect_setequal(
    unique(df_eplets$residue_type),
    c("eplet", "reactivity pattern")
  )

  # database
  expect_setequal(
    unique(df_eplets$locus_group),
    c("ABC", "DRB", "DQ", "DP", "DRDQDP")
  )

  # source
  expect_setequal(
    unique(df_eplets$source),
    c("luminex", "all")
  )
})

test_that("table has no empty eplets/alleles", {
  expect_true(sum(is.na(df_eplets$name) | df_eplets$name == "") == 0)
  expect_true(sum(is.na(df_eplets$alleles) | df_eplets$alleles == "") == 0)
})

test_that("table has no duplicate eplets", {
  eplet_count <- df_eplets |>
    dplyr::distinct(name, description, exposition, confirmation, locus_group) |>
    dplyr::add_count(name) |>
    dplyr::pull(n)

  expect_true(all(eplet_count == 1))
})

test_that("a few randomly selected cells have same value as on the website", {
  expect_equal(
    dplyr::pull(df_eplets[df_eplets$name == "3P", ], "residue_type")[1],
    "eplet"
  )

  expect_equal(
    dplyr::pull(df_eplets[df_eplets$name == "77N+85VG", ], "residue_type")[1],
    "reactivity pattern"
  )

  expect_equal(
    dplyr::pull(df_eplets[df_eplets$name == "37Y", ], "exposition")[1],
    "High"
  )
  expect_equal(
    dplyr::pull(df_eplets[df_eplets$name == "71SA", ], "confirmation")[1],
    "Yes"
  )
  expect_equal(
    dplyr::pull(df_eplets[df_eplets$name == "69AA+65QI", ], "confirmation")[1],
    "N/A"
  )
  expect_equal(
    dplyr::pull(df_eplets[df_eplets$name == "45EV", ], "description")[1],
    "45E46V47Y"
  )
  expect_equal(
    dplyr::pull(df_eplets[df_eplets$name == "11AV", ], "description")[1],
    "11A12V"
  )
  expect_equal(
    dplyr::pull(
      dplyr::filter(df_eplets, name == "9T", source == "luminex"),
      "alleles"
    ),
    c("A*29:01", "A*29:02", "A*31:01", "A*33:01", "A*33:03")
  )
  expect_equal(
    dplyr::pull(
      dplyr::filter(df_eplets, name == "3P", source == "all"),
      "alleles"
    ),
    c(
      "DQB1*06:01", "DQB1*06:103", "DQB1*06:205", "DQB1*06:243",
      "DQB1*06:359", "DQB1*06:382", "DQB1*06:415", "DQB1*06:472"
    )
  )
})

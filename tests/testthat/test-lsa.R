test_that("mock lot files are parsed correctly", {
  # Class I
  df_sa1 <- tibble::tibble(
    LotID = "3012456 3065421-SA1",
    ExpirationDate = "12/12/2024",
    AssayName = "3012456 3065421-SA1",
    AssayVersion = 1,
    LogicName = "SA1-000",
    LogicDescription = "LIFECODES Single Antigen Class I",
    LocusID = 10,
    Extension = "SA1",
    IsDNA = 0,
    Con1Upper = 100,
    Con1Lower = 180,
    PCUpper = 15000,
    PCLower = 13000,
    MFIThreshold = 500,
    Phenotypes = NA,
    Bead = c(1, 2, 3, 4, 5, 5, 6),
    Antigens = c(NA, NA, "A*01:01", "A*02:01", "B*14:01", "Bw6", "C*07:01"),
    Cutoff = c(NA, NA, 3.75, 4.01, 3.79, 3.79, 3.49),
    RAD = c(NA, 4, 1.192, 1.210, 1.304, 1.304, 0.790),
    BackgroundMFI = c(NA, 1500, 180, 170, 180, 180, 245),
    LRA = c(NA, 3, 1, 1, 2, 2, 3),
    Consensus = c("N1", "P1", NA, NA, NA, NA, NA),
    Serology = c(NA, NA, "A1", "A2", "B64(14)", "B64(14)", "Cw7"),
    antigen_id = c("NC", "PC", "103", "104", "105", "105", "106")
  )

  expect_equal(
    read_lotfile(test_path("luminex", "3012456 3065421-SA1.eds")),
    df_sa1)

  # Class II
  df_sa2 <- tibble::tibble(
    LotID = "3062145 3045216-SA2",
    ExpirationDate = "07/25/2027",
    AssayName = "3062145 3045216-SA2",
    AssayVersion = 1,
    LogicName = "SA2-000",
    LogicDescription = "LIFECODES Single Antigen Class II",
    LocusID = 15,
    Extension = "SA2",
    IsDNA = 0,
    Con1Upper = 370,
    Con1Lower = 160,
    PCUpper = 16500,
    PCLower = 16220,
    MFIThreshold = 650,
    Phenotypes = NA,
    Bead = c(1, 2, 3, 4, 5, 5, 6, 6, 7, 7),
    Antigens = c(NA, NA, "DRB1*01:01", "DRB3*02:02", "DQA1*02:01", "DQB1*02:01",
                 "DQA1*01:02", "DQB1*06:04", "DPA1*01:03", "DPB1*01:01"),
    Cutoff = c(NA, NA, 3.90, 4.20, 3.95, 3.95, 5.80, 5.80, 3.90, 3.90),
    RAD = c(NA, 4, 1.006, 0.850, 0.920, 0.920, 0.875, 0.875, 0.930, 0.930),
    BackgroundMFI = c(NA, 1500, 135, 140, 145, 145, 290, 290, 150, 150),
    LRA = c(NA, 3, 1, 1, 3, 3, 3, 3, 2, 2),
    Consensus = c("N1", "P1", NA, NA, NA, NA, NA, NA, NA, NA),
    Serology = c(NA, NA, "DR1", "DR52", "DQ2", "DQ2", "DQ6(1)", "DQ6(1)",
                 "DPw1", "DPw1"),
    antigen_id = c("NC", "PC", "203", "204", "205", "205", "206", "206",
                   "207", "207")
  )

  expect_equal(
    read_lotfile(test_path("luminex", "3062145 3045216-SA2.eds")),
    df_sa2)
})

test_that("mock csv files are parsed correctly", {
  # These rds files are the result of applying read_lum_csv to a simplified
  # version of a Luminex csv file. They were created with a version of the code
  # that produced results that were verified by hand to correspond to the result
  # from Immucor's MATCH IT!® software.
  # Class I
  df_sa1 <- readRDS(test_path("luminex", "LSA1-test.rds"))
  expect_equal(
    read_lum_csv(csv_filepath = test_path("luminex", "LSA1-test.csv"),
                 lots_path = test_path("luminex")),
    df_sa1)
  # Class II
  df_sa2 <- readRDS(test_path("luminex", "LSA2-test.rds"))
  expect_equal(
    read_lum_csv(csv_filepath = test_path("luminex", "LSA2-test.csv"),
                 lots_path = test_path("luminex")),
    df_sa2)
})

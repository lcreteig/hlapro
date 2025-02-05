dl_and_del <- function(url) {
  temp <- tempfile()
  utils::download.file(url, temp) # download xls and save as tempfile
  df_xls <- readxl::read_xls(temp)
  unlink(temp) # delete the file
  return(df_xls)
}

q_title <- paste(
  "Do you want to download some files from GitHub and the NMDP website",
  "that contain the modern (v3) nomenclature equivalents of v2 HLA alleles?"
)

if (dl_permission(q_title) == 1) {
  rlang::check_installed("readr", reason = "to read the data")

  url_nmdp <- "https://bioinformatics.bethematchclinical.org/WorkArea/"

  ### Download files

  # Obsolete allele-specific codes. N.B. unfortunately no longer online
  macs_obsolete <- dl_and_del(paste0(url_nmdp, "DownloadAsset.aspx?id=6551")) |>
    dplyr::select(v2 = `Version 2 type`, v3 = `Version 3 type`)
  # Obsolete DPB1-specific codes N.B. unfortunately no longer online
  macs_dpb1 <- dl_and_del(paste0(url_nmdp, "DownloadAsset.aspx?id=6550")) |>
    dplyr::select(v2 = `Version 2 type`, v3 = `Version 3 type`)
  # v2 to v3 exceptions
  nomencl <- readr::read_table(paste0(url_ipd, "Nomenclature_2009.txt"),
    skip = 2, col_names = c("v2", "v3")
  )

  ### Join

  v2_to_v3 <- dplyr::left_join(nomencl, deleted_changed,
    by = dplyr::join_by(v2 == allele_old)
  ) |>
    # if no v3 assigned, look it up in deleted alleles list
    dplyr::mutate(v3 = dplyr::if_else(v3 == "None", allele_new, v3)) |>
    dplyr::select(!c(allele_new, date_changed)) |>
    # Cw*0422 is corrected to Cw*0421 (still v2), which is "C*04:15:02" in v3
    dplyr::mutate(v3 = dplyr::if_else(v3 == "Cw*0421", "C*04:15:02", v3)) |>
    dplyr::bind_rows(macs_obsolete, macs_dpb1)

  rlang::check_installed("usethis", reason = "to save the data")
  usethis::use_data(v2_to_v3, overwrite = TRUE) # save dataframe
} else {
  stop("Download aborted")
}

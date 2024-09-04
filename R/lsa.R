read_lum_csv <- function(csv_filepath, lots_path) {
  # read luminex csv line-by-line
  dat_lines <- readLines(csv_filepath)

  # get lot number
  lot_id <- dat_lines |>
    stringr::str_subset("^\"ProtocolName\",\"") |>
    stringr::str_extract(r"(\d+\s\d+\-\w+)")

  # read corresponding lotfile and extract its data
  df_eds <- read_lotfile(file.path(lots_path, paste0(lot_id, ".eds")))

  # get run date
  batch_date <- dat_lines |>
    stringr::str_subset("^\"Date\",\"") |>
    stringr::str_extract(r"(\d{2}/\d{2}/\d{4})") |>
    as.Date(format = "%m/%d/%Y")

  # get all indices where a data rectangle starts (excluding logs)
  table_starts <- stringr::str_which(
    dat_lines,
    "\"DataType:\",\"(?!Audit Logs|Warnings/Errors)"
  )
  pivot_cols <- c("Location", "Sample", "Analyte:", "Dilution Factor")
  dfs <- vector("list", length(table_starts) - 1) # initialize vec to store dfs

  # for each table
  for (ii in seq_along(utils::head(table_starts, -1))) {
    line_start <- table_starts[ii] + 1
    line_end <- table_starts[ii + 1] - 2
    # type of data
    measure <- stringr::str_remove_all(
      dat_lines[line_start - 1],
      r"([\\:,""]|DataType)"
    )

    # convert to dataframe
    df <- utils::read.csv(
      text = dat_lines[line_start:line_end],
      check.names = FALSE
    )

    # if there's a value for each bead
    if (!all(colnames(df) %in% pivot_cols)) {
      df <- df |>
        # redundant row throws off pivot
        dplyr::filter(dplyr::if_any(
          dplyr::any_of("Analyte:"), ~ . != "BeadID:"
        )) |>
        # reshape to one row per bead
        tidyr::pivot_longer(
          cols = -dplyr::any_of(pivot_cols),
          names_to = "antigen_id",
          values_to = measure
        ) |>
        dplyr::select(-dplyr::any_of("Analyte:")) |> # drop redundant columns
        dplyr::filter(.data$antigen_id != "Total Events") # delete sum row
    }
    dfs[[ii]] <- df
  }

  df_lsa <- purrr::reduce(dfs, dplyr::full_join) |> # concat all dataframes
    # add single values we got earlier
    tibble::add_column(
      batch_date = batch_date,
      lot_id = lot_id,
      .before = "antigen_id"
    ) |>
    dplyr::left_join(df_eds) |> # add lot info
    dplyr::group_by(.data$Sample, .data$LRA) |>
    # add computed colummns
    dplyr::mutate(
      mfi_lra = .data$Median / min(.data$Median),
      assignment = dplyr::if_else(
        (.data$mfi_lra > .data$Cutoff) & (.data$Median > .data$MFIThreshold),
        "Positive", "Negative"
      ),
      bg_adjusted = .data$Median - .data$BackgroundMFI,
      ad_mfi = .data$Median / .data$RAD,
      ad_bg_adjusted = .data$bg_adjusted / .data$RAD
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(row = dplyr::row_number(), .by = c("Sample", "antigen_id")) |>
    # format according to matchit view
    tidyr::pivot_wider(
      names_from = "row", values_from = "Antigens",
      names_prefix = "Antigens"
    ) |> # separate 1st, 2nd row per bead
    tidyr::pivot_wider(
      names_from = "LRA", values_from = c("Antigens1", "Serology"),
      names_prefix = "locus"
    ) |> # separate loci
    dplyr::arrange(dplyr::desc(.data$mfi_lra))

  # Format depending on class I or II assays
  if (stringr::str_ends(lot_id, "-SA1")) {
    df_lsa |>
      dplyr::select("Sample", "Count", "antigen_id", "Cutoff", "Median",
        "mfi_lra", "assignment", "bg_adjusted", "ad_mfi", "ad_bg_adjusted",
        A = "Antigens1_locus1", B = "Antigens1_locus2", C = "Antigens1_locus3",
        Bw = "Antigens2",
        A_serology = "Serology_locus1",
        B_serology = "Serology_locus2",
        C_serology = "Serology_locus3",
        "RAD"
      )
  } else if (stringr::str_ends(lot_id, "-SA2")) {
    df_lsa |>
      dplyr::mutate(
        DPB = dplyr::if_else(stringr::str_starts(.data$Antigens2, "DPB"),
          .data$Antigens2, NA
        ),
        DQB = dplyr::if_else(stringr::str_starts(.data$Antigens2, "DQB"),
          .data$Antigens2, NA
        )
      ) |>
      dplyr::select("Sample", "Count", "antigen_id", "Cutoff",
        "Median", "mfi_lra", "assignment",
        "bg_adjusted", "ad_mfi", "ad_bg_adjusted",
        DR_DR5X = "Antigens1_locus1", DQA = "Antigens1_locus3", "DQB",
        DPA = "Antigens1_locus2", "DPB", DR_serology = "Serology_locus1",
        DQ_serology = "Serology_locus3", DP_serology = "Serology_locus2", "RAD"
      )
  }
}

read_lotfile <- function(filepath) {
  xml_doc <- xml2::read_xml(filepath) |>
    xml2::xml_ns_strip() # strip default namespace to find nodes more quickly

  # Get all singular fields (e.g. lot ID)
  single_nodes <- xml_doc |>
    xml2::xml_children() |>
    xml2::xml_path() |>
    stringr::str_subset("AntigenIDs", negate = TRUE) # deal with beads later

  # Extract their text into a table
  lot_info <- single_nodes |>
    # loop, add to list
    purrr::map(\(x) xml2::xml_text(xml2::xml_find_all(xml_doc, x))) |>
    # clean yup and restore names
    purrr::set_names(stringr::str_remove_all(single_nodes, "/.*/")) |>
    # convert to data frame
    tibble::as_tibble_row() |>
    utils::type.convert()

  antigen_ids <- xml2::xml_find_all(xml_doc, "//AntigenID")

  # Get bead/antigen IDs
  ids <- antigen_ids |>
    purrr::map(xml2::xml_attrs) |>
    purrr::list_c() |>
    tibble::tibble("antigen_id" = _)

  # Get lists of bead-related data, and convert to dataframe
  antigens <- antigen_ids |>
    xml2::as_list() |>
    tibble::as_tibble_col() |> # single column, with list for each bead
    tidyr::unnest_wider("value") |> # expand into several columns
    dplyr::bind_cols(ids) |> # add antigen id column
    tidyr::unnest(cols = !"antigen_id") |>
    tidyr::unnest(cols = !"antigen_id") |> #  (double) unlist all columns
    tidyr::unnest(cols = "Antigens", keep_empty = TRUE) |> # unlist HLAs
    utils::type.convert(as.is = TRUE)

  # combine singular fields with bead data
  invisible(dplyr::bind_cols(lot_info, antigens))
}

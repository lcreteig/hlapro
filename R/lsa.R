read_lum_csv <- function(csv_filepath, lots_path) {
  # read luminex csv line-by-line
  dat_lines <- readLines(csv_filepath)

  # get lot number
  lot_id <- dat_lines |>
    stringr::str_subset("^\"ProtocolName\",\"") |>
    stringr::str_extract(r"(\d+\s\d+\-\w+)")

  # read corresponding lotfile and extract its data
  df_eds <- read_lotfile(file.path(lots_path, lot_id, ".eds"))
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

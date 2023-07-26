# TODO: extract_alleles_df() and gl_to_df() have a different API;
# consider equalizing (i.e., take in a df or a vector)

#' Split an HLA typing string into alleles
#'
#' @description Takes in a space-separated HLA typing string and splits it into
#'   its constituent loci and alleles ("A_1", "A_2", "DRB1_1").
#'
#'   `extract_alleles_str()` takes in a single string, and returns a named
#'   character vector of alleles.
#'
#'   `extract_alleles_df()` takes in a data frame, where one column contains the
#'   typing string, and returns the same data frame along with a new column for
#'   each allele.
#'
#' @param string String, space-separated HLA typing.
#' @param df A data frame.
#' @param col_typing The column in `df` that contains a space-separated HLA
#'   typing `string` for each row.
#' @param loci A string or character vector with the loci you are interested in.
#'   Only these alleles will be returned. Defaults to all. `DRB.` is used for
#'   DRB3, DRB4, and DRB5.
#' @param strip_locus Include the locus in the output or remove it?
#'   - If `TRUE` (the default), the locus will be removed from the extracted
#'     alleles.
#'   - If `FALSE`, will retain the locus as it was in the original typing.
#'
#' @return Either a character vector or a data frame with the named alleles.
#'  A warning will be shown if any loci in the input have more than two alleles.
#' @export
#'
#' @examples
#' typing <- "A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53"
#' extract_alleles_str(typing, loci = "A")
#' extract_alleles_str(typing)
#'
#' df <- tidyr::tibble(typing = typing)
#' extract_alleles_df(df, typing, loci = c("A", "B", "C"))
#'
#' # Can also handle newer nomenclature
#' extract_alleles_str("DQB1*03:01 DQB1*05:01 DRB1*04:AMR",
#'   loci = c("DRB1", "DQB1")
#' )
extract_alleles_str <- function(string,
                                loci = c(
                                  "A", "B", "C",
                                  "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB."
                                ),
                                strip_locus = TRUE) {
  loci <- rlang::arg_match(loci, multiple = TRUE)
  stopifnot(is.logical(strip_locus))

  # get rid of any leading/trailing/double spaces
  string <- stringr::str_squish(string)

  if (any(count_alleles(string, loci) > 2)) {
    rlang::warn(c("One or more loci found with more than 2 alleles.",
      "x" = "`extract_alleles_str()` will only pick the first two.",
      "i" = "Use `hlapro::count_alleles()` to find out more."
    ))
  }

  extract_str <- function(locus, string) {
    pattern <- build_pattern(locus, strip_locus)
    allele_names <- stringr::str_c(locus, c("_1", "_2"))

    # return named list (locus_allele) with matches
    rlang::set_names(
      stringr::str_match(string, pattern)[, 2:3],
      allele_names
    )
  }

  purrr::map(loci, \(x) extract_str(x, string)) |>
    purrr::flatten_chr()
}

#' @rdname extract_alleles_str
#' @export
extract_alleles_df <- function(df,
                               col_typing,
                               loci = c(
                                 "A", "B", "C",
                                 "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB."
                               ),
                               strip_locus = TRUE) {
  loci <- rlang::arg_match(loci, multiple = TRUE)

  # get rid of any leading/trailing/double spaces
  df <- dplyr::mutate(df, dplyr::across({{ col_typing }}, stringr::str_squish))

  warn_flag <- df |>
    dplyr::pull({{ col_typing }}) |>
    append(NA) |> # force this to be a list
    count_alleles(loci) |>
    purrr::list_c() |>
    purrr::some(\(x) !is.na(x) & x > 2)

  if (warn_flag) {
    rlang::warn(c("One or more loci found with more than 2 alleles.",
      "x" = "`extract_alleles_df()` will only pick the first two.",
      "i" = "Use `hlapro::count_alleles()` to find out more."
    ))
  }

  extract_df <- function(locus, df, col_typing) {
    pattern <- build_pattern(locus, strip_locus)
    allele_names <- stringr::str_c(locus, c("_1", "_2"))

    tidyr::extract(df, {{ col_typing }},
      into = allele_names,
      regex = pattern,
      remove = FALSE
    )
  }

  # Split typing for each locus into two columns; combine resulting data frames
  purrr::map(loci, \(x) extract_df(x, df, {{ col_typing }})) |>
    purrr::reduce(dplyr::full_join)
}

#' Count the number of alleles for each locus in an HLA typing string
#'
#' @description
#'
#' `count_alleles()` takes in a character vector or string of HLA typings, and
#' returns the number of alleles that the string contains for each locus.
#'
#' This can be useful when validating whether a typing string contains a typing
#' for each locus. Also, this function can alert you when a typing string
#' contains more than two alleles for each locus, which can be intentional
#' (e.g. when the typing contains both the broad and the split) or a mistake.
#'
#' @inheritParams extract_alleles_str
#' @param typings A(n) (character vector of) HLA typing string(s).
#'
#' @return A (list of) named integer vector(s), with the loci as names, and
#'   the number of found alleles per locus as values.
#' @export
#'
#' @examples
#' typing <- "A1 A2 B7 B8 Cw3 DQ5 DQ8 DR4 DR11 DR52 DR53"
#' count_alleles(typing, loci = "A")
#' count_alleles(typing)
#'
#' # Also works with character vectors
#' typing <- c("A1 A2 B7 B8 Cw3", "Cw3 Cw7", NA)
#' count_alleles(typing, loci = c("A", "B", "C"))
count_alleles <- function(typings,
                          loci = c(
                            "A", "B", "C",
                            "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB."
                          )) {
  loci <- rlang::arg_match(loci, multiple = TRUE)

  # locus should not be preceded by / (an ambiguity) or a capital,
  # and should be followed either by a digit (serology) or * (molecular)
  base_pattern <- r"((?<!\/|[A-Z]){locus}(\d+|\*))"

  # count alleles per locus and return as named integer vector
  count_alleles_1 <- function(string, loci) {
    counts <- stringr::str_count(string,
      pattern = stringr::str_glue(
        base_pattern,
        locus = locus_patterns[loci]
      )
    )
    rlang::set_names(counts, names(locus_patterns[loci]))
  }

  if (length(typings) > 1) {
    return(purrr::map(typings, \(x) count_alleles_1(x, loci)))
  }
  count_alleles_1(typings, loci)
}

build_pattern <- function(locus, strip_locus) {
  if (strip_locus) { # do not include the locus and * in the match
    allele_pattern <- r"((?:{locus}[-\*]?(\S+)))"
  } else { # match locus the locus and following non-spaces
    allele_pattern <- r"(({locus}\S+))" # match locus and following non-spaces
  }

  regexps <- list(
    # don't match if locus is preceded by : or other capital letter. This
    # prevents matching in NMDP Multiple Allele Codes (":AABJE") or other loci
    # ("B" in "DRB1") or other prefixes (e.g. "A" in "HLA-")
    neg = "(?<![:A-Z])",
    allele = allele_pattern, # match locus and following non-spaces
    # don't capture other alleles in between: i.e. any that aren't current locus
    inbetween = r"((?:\s(?!{locus})\S+)*\s?)"
  )

  # build the regular expression from its sub-parts
  stringr::str_c(
    regexps$neg,
    stringr::str_glue(regexps$allele, locus = locus_patterns[locus]),
    stringr::str_glue(regexps$inbetween, locus = locus_patterns[locus]),
    stringr::str_glue(regexps$allele, locus = locus_patterns[locus]),
    "?" # second allele is optional
  )
}

locus_patterns <- c(
  A = "A",
  B = "B(?![Ww])", # B cannot be followed by "W" (public)
  C = "C[Ww]?", # in serological notation, C is followed by "w"
  DPA1 = "DPA[-1]", # either "DQA-" or DQA1
  DPB1 = "DP(?:w|-|B1)", # either "DPw" "DP-" or "DPB1"
  DQA1 = "DQA[-1]", # either "DQA-" or DQA1
  DQB1 = "DQ(?!A)(?:B1)?", # could be "DQB1" but also "DQ" with no A after
  DRB1 = "DR(?!A|B[2-9]|5[1-3])(?:B1)?", # either "DRB1"
  # or "DR" (excluding DRA, DR51-53)
  # TODO: now what we have count_alleles(), remove *Neg
  DRB. = r"(DR(?:(?=5[1-3])|B(?=[3-5](?!\*Neg))))" # either DR51/52/53 or
  # DRB3/4/5. Sometimes all 3 are specified, with 1-2 having suffix "*Neg".
  # These should be skipped, in order to still retrieve full typing)
)

df_to_gl <- function(df,
                     namespace = "hla",
                     version_or_date = NULL,
                     col_typing = "glstring",
                     loci = c(
                       "A", "B", "C",
                       "DPA1", "DPB1", "DQA1", "DQB1", "DRB1", "DRB."
                     ),
                     suffixes = c("1", "2"),
                     sep = "_") {
  loci <- rlang::arg_match(loci, multiple = TRUE)

  # build names of columns containing HLAs (e.g. "A_1", "A_2", "B_1", etc.)
  typing_cols <- paste0(
    rep(loci, each = length(suffixes)),
    sep, suffixes
  )

  df |>
    tidyr::pivot_longer(
      cols = tidyr::any_of(typing_cols),
      names_to = c("locus", "allele"),
      names_sep = sep,
      values_to = col_typing
    ) |>
    dplyr::summarise(dplyr::across(col_typing, ~ vec_to_gl(.x,
      namespace = namespace,
      version_or_date = version_or_date
    )))
}

#' @importFrom rlang .data
gl_to_df <- function(glstrings) {
  glsc_lst <- purrr::map(glstrings, gl_to_vec)
  tidyr::tibble(glsc = glsc_lst) |>
    tidyr::unnest_wider(.data$glsc) |>
    dplyr::mutate(
      glstring_id = dplyr::row_number(),
      .before = dplyr::everything()
    ) |>
    tidyr::unnest_longer(.data$allele_list, values_to = "typing") |>
    dplyr::group_by(.data$glstring_id) |>
    dplyr::mutate(
      locus = get_loci(.data$typing),
      allele = dplyr::row_number(.data$locus)
    ) |>
    tidyr::pivot_wider(
      id_cols = .data$glstring_id,
      values_from = .data$typing,
      names_from = c(.data$locus, .data$allele),
      names_sep = "_",
      names_expand = TRUE
    ) |>
    dplyr::ungroup()
}

vec_to_gl <- function(allele_list, namespace = "hla", version_or_date = NULL) {
  if (is.null(version_or_date)) {
    version_or_date <- Sys.Date()
  }
  if (all(is.na(allele_list))) {
    return(NA)
  }

  allele_list <- allele_list[!is.na(allele_list)] |>
    stringr::str_sort() |>
    add_hla_prefix()
  loci <- get_loci(allele_list)

  for (ii in seq_along(loci)[-1]) { # for each allele
    sep <- "^" # assume they are different loci, use locus delimiter
    if (loci[ii] == loci[ii - 1]) { # but if locus is same as previous allele
      sep <- "+" # use genotype delimiter
    }
    allele_list[ii] <- paste0(sep, allele_list[ii]) # add delim to allele
  }

  # put all alleles and their delims together, prefix with metadata
  stringr::str_c(
    namespace, "#",
    version_or_date, "#",
    stringr::str_flatten(allele_list)
  )
}

gl_to_vec <- function(glstring) {
  if (is.na(glstring)) {
    return(NA)
  }

  glsc <- stringr::str_split_1(glstring, "#")

  list(
    namespace = glsc[1],
    version_or_date = glsc[2],
    allele_list = stringr::str_split_1(glsc[3], r"([\+\^"])")
  )
}

get_loci <- function(allele_list) {
  # get every letter/digit before a "*"
  stringr::str_extract(allele_list, "\\w+(?=\\*)")
}

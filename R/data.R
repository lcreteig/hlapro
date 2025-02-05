#' ETRL HLA Tables
#'
#' EuroTransplant Reference Laboratory tables with Human Leukocyte Antigen
#' alleles. These can be used as lookup tables to get the serological equivalent
#' of a given allele.
#'
#' @format ## `etrl_hla`
#' A data frame with 737 rows and 4 columns:
#' \describe{
#'   \item{Allele}{HLA allele (in modern allele notation)}
#'   \item{ET MatchDeterminantSplit}{serological equivalent at the split level}
#'   \item{ET MatchDeterminantBroad}{serological equivalent at the broad level}
#'   \item{Public}{Bw4 or Bw6 (if the allele has that epitope)}
#' }
#' @source <https://etrl.eurotransplant.org/resources/hla-tables/>
"etrl_hla"

#' v2 to v3 HLA nomenclature conversion table
#' @description
#'
#' In April 2010, the HLA nomenclature officially changed from "v2" to "v3",
#' which mandates the use of colons ("`:`") to delimit fields. For example:
#' `A*01010101` (v2) is since written as `A*01:01:01:01`. Note though that not
#' all conversions are this straightforward, there are a number of exceptions.
#'
#' This conversion table tries to collect all standard v2 --> v3 conversions
#' from different data sources (see below).
#' @format ## `v2_to_v3`
#' A data frame with 5,151 rows and 2 columns:
#' \describe{
#'   \item{v2}{v2 version of the allele}
#'   \item{v3}{v3 version of the allele}
#' }
#' @source
#' - Most pre-2010 alleles:
#' <https://github.com/ANHIG/IMGTHLA/blob/Latest/Nomenclature_2009.txt>
#' - Additional mappings for deleted alleles:
#' <https://github.com/ANHIG/IMGTHLA/blob/Latest/Deleted_alleles.txt>
#' - Allele-specific & DPB1-specific multiple-allele codes (see 4. and 5. under
#' "Allele Code Update to Version 3 WHO Nomenclature")
#' N.B. Unfortunately this link no longer works
#' <https://bioinformatics.bethematchclinical.org/hla-resources/allele-codes/>
"v2_to_v3"

#' Table of changed/deleted HLA alleles
#'
#' The IPD-IMGT/HLA database contains a list of alleles in the database that
#' were once deleted (because, for example, they were later discovered to be
#' identical to an existing allele) or changed (for example, their sequence was
#' renamed). This conversion table lists all the old alleles and their new
#' equivalents.
#'
#' The new names have been parsed from a description of the change contained in
#' the IPD-IMGT/HLA allele list. This makes it easier to programmatically
#' convert deleted/changed alleles without having to manually look up the new
#' designation.
#'
#' N.B. The deleted alleles list is updated from time to time with new releases
#' of the IPD-IMGT/HLA database (though relatively infrequently; definitely
#' not every release). It's possible the conversion table in this package is
#' based on an older version.
#'
#' @format ## `deleted_changed`
#' A data frame with 275 rows and 2 columns:
#' \describe{
#'   \item{allele_old}{allele name that was deleted/changed}
#'   \item{allele_new}{current allele name}
#' }
#' @source
#' - <https://hla.alleles.org/alleles/deleted.html>
#' - csv version of this data (used as the source)
#'   <https://github.com/ANHIG/IMGTHLA/blob/Latest/Nomenclature_2009.txt>
"deleted_changed"

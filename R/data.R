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

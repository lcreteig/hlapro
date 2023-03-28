validate_allele <- function(allele) {
  stringr::str_detect(allele, r"(^(?:HLA-)?(?<locus>A|Bw?|Cw?|DRA|DR(?:B[1-9])?|DQA1?|DQ(?:B1)?|DPA1?|DP(?:B1)?)[\*-]?(?<allele>\d{1,4}):?(?<protein>[A-Z]{2,5}$|\d{2,3}P?)?(?<coding>:\d{2,3}G?)?(?<noncoding>:\d{2,3})?(?<suffix>[NLSCAQ]?)(\/(?:\1\*)?\d{2,3}(?::\d{2,3}[NLSCAQ]?)*)*$)")
}

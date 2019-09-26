# utilities for the project 190108-mzt-rna-stability

#' Load the codon composition anf the decay rate data
#' function to process and load the codon composition and decay rates for fish
#'
#' @param cc_dp path to 190108-mzt-rna-stability/data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/zfish_codon_composition.csv
#' @param decay_dp path to 190108-mzt-rna-stability/results/19-01-11-GetDecayRateFromTimeCourse/results_data/estimated_decay_rates.csv
#'
#' @return
#' @export
#' @importFrom magrittr %>%
load_decay_aa_codon_composition_data <- function(cc_dp, decay_dp) {
  cc <- readr::read_csv(cc_dp) %>%
    dplyr::rename(Gene_ID = id_col)
  decay <- readr::read_csv(decay_dp)

  # apply filter alpha > 0 (drop genes low expression)
  decay <- decay %>%
    dplyr::select(.data$Gene_ID:.data$estimate) %>%
    tidyr::spread(key = .data$term, value = .data$estimate) %>%
    dplyr::filter(.data$alpha > 0) %>%
    dplyr::select(-.data$alpha) %>%
    dplyr::rename(decay_rate = beta) %>%
    dplyr::inner_join(cc)

  return(decay)
}

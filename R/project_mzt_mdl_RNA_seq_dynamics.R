# functions to model mRNA-seq dynamics ------------------------------------
# create model matrix

#' Make Model Frame
#'
#' @param log2FCtimecourse_dp character to file path: data/19-02-05-FoldChangeData/data/log2FC_earlyVSlate_tidytimecourse.csv
#' @param pls_optimality_dp character to file path: results/19-01-18-ObtainEstimateOfCodonOptimalityPLS/results_data/pls_components_fish_genes.csv
#' @param utrs_dp character to file path: data/19-01-17-Get-ORFS-UTRS-codon-composition/sequence-data/fish_seq_data_cds_3utr.csv
#' @param kmers_to_keep character vector of k-mers to add, this predictor
#' represents the number of times that a k-mer is present in the 3' UTR seqs
#' @param .sample_condition character any of aamanit2in_polya wt_polya wt_ribo
#' @param .maternal logical, True only matenal genes, assume there is column is_maternal
#' @param minimum_time_point numeric lowest time point used for modeling
#'
#' @return tibble model frame
#' @export
#' @importFrom dplyr %>%
#' @examples
make_mdl_frame <- function(log2FCtimecourse_dp,
                           pls_optimality_dp,
                           utrs_dp,
                           kmers_to_keep,
                           .sample_condition = "wt_ribo",
                           .maternal = TRUE,
                           minimum_time_point = -1) {

  purrr::map_lgl(
    .x = c(pls_optimality_dp, utrs_dp, log2FCtimecourse_dp),
    ~!file.exists(.)
  ) %>%
    any() %>%
    if (.) stop("supplied file(s) not found")

  # tidy the data
  opt <- readr::read_csv(pls_optimality_dp) %>%
    dplyr::select(-decay_rate)

  utrs <- readr::read_csv(utrs_dp) %>%
    dplyr::select(-.data$coding) %>%
    dplyr::rename(Gene_ID = ensembl_gene_id)

  fc_tc <- readr::read_csv(log2FCtimecourse_dp)

  # add the kmer counts as predictors
  for (kmer in kmers_to_keep) {
    utrs[, kmer] <- stringr::str_count(utrs$`3utr`, kmer)
  }

  utrs <- dplyr::select(utrs, -`3utr`)

  if (.maternal) fc_tc <- dplyr::filter(fc_tc, .data$is_maternal)
  fc_tc <-
    fc_tc %>%
    dplyr::filter(
      .data$time > minimum_time_point,
      .data$sample_condition == .sample_condition,
      !is.infinite(.data$log2FC) # some log2FC are infinite I will drop them
    ) %>%
    dplyr::select(-.data$is_maternal, -.data$sample_condition)

  mdl_frame <- dplyr::left_join(
    fc_tc,
    dplyr::inner_join(opt, utrs, by = "Gene_ID"),
    by = "Gene_ID"
  ) %>%
    .[stats::complete.cases(.), ]

  mdl_frame

}

#' Count k-mers in DNA sequence
#'
#' @param sequence str dna sequence
#' @param k integer k-mer length
#'
#' @return tibble with 2 columns: k-mer (string), n (number of sites)
#' @export
#' @importFrom dplyr %>%
#' @examples
#' count_kmers("AAGGTTCC", k = 4)
count_kmers <- function(sequence, k) {
  sequence <- stringr::str_to_upper(sequence)
  if (nchar(sequence) < k) {
    stop("sequence shorter than k")
  }
  seq(from = 1, to = nchar(sequence) - k + 1) %>%
    purrr::map_chr(~ substr(sequence, ., . + k - 1)) %>%
    table() %>%
    dplyr::as_tibble() %>%
    dplyr::rename_("kmer" = ".")
}

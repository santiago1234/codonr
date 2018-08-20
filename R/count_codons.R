#' Count codons in DNA sequence
#'
#' @param seq A character.
#'
#' @return The frequency of the codons in \code{seq}
#' @export
#' @importFrom magrittr %>%
#' @examples
#' count_codons('ACGGGG')
#'
count_codons <- function(seq) {
  if (nchar(seq) %% 3 != 0) {
    stop("sequence not a multiple of 3")
  }
  seq <- toupper(seq)
  seq.int(from = 1, to = nchar(seq) - 2, by = 3) %>%
    purrr::map_chr(function(x) substr(
      seq,
      x, x + 2
    )) %>%
    table() %>%
    tibble::as_tibble() %>%
    dplyr::rename_(codon = ".")
}

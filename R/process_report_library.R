#' find the position of the 1st stop codon that occurs in frame
#'
#' @param seq dna string
#'
#' @return int of position if not position found -1
find_first_stop_codon_position <- function(seq) {
  for (i in seq(from = 1, to = nchar(seq), by = 3)) {
    codon <- substr(seq, i, i + 2)
    if (codon %in% c("TAG", "TAA", "TGA")) {
      return(as.integer(i))
    } # stop codon found
  }

  # no premature stop codon
  return(as.integer(-1))
}

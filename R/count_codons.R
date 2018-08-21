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


#' Codon Composition of Sequences
#'
#' @param orfs tibble, data with cds sequences
#' @param orf_col character, column name containing cds
#' @param id_col character, column name with gene/name id
#' @param cores integer, number of course to use
#'
#' @return tible n by x, codon counts for each sequence, n is the number
#' of sequences and x is the number of codons (x=61) it may contain more thant
#' 61 in case there is a codon loke ACN
#' @export
#' @importFrom magrittr %>%
#' @examples
#' orf <- tibble::tibble(cds= c("ACGTTT", "TTTCCC"), id=c("s0","s1"))
#' codon_composition(orf, orf_col = "cds", id_col = "id")
codon_composition <- function(orfs, orf_col, id_col, cores=4) {

  # rename the columns for easy manipulation
  orfs <- dplyr::rename_(orfs, "cds_seq" = orf_col, "id_col" = id_col)
  n_ids <- length(unique(orfs$id_col))
  if (n_ids != nrow(orfs)){
    stop("id_col should contain unique identifiers")
  }
  composition_tb <-
    orfs %>%
    split(.$id_col) %>%
    parallel::mclapply(
      function(x) {
        count_codons(x$cds_seq) %>%
          dplyr::mutate(id_col = x$id_col)
      },
      mc.cores = cores
    ) %>%
    do.call(what = dplyr::bind_rows) %>%
    tidyr::spread(key = codon, value = n) %>%
    dplyr::mutate_if(is.numeric, dplyr::funs(replace(., is.na(.), 0)))

  return(composition_tb)
}

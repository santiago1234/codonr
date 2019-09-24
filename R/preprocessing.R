#' Count codons
#'
#' @param orf character DNA sequence
#'
#' @return tibble with codon counts 1x64
#'
#' @examples
#' count_codons("AAA")
count_codons2 <- function(orf) {
  if (nchar(orf) %% 3 != 0) {
    stop("sequence not a multiple of 3")
  }

  composition <-
    orf %>%
    Biostrings::DNAString() %>%
    Biostrings::oligonucleotideFrequency(width = 3, step = 3) %>%
    base::as.list() %>%
    dplyr::as_tibble()
  # the c_ prefix is used in the preprocessing step
  colnames(composition) <- paste0("c_", colnames(composition))

  composition
}


#' Initial Preprocessing function
#'
#' This function computes the codon composition and an addition variable for
#' the log length of the open reading frame
#' @param data: data to transform (i.e. train_set, test_set)
#'
#' @return A tibble: 6 x 71
#' @export
#'
#' @examples
initial_preproc_codon_composition <- function(data) {
  data %>%
    dplyr::mutate(
      len_log10_coding = log10(nchar(coding)),
      cc = purrr::map(coding, count_codons2)
    ) %>%
    tidyr::unnest() %>%
    dplyr::select(-coding)
}


#' Preprate train and test set for modeling
#'
#' @return list with train and test data
#' @export
#'
#' @examples
prepare_train_and_test_sets <- function() {
  list(
    train_set = initial_preproc_codon_composition(train_set),
    test_set = initial_preproc_codon_composition(test_set)
  )
}


#' Data preprocessing
#'
#' Preprocess the train_set data for trainin/testing  models
#' The following preprocessing steps are implemented
#' \enumerate{
#'   \item Impute the 3' UTR length with median
#'   \item Apply \href{https://www.ncbi.nlm.nih.gov/pubmed/16711760}{spatial sign} transformation to codon composition
#'   \item Normalize all numeric variables
#'   \item Dummy transformation to catergorical variables: `cell_type`, `specie`, and `datatype`
#' }
#'
#' @param train_set_prepared training_set, the train_set output of \code{\link{prepare_train_and_test_sets}}
#'
#' @return trained Data Recipe
#' @export
#'
#' @examples
#' # first use the function \code{\link{prepare_train_and_test_sets}} to add the
#' codon composition and 3' UTR length
#' data_prepared <- prepare_train_and_test_sets()
#' rcipe <- preprocessing(data_prepared$train_set)
#' # apply the pre-processing to test_set
#' recipes::bake(rcipe, data_prepared$test_set)
preprocessing <- function(train_set_prepared) {

  rcipe <-
    recipes::recipe(decay_rate ~ ., data = train_set_prepared) %>%
    recipes::update_role(gene_id, new_role = "id variable") %>%
    recipes::step_medianimpute(utrlenlog) %>%
    recipes::step_spatialsign(dplyr::starts_with("c_")) %>%
    recipes::step_normalize(recipes::all_numeric()) %>%
    recipes::step_dummy(specie, cell_type, datatype, one_hot = FALSE)

  # train the recipe with the training data
  recipes::prep(rcipe, training = train_set_prepared)

}


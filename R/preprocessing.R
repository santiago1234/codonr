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
  message("couting codons ... ")
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
#' @return list with train and test data (X and y)
#' @export
#'
#' @examples
prepare_train_and_test_sets <- function() {

  train_set = initial_preproc_codon_composition(train_set)
  test_set = initial_preproc_codon_composition(test_set)

  list(
    X_train = dplyr::select(train_set, -decay_rate),
    X_test = dplyr::select(test_set, -decay_rate),
    y_train = train_set$decay_rate,
    y_test = test_set$decay_rate
  )
}


#' Data preprocessing
#'
#' Preprocess the predictors matrix X (train_set) data for trainin/testing  models
#' The following preprocessing steps are implemented
#' \enumerate{
#'   \item Impute the 3' UTR length with median
#'   \item Apply \href{https://www.ncbi.nlm.nih.gov/pubmed/16711760}{spatial sign} transformation to codon composition
#'   \item Normalize all numeric variables
#'   \item Dummy transformation for catergorical variables: `cell_type`, `specie`, and `datatype`
#' }
#'
#' @param X_train training_set predictors, the X_train output of \code{\link{prepare_train_and_test_sets}}
#'
#' @return trained Data Recipe
#' @export
#'
#' @examples
#' # first use the function \code{\link{prepare_train_and_test_sets}} to add the
#' codon composition and 3' UTR length
#' data_prepared <- prepare_train_and_test_sets()
#' rcipe <- preprocessing(data_prepared$X_train)
#' # apply the pre-processing to test_set
#' recipes::bake(rcipe, data_prepared$X_test)
preprocessing <- function(X_train) {

  rcipe <-
    recipes::recipe(X_train) %>%
    recipes::step_rm(gene_id) %>%
    recipes::step_medianimpute(utrlenlog) %>%
    recipes::step_spatialsign(dplyr::starts_with("c_")) %>%
    recipes::step_dummy(specie, cell_type, datatype, one_hot = FALSE) %>%
    recipes::step_normalize(recipes::all_numeric())

  # train the recipe with the training data
  recipes::prep(rcipe, training = X_train)

}


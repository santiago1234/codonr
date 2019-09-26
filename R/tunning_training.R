#' Cross fold validation technique to train mRNA stability predictor
#'
#' This function generates the trainControl object to be passed to the
#' caret train function.
#' Here, I use grouped 10 fold cross validation, the motivation for this
#' is that the same gene (gene_id) cannot be in both groups, train and test,
#' at the same time. The same gene can appear multiple times, for example the human
#' data that contains profiles for different cell types and with differen techniques
#' The group variable is the `gene_id``
#' @param k integer number of folds
#'
#' @return trainControl object to use in \code{caret::train } function
#' @export
#'
#' @examples
#' cross_fold_validation()
cross_fold_validation <- function(k = 10) {
  grouped_folds_in_geneid <- caret::groupKFold(train_set$gene_id, k = k)
  caret::trainControl(method = "cv", index = grouped_folds_in_geneid)
}

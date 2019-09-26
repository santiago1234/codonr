#' Cross fold validation technique to train mRNA stability predictor
#'
#' This function generates the trainControl object to be passed to the
#' caret train function.
#' Here, I use grouped 10 fold cross validation, the motivation for this
#' is that the same gene (gene_id) cannot be in both groups, train and test,
#' at the same time. The same gene can appear multiple times in the data,
#' for example the human data that contains profiles for different cell types
#' and with differen techniques
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


#' Train Machine Learning model to predict mRNA stability
#'
#' @param ncores
#'
#' @return
#' @export
#'
#' @examples
tune_models <- function(ncores = 32) {
  # parallel processing for grid search
  cl <- parallel::makePSOCKcluster(ncores)
  doParallel::registerDoParallel(cl)

  message("group 10-fold cv")
  control_object <- cross_fold_validation()

  # LINEAR MODELS
  message("linear regression ...")
  set.seed(669)
  linear_reg <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "lm",
    trControl = control_object,
    metric = "Rsquared"
  )

  message("pls regression ...")
  set.seed(669)
  pls_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "pls",
    tuneLength = 15,
    trControl = control_object,
    metric = "Rsquared"
  )

  message("elastic net regression ...")
  set.seed(669)
  enet_grid <- expand.grid(
    lambda = c(0, .001, .01, .1),
    fraction = seq(0.05, 1, length.out = 20)
  )

  enet_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "pls",
    tuneLength = 15,
    trControl = control_object,
    metric = "Rsquared"
  )


  # collect results
  all_resamples <- caret::resamples(list(
    "linear regression" = linear_reg,
    "partial least squares" = pls_model,
    "elastic net" = enet_model
  ))

  parallel::stopCluster(cl)
  all_resamples
}

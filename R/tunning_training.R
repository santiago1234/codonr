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
cross_fold_validation <- function(k = 10) {
  grouped_folds_in_geneid <- caret::groupKFold(train_set$gene_id, k = k)
  caret::trainControl(method = "cv", index = grouped_folds_in_geneid)
}


#' Train Machine Learning models to predict mRNA stability
#'
#' Evaluates a set of machine learning models, also perfoms hyperparameter tunning.
#' I follow the same approach as in chapter 10 of Applied Predictive Modeling book.
#' This function uses the \code{datos_preprocessed} training data for tunning,
#' @param ncores integer number of cores for parallel processing. set to the
#' numbers of cores in your machine
#'
#' @return list of trained models
#' @export
#' @note This function is computational expensive and it requires parallel processing
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
    method = "enet",
    tuneGrid = enet_grid,
    trControl = control_object,
    metric = "Rsquared"
  )

  # NON LINEAR MODELS

  message("earth regression ...")
  set.seed(669)
  earth_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "earth",
    tuneGrid = expand.grid(degree = 1, nprune = 2:25),
    trControl = control_object,
    metric = "Rsquared"
  )


  message("svm regression ...")
  set.seed(669)
  svmr_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "svmLinear",
    tuneLength = 10,
    trControl = control_object,
    metric = "Rsquared"
  )

  message("neural network regression ...")
  set.seed(669)
  nnet_grid <- expand.grid(
    decay = c(.001, .01, .1),
    size = seq(1, 27, by = 2),
    bag = FALSE
  )
  nnet_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "avNNet",
    tuneGrid = nnet_grid,
    trControl = control_object,
    metric = "Rsquared"
  )

  message("k-nearest neighbors regression ...")
  set.seed(669)
  knn_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "knn",
    tuneGrid = base::data.frame(k = 1:20),
    trControl = control_object,
    metric = "Rsquared"
  )


  # REGRESSION AND MODEL TRESS

  message("rpart  regression ...")
  set.seed(669)
  rpart_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "rpart",
    tuneLength = 30,
    trControl = control_object,
    metric = "Rsquared"
  )

  message("ctree  regression ...")
  set.seed(669)
  ctree_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "ctree",
    tuneLength = 10,
    trControl = control_object,
    metric = "Rsquared"
  )

  message("mtree  regression ...")
  set.seed(669)
  mt_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "M5",
    trControl = control_object,
    metric = "Rsquared"
  )

  # ENSEMBLE METHODS

  message("random forest regression ...")
  set.seed(669)
  rf_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "rf",
    tuneLength = 10,
    ntrees = 1000,
    importance = TRUE,
    trControl = control_object,
    metric = "Rsquared"
  )

  message("gradient boosting regression ...")
  set.seed(669)
  gbm_grid <- expand.grid(
    interaction.depth = seq(1, 7, by = 2),
    n.trees = seq(100, 1000, by = 50),
    shrinkage = c(0.01, 0.1, 0.5),
    n.minobsinnode = 10
  )
  gbm_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "gbm",
    tuneGrid = gbm_grid,
    verbose = FALSE,
    trControl = control_object,
    metric = "Rsquared"
  )

  message("cubist regression ...")
  set.seed(669)
  cubist_grid <- expand.grid(
    commmitees = c(1, 5, 10, 50, 75, 100),
    neighbors = c(0, 1, 3, 5, 7, 9)
  )
  cubist_model <- caret::train(
    x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
    method = "cubist",
    tuneGrid = cubist_grid,
    trControl = control_object,
    metric = "Rsquared"
  )

  # collect results
  all_models <- list(
    "linear regression" = linear_reg,
    "partial least squares" = pls_model,
    "elastic net" = enet_model,
    "adaptative regression spline" = earth_model,
    "support vector machine" = svmr_model,
    "neural network" = nnet_model,
    "k-nearest neighbors" = knn_model,
    "single tree" = rpart_model,
    "model tree" = mt_model,
    "c tree" = ctree_model,
    "random forest" = rf_model,
    "boosted trees" = gbm_model,
    "cubist" = cubist_model
  )

  parallel::stopCluster(cl) # stop cluster
  all_models
}


train_final_model <- function(X_train = datos_preprocessed$X_train, y_train = datos_preprocessed$y_train) {
  NULL
}

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
#' @param output_dir_to_save_models character: output directory to save trained models
#'
#' @return None, all the model if not trained will be saved to \code{output_dir_to_save_models}
#' @export
#' @note This function is computational expensive and it requires parallel processing
tune_models <- function(ncores = 32, output_dir_to_save_models = "results-trained-models") {
  # helper functions to check if model exists and save them

  make_model_id_name <- function(mdl_name) {
    file.path(output_dir_to_save_models, paste0(mdl_name, ".rds"))
  }

  has_model_been_trained <- function(mdl_name) {
    # returns true if the model exits
    file.exists(make_model_id_name(mdl_name))
  }

  save_model <- function(mdl, mdl_name) {
    message(paste0("saving model to ", make_model_id_name(mdl_name), " ..."))
    readr::write_rds(mdl, make_model_id_name(mdl_name))
  }

  # parallel processing for grid search
  cl <- parallel::makePSOCKcluster(ncores)
  doParallel::registerDoParallel(cl)

  message("group 10-fold cv")
  control_object <- cross_fold_validation()

  # LINEAR MODELS
  mdl_name <- "linear"
  if (!has_model_been_trained(mdl_name)) {
    message("linear regression ...")
    set.seed(669)
    linear_reg <- caret::train(
      x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
      method = "lm",
      trControl = control_object,
      metric = "Rsquared"
    )
    save_model(linear_reg, mdl_name)
  }


  mdl_name <- "pls"
  if (!has_model_been_trained(mdl_name)) {
    message("pls regression ...")
    set.seed(669)
    pls_model <- caret::train(
      x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
      method = "pls",
      tuneLength = 15,
      trControl = control_object,
      metric = "Rsquared"
    )
    save_model(pls_model, mdl_name)
  }

  mdl_name <- "enet"
  if (!has_model_been_trained(mdl_name)) {
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
    save_model(enet_model, mdl_name)
  }


  mdl_name <- "earth"
  if (!has_model_been_trained(mdl_name)) {
    message("earth regression ...")
    set.seed(669)
    earth_model <- caret::train(
      x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
      method = "earth",
      tuneGrid = expand.grid(degree = 1, nprune = 2:25),
      trControl = control_object,
      metric = "Rsquared"
    )
    save_model(earth_model, mdl_name)
  }


  mdl_name <- "svm"
  if (!has_model_been_trained(mdl_name)) {
    message("svm regression ...")
    set.seed(669)
    svmr_model <- caret::train(
      x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
      method = "svmLinear",
      tuneLength = 10,
      trControl = control_object,
      metric = "Rsquared"
    )
    save_model(svmr_model, mdl_name)
  }


  mdl_name <- "nnet"
  if (!has_model_been_trained(mdl_name)) {
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
    save_model(nnet_model, mdl_name)
  }


  mdl_name <- "rpart"
  if (!has_model_been_trained(mdl_name)) {
    message("rpart  regression ...")
    set.seed(669)
    rpart_model <- caret::train(
      x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
      method = "rpart",
      tuneLength = 30,
      trControl = control_object,
      metric = "Rsquared"
    )
    save_model(rpart_model, mdl_name)
  }


  mdl_name <- "ctree"
  if (!has_model_been_trained(mdl_name)) {
    message("ctree  regression ...")
    set.seed(669)
    ctree_model <- caret::train(
      x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
      method = "ctree",
      tuneLength = 10,
      trControl = control_object,
      metric = "Rsquared"
    )
    save_model(ctree_model, mdl_name)
  }

  mdl_name <- "mtree"
  if (!has_model_been_trained(mdl_name)) {
    message("mtree  regression ...")
    set.seed(669)
    mt_model <- caret::train(
      x = datos_preprocessed$X_train, y = datos_preprocessed$y_train,
      method = "M5",
      trControl = control_object,
      metric = "Rsquared"
    )
    save_model(mt_model, mdl_name)
  }


  mdl_name <- "random_forest"
  if (!has_model_been_trained(mdl_name)) {
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
    save_model(rf_model, mdl_name)
  }

  mdl_name <- "gbm"
  if (!has_model_been_trained(mdl_name)) {
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
    save_model(gbm_model, mdl_name)
  }

  mdl_name <- "cubist"
  if (!has_model_been_trained(mdl_name)) {
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
    save_model(cubist_model)
  }

  # stop cluster
  parallel::stopCluster(cl)
}


#' Train Final Model
#'
#' This function trains the final model with
#' the final hyperparameters found after model
#' tunning
#'
#' @param X_train train predictors
#' @param y_train test response
#'
#' @return trained gradient boosted model with gaussian loss function.
#' @export
#' @note when using the output model for prediction, the
#' data should be a data.frame object and not a tibble, if not
#' the function will throw an error
train_final_model <- function(
                              X_train = datos_preprocessed$X_train,
                              y_train = datos_preprocessed$y_train) {

  # the X_train has be be a data.fram in order to preict
  # and train with gbm
  X_train <- base::as.data.frame(X_train)

  gbm::gbm.fit(
    # this model does no accept the tibble only the data frame
    x = X_train,
    y = y_train,
    distribution = "gaussian",
    n.trees = 750,
    interaction.depth = 5,
    shrinkage = 0.1,
    n.minobsinnode = 10,
    verbose = TRUE
  )
}

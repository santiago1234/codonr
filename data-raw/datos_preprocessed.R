# this code generates the preprocessing pipeline as well as the
# preprocessed train and test data

datos <- prepare_train_and_test_sets()
preprocessing_pipeline <- preprocessing(datos$X_train)

datos$X_train <- recipes::bake(preprocessing_pipeline, datos$X_train)
datos$X_test <- recipes::bake(preprocessing_pipeline, datos$X_test)
datos_preprocessed <- datos
usethis::use_data(datos_preprocessed)
usethis::use_data(preprocessing_pipeline)

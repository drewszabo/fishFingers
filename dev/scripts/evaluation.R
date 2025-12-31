# ---------------- Setup ------------------------------------------------------------------
library(tidyverse)
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(xgboost)

set.seed(42)

# ---------------- Set Model Control Preferences --------------------------------------

# caret training control
ctrl <- trainControl(method = "cv", number = 10)

# ---------------- Model 1: LASSO Regression (glmnet) --------------------------------------
lasso_grid <- expand.grid(alpha = 1, lambda = 10^seq(-3, 1, length = 10))
lasso_model <- train(
  x = train_x, y = train_y,
  method = "glmnet",
  trControl = ctrl,
  tuneGrid = lasso_grid
)
print(lasso_model)
plot(lasso_model)

# ---------------- Model 2: Random Forest --------------------------------------------------
rf_grid <- expand.grid(mtry = c(5, 10, 20, 50, 200))
rf_model <- train(
  x = train_x, y = train_y,
  method = "rf",
  trControl = ctrl,
  tuneGrid = rf_grid,
  ntree = 500
)
print(rf_model)
plot(rf_model)

# ---------------- Model 3: Support Vector Regression -------------------------------------
svm_grid <- expand.grid(C = 2^seq(-2, 5, 2), sigma = 0.01)
svm_model <- train(
  x = train_x, y = train_y,
  method = "svmRadial",
  trControl = ctrl,
  tuneGrid = svm_grid
)
print(svm_model)
plot(svm_model)

# ---------------- Model 4: k-NN -----------------------------------------------------------
knn_grid <- expand.grid(k = seq(3, 21, 2))
knn_model <- train(
  x = train_x, y = train_y,
  method = "knn",
  trControl = ctrl,
  tuneGrid = knn_grid
)
print(knn_model)
plot(knn_model)

# ---------------- Model 5: XGBoost --------------------------------------------------------

xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 6, 9),
  eta = c(0.01, 0.1, 0.3),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)
xgb_model <- train(
  x = train_x, y = train_y,
  method = "xgbTree",
  trControl = ctrl,
  tuneGrid = xgb_grid
)
print(xgb_model)
plot(xgb_model)

# ---------------- Model 6: Linear Regression ---------------------------------------------

lm_model <- train(
  x = train_x, y = train_y,
  method = "lm",
  trControl = ctrl
)

# ---------------- Model 7: Neural Networking -----------------------------------------


# ---------------- Evaluation -------------------------------------------------------------
predict_and_eval <- function(model, x, y, setname) {
  preds <- predict(model, newdata = x)
  rmse <- RMSE(preds, y)
  mae  <- MAE(preds, y)
  r2   <- R2(preds, y)
  tibble(Set = setname, RMSE = rmse, MAE = mae, R2 = r2)
}

results <- bind_rows(
  predict_and_eval(lm_model, train_x, train_y, "LM_train"),
  predict_and_eval(lm_model, test_x, test_y, "LM_test"),
  predict_and_eval(lasso_model, train_x, train_y, "LASSO_train"),
  predict_and_eval(lasso_model, test_x, test_y, "LASSO_test"),
  predict_and_eval(rf_model, train_x, train_y, "RF_train"),
  predict_and_eval(rf_model, test_x, test_y, "RF_test"),
  predict_and_eval(svm_model, train_x, train_y, "SVM_train"),
  predict_and_eval(svm_model, test_x, test_y, "SVM_test"),
  predict_and_eval(knn_model, train_x, train_y, "KNN_train"),
  predict_and_eval(knn_model, test_x, test_y, "KNN_test"),
  predict_and_eval(xgb_model, train_x, train_y, "XGB_train"),
  predict_and_eval(xgb_model, test_x, test_y, "XGB_test")
)

write_csv(results, "baseline_results.csv")
print(results)


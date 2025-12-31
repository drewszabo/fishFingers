# ---------------- Setup ------------------------------------------------------------------
library(tidyverse)
library(caret)
library(xgboost)
library(ModelMetrics)
library(ggplot2)
library(extrafont)
library(ggtext)

set.seed(42)

# ---------------- Set Model Control Preferences --------------------------------------

# caret training control
ctrl <- trainControl(method = "cv", number = 10)

# initial grid search parameters
xgb_grid1 <- expand.grid(
  nrounds = seq(from = 200, to = 1000, by = 50),
  eta = c(0.025, 0.05, 0.1, 0.3),
  max_depth = c(2, 3, 4, 5, 6),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

xgb_tune1 <- train(
  x = train_x, y = train_y,
  method = "xgbTree",
  trControl = ctrl,
  tuneGrid = xgb_grid1
)

print(xgb_tune1)
plot(xgb_tune1)
xgb_tune$bestTune

# tune grid parameters 2
xgb_grid2 <- expand.grid(
  nrounds = seq(from = 50, to = 1000, by = 50),
  eta = xgb_tune1$bestTune$eta,
  max_depth = ifelse(xgb_tune1$bestTune$max_depth == 2,
                     c(xgb_tune1$bestTune$max_depth:4),
                     xgb_tune1$bestTune$max_depth - 1:xgb_tune1$bestTune$max_depth + 1),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = c(1, 2, 3),
  subsample = 1
)

xgb_tune2 <- caret::train(
  x = train_x,
  y = train_y,
  trControl = ctrl,
  tuneGrid = xgb_grid2,
  method = "xgbTree"
)

plot(xgb_tune2)
xgb_tune2$bestTune

# tune grid parameters 3
xgb_grid3 <- expand.grid(
  nrounds = seq(from = 50, to = 1000, by = 50),
  eta = xgb_tune1$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = 0,
  colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = c(0.5, 0.75, 1.0)
)

xgb_tune3 <- caret::train(
  x = train_x,
  y = train_y,
  trControl = ctrl,
  tuneGrid = xgb_grid3,
  method = "xgbTree"
)

plot(xgb_tune3, probs = .95)
xgb_tune3$bestTune

save.image(".RData")

# tune grid parameters 4
tune_grid4 <- expand.grid(
  nrounds = seq(from = 50, to = 1000, by = 50),
  eta = xgb_tune1$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0),
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune4 <- caret::train(
  x = train_x,
  y = train_y,
  trControl = ctrl,
  tuneGrid = tune_grid4,
  method = "xgbTree"
)

plot(xgb_tune4)
xgb_tune4$bestTune

save.image(".RData")

# tune grid parameters 5

tune_grid5 <- expand.grid(
  nrounds = seq(from = 100, to = 5000, by = 100),
  eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = xgb_tune4$bestTune$gamma,
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune5 <- caret::train(
  x = train_x,
  y = train_y,
  trControl = ctrl,
  tuneGrid = tune_grid5,
  method = "xgbTree"
)

plot(xgb_tune5)
xgb_tune5$bestTune

evaluate_model(xgb_tune5)

save.image(".RData")

# Manual tuning to decrease complexity

manual_grid1 <- expand.grid(
  nrounds = xgb_tune5$bestTune$nrounds,
  eta = xgb_tune5$bestTune$eta,
  max_depth = xgb_tune5$bestTune$max_depth,
  gamma = c(0.5,1,2,3,4,5),
  colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
  min_child_weight = xgb_tune5$bestTune$min_child_weight,
  subsample = xgb_tune5$bestTune$subsample
)

xgb_tune6 <- caret::train(
  x = train_x,
  y = train_y,
  trControl = ctrl,
  tuneGrid = manual_grid1,
  method = "xgbTree"
)

plot(xgb_tune6)


# Force max depth lower and min child higher

manual_grid2 <- expand.grid(
  nrounds = xgb_tune5$bestTune$nrounds,
  eta = xgb_tune5$bestTune$eta,
  max_depth = c(3, 4, 5),
  gamma = xgb_tune6$bestTune$gamma,
  colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
  min_child_weight = c(4, 5, 6, 7, 8),
  subsample = xgb_tune5$bestTune$subsample
)

xgb_tune7 <- caret::train(
  x = train_x,
  y = train_y,
  trControl = ctrl,
  tuneGrid = manual_grid2,
  method = "xgbTree"
)

plot(xgb_tune7)

evaluate_model(xgb_tune7)

# final model training

final_grid <- expand.grid(
  nrounds = xgb_tune5$bestTune$nrounds,
  eta = xgb_tune5$bestTune$eta,
  max_depth = xgb_tune7$bestTune$max_depth,
  gamma = xgb_tune6$bestTune$gamma,
  colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
  min_child_weight = xgb_tune7$bestTune$min_child_weight,
  subsample = xgb_tune5$bestTune$subsample
)

xgb_model <- caret::train(
  x = train_x,
  y = train_y,
  trControl = ctrl,
  tuneGrid = final_grid,
  method = "xgbTree"
)

save.image(".RData")
saveRDS(xgb_model, "final_model.rds")

train_cols <- colnames(train_df)  # after dummyVars or preprocessing
xgb_model$finalModel$feature_names <- train_cols
saveRDS(xgb_model, "final_model.rds")

evaluate_model(xgb_model, cv = TRUE)

evaluate_model <- function(
    model,
    train_x = NULL,
    train_y = NULL,
    cv = FALSE
) {
  
  # Use global objects if not provided
  if (is.null(train_x)) train_x <- get("train_x", envir = parent.frame())
  if (is.null(train_y)) train_y <- get("train_y", envir = parent.frame())
  
  # Predictions
  pred <- predict(model, train_x)

  # Metrics
  pred_rmse <- rmse(train_y, pred)

  pred_mae  <- mae(train_y, pred)

  pred_r2   <- cor(train_y, pred)^2

  if(cv == TRUE) {
    
    # Extract best CV RMSE from caret model
    best_params <- model$bestTune
    
    cv_rmse <- model$results %>%
      dplyr::filter(
        nrounds == best_params$nrounds,
        max_depth == best_params$max_depth,
        eta == best_params$eta,
        gamma == best_params$gamma,
        colsample_bytree == best_params$colsample_bytree,
        min_child_weight == best_params$min_child_weight,
        subsample == best_params$subsample
      ) %>%
      dplyr::pull(RMSE)
    
  }
  
  # Return everything as a list
  return(list(
    rmse = pred_rmse,
    rmse_cv = if(cv == TRUE){cv_rmse}else{NULL},
    mae  = pred_mae,
    r2   = pred_r2
  ))
}

# Predictions
train_pred <- predict(xgb_model, train_x)
test_pred  <- predict(xgb_model, test_x)

plot_train <- data.frame(measured = train_y, predicted = train_pred) %>%
  bind_cols(select(train_df, species))

# annotate measured and predicted values
plot_train <- plot_train %>%
  mutate(measured_ann = case_when(measured < log10(2000) ~ "nB",
                                  measured >= log10(5000) ~ "vB",
                                  TRUE ~ "B"),
         predicted_ann = case_when(predicted < log10(2000) ~ "nB",
                                   predicted >= log10(5000) ~ "vB",
                                   TRUE ~ "B"))

# Convert to ordinal factors (numerical for ranking)
plot_train$measured_ann <- factor(plot_train$measured_ann,
                                  levels = c("nB", "B", "vB"),
                                  ordered = TRUE)

plot_train$predicted_ann <- factor(plot_train$predicted_ann,
                                   levels = c("nB", "B", "vB"),
                                   ordered = TRUE)

# Convert to numeric ranks
true  <- as.numeric(plot_train$measured_ann)
pred  <- as.numeric(plot_train$predicted_ann)

# Calculate confustion matrix
cm <- caret::confusionMatrix(plot_train$predicted_ann,
                      plot_train$measured_ann)

cm

# Calculate the quadratic weighted kappa (QWK)
qwk <- irr::kappa2(
  data.frame(true, pred),
  weight = "squared"
)

qwk

# Find the top 3 most common species
top3 <- plot_train %>%
  count(species, sort = TRUE) %>%
  slice_head(n = 3) %>%
  pull(species)

# Create the new column
plot_train <- plot_train %>%
  mutate(top_species = if_else(species %in% top3, species, "Other"))

label_text <- paste0(
  "RMSE = ", round(model_metrics$train_rmse, 2), "<br>",
  "RMSE (CV) = ", round(model_metrics$cv_rmse, 2), "<br>",
  "MAE = ", round(model_metrics$train_mae, 2), "<br>",
  "<i>R</i><sup>2</sup> = ", round(model_metrics$train_r2, 2)
)

p1 <- ggplot(data = plot_train, aes(x = measured, y = predicted, colour = factor(top_species, levels = c("Cyprinus_carpio", "Oncorhynchus_mykiss", "Poecilia_reticulata", "Other")))) +
  geom_point(size = 1, shape = 1, stroke = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_abline(slope = 1, intercept = 1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = 2, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -2, linetype = "dotted", color = "grey40") +
  scale_x_continuous(name = "Measured logBCF",
                     limits = c(-2, 6)) +
  scale_y_continuous(name = "Predicted logBCF",
                     limits = c(-1, 6)) +
  scale_color_discrete(
    name = "",
    labels = list(
      expression(italic("C. carpio")),
      expression(italic("O. mykiss")),
      expression(italic("P. reticulata")),
      "Other"
    )
  ) +
  labs(title = "Train") +
  theme_bw() +
  theme(legend.position = "left",
        #axis.title.x = element_blank(),
        #axis.text.x = element_blank()
        ) +
  guides(colour = guide_legend(nrow = 2)) +
  coord_equal(ratio = 1)

p1

ggsave("train_plot.svg",
       plot = p1,
       dpi = 1200,
       height = 4,
       width = 6,
       unit = "in")


plot_test <- data.frame(measured = test_y, predicted = test_pred) %>%
  bind_cols(select(test_df, species))

# Create the new column
plot_test <- plot_test %>%
  mutate(top_species = if_else(species %in% top3, species, "Other"))

# annotate measured and predicted values
plot_test <- plot_test %>%
  mutate(measured_ann = case_when(measured < log10(2000) ~ "nB",
                                  measured >= log10(5000) ~ "vB",
                                  TRUE ~ "B"),
         predicted_ann = case_when(predicted < log10(2000) ~ "nB",
                                   predicted >= log10(5000) ~ "vB",
                                   TRUE ~ "B"))

# Convert to ordinal factors (numerical for ranking)
plot_test$measured_ann <- factor(plot_test$measured_ann,
                                  levels = c("nB", "B", "vB"),
                                  ordered = TRUE)

plot_test$predicted_ann <- factor(plot_test$predicted_ann,
                                   levels = c("nB", "B", "vB"),
                                   ordered = TRUE)

# Convert to numeric ranks
true  <- as.numeric(plot_test$measured_ann)
pred  <- as.numeric(plot_test$predicted_ann)

# Calculate confustion matrix
cm <- caret::confusionMatrix(plot_test$predicted_ann,
                             plot_test$measured_ann)

cm

# Calculate the quadratic weighted kappa (QWK)
qwk <- irr::kappa2(
  data.frame(true, pred),
  weight = "squared"
)

qwk

evaluate_model(xgb_model, test_x, test_y)

p2 <- ggplot(data = plot_test, aes(x = measured, y = predicted, colour = factor(top_species, levels = c("Cyprinus_carpio", "Oncorhynchus_mykiss", "Poecilia_reticulata", "Other")))) +
  geom_point(size = 1, shape = 1, stroke = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_abline(slope = 1, intercept = 1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -1, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = 2, linetype = "dotted", color = "grey40") +
  geom_abline(slope = 1, intercept = -2, linetype = "dotted", color = "grey40") +
  scale_x_continuous(name = "Measured logBCF",
                     limits = c(-2, 6)) +
  scale_y_continuous(name = "Predicted logBCF",
                     limits = c(-1, 6)) +
  scale_color_discrete(name = "",
                       labels = c(expression(italic("C. carpio")), expression(italic("O. mykiss")), expression(italic("P. reticulata")), "Other")) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Test") +
  coord_equal(ratio = 1)

p2


ggsave("test_plot.svg",
       plot = p2,
       dpi = 1200,
       height = 4,
       width = 4,
       unit = "in")

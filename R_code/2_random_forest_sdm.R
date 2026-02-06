# ============================================================
# Random Forest SDM w/ spatial CV + tuning (ranger)
# + variable importance, PDPs, and prediction raster
# ============================================================

# -------------------------------
# Packages
# -------------------------------
library(dplyr)
library(tidyr)
library(ranger)
library(ggplot2)
library(raster)
library(pdp)
library(pROC)

# -------------------------------
# Load data
# -------------------------------
analysis_data <- read.csv("data/a_chamek_ter_mammals_finalData_Oct22.csv")

# Keep only columns needed (including fold for CV)
analysis_data_v2 <- analysis_data %>%
  dplyr::select(
    presence,
    fold,
    bio13_precip_wettest_month,
    cmi_min,
    mean_forest,
    mean_farming,
    diff_forest_formation
  ) %>%
  # Shuffle rows (helps avoid ordering artifacts)
  dplyr::slice_sample(prop = 1)

# -------------------------------
# Cross-validation + tuning setup
# -------------------------------
fold_ids <- sort(unique(analysis_data_v2$fold))  # works even if folds aren’t 1:3
n_folds  <- length(fold_ids)

rf_performance <- data.frame(
  model      = rep("RF", n_folds),
  fold_id    = fold_ids,
  auc        = rep(NA_real_, n_folds),
  presence   = rep(NA_integer_, n_folds),
  background = rep(NA_integer_, n_folds)
)

hypergrid_final <- data.frame(
  mtry          = rep(NA_integer_, n_folds),
  node_size     = rep(NA_integer_, n_folds),
  sample_frac   = rep(NA_real_,    n_folds)
)

# Hyperparameter grid (typos fixed: sample_size -> sample_frac)
hyper_grid <- expand.grid(
  mtry        = seq(1, 3, by = 1),
  node_size   = seq(1, 4, by = 1),
  sample_frac = c(0.60, 0.70, 0.80),
  OOB_RMSE    = 0
)

# -------------------------------
# CV loop (robust to masking)
# -------------------------------
fold_ids <- sort(unique(analysis_data_v2$fold))

for (k in seq_along(fold_ids)) {
  
  i <- fold_ids[k]
  
  train <- analysis_data_v2 %>%
    dplyr::filter(fold != i) %>%
    dplyr::select(-fold)
  
  test <- analysis_data_v2 %>%
    dplyr::filter(fold == i) %>%
    dplyr::select(-fold)
  
  train_complete <- train[complete.cases(train), ]
  test_complete  <- test[complete.cases(test), ]
  
  hyper_grid$OOB_RMSE <- 0
  
  for (j in seq_len(nrow(hyper_grid))) {
    
    model_j <- ranger(
      presence ~ .,
      data            = train_complete,
      num.trees       = 2000,
      mtry            = hyper_grid$mtry[j],
      min.node.size   = hyper_grid$node_size[j],
      sample.fraction = hyper_grid$sample_frac[j],
      probability     = TRUE,
      replace         = TRUE,
      splitrule       = "hellinger",
      seed            = 123
    )
    
    hyper_grid$OOB_RMSE[j] <- sqrt(model_j$prediction.error)
  }
  
  hyper_grid2 <- hyper_grid %>% dplyr::arrange(OOB_RMSE)
  
  train_model <- ranger(
    presence ~ .,
    data            = train_complete,
    num.trees       = 2000,
    mtry            = hyper_grid2$mtry[1],
    min.node.size   = hyper_grid2$node_size[1],
    sample.fraction = hyper_grid2$sample_frac[1],
    probability     = TRUE,
    replace         = TRUE,
    splitrule       = "hellinger",
    seed            = 123
  )
  
  pred_obj <- predict(train_model, data = test_complete)
  pred_mat <- pred_obj$predictions
  
  # choose probability for class "1" if present; else fall back to column 2
  if (is.matrix(pred_mat) && ("1" %in% colnames(pred_mat))) {
    pred <- pred_mat[, "1"]
  } else if (is.matrix(pred_mat)) {
    pred <- pred_mat[, 2]
  } else {
    pred <- pred_mat
  }
  
  roc_obj <- pROC::roc(
    response  = test_complete$presence,
    predictor = pred,
    levels    = c(0, 1),
    quiet     = TRUE
  )
  
  rf_performance$auc[k]        <- as.numeric(roc_obj$auc)
  rf_performance$presence[k]   <- sum(test$presence == 1, na.rm = TRUE)
  rf_performance$background[k] <- sum(test$presence == 0, na.rm = TRUE)
  
  hypergrid_final$mtry[k]        <- hyper_grid2$mtry[1]
  hypergrid_final$node_size[k]   <- hyper_grid2$node_size[1]
  hypergrid_final$sample_frac[k] <- hyper_grid2$sample_frac[1]
}

# Mean CV AUC
mean_auc <- mean(rf_performance$auc, na.rm = TRUE)
mean_auc

# -------------------------------
# Train final model (full data)
# -------------------------------
final_df <- analysis_data_v2[complete.cases(analysis_data_v2), ] %>%
  dplyr::select(-fold)

# Use average best hyperparams from CV
final_mtry        <- round(mean(hypergrid_final$mtry, na.rm = TRUE))
final_node_size   <- round(mean(hypergrid_final$node_size, na.rm = TRUE))
final_sample_frac <- mean(hypergrid_final$sample_frac, na.rm = TRUE)

final_model <- ranger(
  presence ~ .,
  data            = final_df,
  num.trees       = 2000,
  mtry            = final_mtry,
  min.node.size   = final_node_size,
  sample.fraction = final_sample_frac,
  probability     = TRUE,
  replace         = TRUE,
  splitrule       = "hellinger",
  importance      = "permutation",
  seed            = 123
)

# In-sample AUC (sanity check for overfitting)
pred_obj <- predict(final_model, data = final_df)
pred_mat <- pred_obj$predictions
if (is.matrix(pred_mat) && ("1" %in% colnames(pred_mat))) {
  pred <- pred_mat[, "1"]
} else if (is.matrix(pred_mat)) {
  pred <- pred_mat[, 2]
} else {
  pred <- pred_obj$predictions
}

in_sample_auc <- pROC::roc(
  response  = final_df$presence,
  predictor = pred,
  levels    = c(0, 1),
  quiet     = TRUE
)$auc

in_sample_auc

# -------------------------------
# Variable importance plot
# -------------------------------
permutation_importance <- data.frame(
  variable   = names(final_model$variable.importance),
  importance = as.vector(final_model$variable.importance)
) %>%
  arrange(desc(importance))

ggplot(permutation_importance, aes(x = reorder(variable, importance), y = importance)) +
  geom_col() +
  coord_flip() +
  labs(title = "Permutation importance", x = NULL, y = "Importance") +
  theme_classic(base_size = 14)

# -------------------------------
# Partial Dependence Plots (PDPs)
# -------------------------------
var_names <- setdiff(names(final_df), "presence")

pd_df <- data.frame(variable = character(), value = numeric(), yhat = numeric())

for (v in var_names) {
  out <- as.data.frame(
    pdp::partial(
      final_model,
      pred.var = v,
      prob     = TRUE,
      train    = final_df
    )
  )
  
  # pdp output: first col = grid for v; second col = yhat
  loop_df <- data.frame(
    variable = v,
    value    = out[[1]],
    yhat     = out[[2]]
  )
  
  pd_df <- dplyr::bind_rows(pd_df, loop_df)
}

ggplot(pd_df, aes(x = value, y = yhat)) +
  geom_smooth() +
  facet_wrap(~variable, scales = "free") +
  labs(y = "Predicted probability", x = NULL) +
  theme_bw(base_size = 14)

# -------------------------------
# Prediction map (raster stack -> GeoTIFF)
# -------------------------------
env_data <- list.files(
  path = "env_data",
  pattern = "tif$",
  full.names = TRUE,
  recursive = TRUE
)

e <- raster::stack(env_data)

prediction_df <- as.data.frame(raster::rasterToPoints(e)) %>%
  # rename these if the stack layer names don’t match your covariate names
  rename(
    bio13_precip_wettest_month = mean.1,
    cmi_min                    = mean.2
  )

prediction_df_complete <- prediction_df[complete.cases(prediction_df), ]

predictions <- predict(final_model, data = prediction_df_complete)

prediction_df_full <- cbind(prediction_df_complete, probability = predictions)

rf_tiff_df <- prediction_df_full[, c("x", "y", "probability.X1")]

rf_sdm_raster <- raster::rasterFromXYZ(rf_tiff_df)

# quick visualization
plot(rf_sdm_raster)

# if that doesn't work, tell R you want to use the raster library for plotting
# raster::plot(rf_sdm_raster)

writeRaster(
  rf_sdm_raster,
  filename  = "final_figures/rf_sdm_example_predictions.tif",
  format    = "GTiff",
  options   = c("INTERLEAVE=BAND", "COMPRESS=LZW"),
  overwrite = TRUE
)

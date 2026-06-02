run_lasso_nested_cv <- function(data, outcome_bin_col, covariate_cols, metabolite_cols,
                                alpha = 1, n_outer = 10, n_inner = 10, seed = 123) {
  
  library(glmnet)
  library(pROC)
  library(caret)
  
  set.seed(seed)
  
  n_samples <- nrow(data)
  if (length(unique(data[[outcome_bin_col]])) != 2) {
    stop("Outcome must be binary with exactly 2 levels.")
  }
  
  # ---- Covariates and metabolites model matrices (no folds yet) ----
  x_cov_raw <- model.matrix(~ . - 1, data[, covariate_cols, drop = FALSE])
  x_met_raw <- model.matrix(~ . - 1, data[, metabolite_cols, drop = FALSE])
  y_var <- as.numeric(as.factor(data[[outcome_bin_col]])) - 1
  
  # Identify numeric covariates (for scaling)
  cov_is_numeric <- sapply(as.data.frame(x_cov_raw), is.numeric)
  
  # ---- Global scaling (once, outside folds) ----
  x_cov <- x_cov_raw
  if (any(cov_is_numeric)) {
    sc_cov <- scale(x_cov_raw[, cov_is_numeric, drop = FALSE])
    x_cov[, cov_is_numeric] <- sc_cov
  }
  
  x_met <- scale(x_met_raw)
  
  # Outer folds
  outer_folds <- createFolds(y_var, k = n_outer, list = TRUE, returnTrain = FALSE)
  
  outer_pred <- rep(NA_real_, n_samples)
  selected_features_all <- vector("list", n_outer)
  
  # Penalty factors: 0 for covariates, 1 for metabolites
  penalty_factors <- c(rep(0, ncol(x_cov)), rep(1, ncol(x_met)))
  
  for (i in seq_along(outer_folds)) {
    test_idx <- outer_folds[[i]]
    train_idx <- setdiff(seq_len(n_samples), test_idx)
    
    # Split (already globally scaled; no refit of scaling inside folds)
    x_train <- cbind(x_cov[train_idx, ], x_met[train_idx, ])
    x_test  <- cbind(x_cov[test_idx,  ], x_met[test_idx,  ])
    
    # Inner CV
    inner_cv <- cv.glmnet(
      x           = x_train,
      y           = y_var[train_idx],
      family      = "binomial",
      alpha       = alpha,
      type.measure = "deviance",
      nfolds      = n_inner,
      penalty.factor = penalty_factors
    )
    
    # Predict using lambda.1se (consistent with final model)
    pred_test <- predict(inner_cv, newx = x_test, s = "lambda.1se", type = "response")
    outer_pred[test_idx] <- as.vector(pred_test)
    
    # Collect selected metabolite features at lambda.1se
    coef_1se <- coef(inner_cv, s = "lambda.1se")
    sel_all  <- rownames(coef_1se)[which(coef_1se != 0)]
    sel_met  <- setdiff(sel_all, c("(Intercept)", colnames(x_cov)))
    selected_features_all[[i]] <- sel_met
  }
  
  # ROC and metrics
  roc_obj <- pROC::roc(
    response   = y_var,
    predictor  = outer_pred,
    levels     = c(0, 1),
    direction  = "<"
  )
  
  opt_thresh <- as.numeric(
    pROC::coords(roc_obj, "best", input = "threshold", best.method = "youden")["threshold"]
  )
  
  pred_class <- ifelse(outer_pred > opt_thresh, 1, 0)
  
  valid_idx <- complete.cases(data.frame(pred_class, y_var))
  conf_matrix <- caret::confusionMatrix(
    factor(pred_class[valid_idx], levels = c(0, 1)),
    factor(y_var[valid_idx], levels = c(0, 1)),
    positive = "1"
  )
  
  metrics <- data.frame(
    Metric = c("AUC", "Sensitivity", "Specificity", "Accuracy", "Threshold"),
    Value = c(
      as.numeric(pROC::auc(roc_obj)),
      as.numeric(conf_matrix$byClass["Sensitivity"]),
      as.numeric(conf_matrix$byClass["Specificity"]),
      as.numeric(conf_matrix$overall["Accuracy"]),
      opt_thresh
    )
  )
  
  # ---- Final model on all data (same global scaling) ----
  final_x <- cbind(x_cov, x_met)
  
  final_cv <- cv.glmnet(
    x           = final_x,
    y           = y_var,
    family      = "binomial",
    alpha       = alpha,
    type.measure = "deviance",
    nfolds      = n_inner,
    penalty.factor = penalty_factors
  )
  
  # Use lambda.1se in final model (not lambda.min)
  final_coef <- coef(final_cv, s = "lambda.1se")
  
  final_selected_met <- setdiff(
    setdiff(rownames(final_coef)[which(final_coef != 0)], "(Intercept)"),
    colnames(x_cov)
  )
  
  list(
    metrics = metrics,
    final_selected_metabolites = final_selected_met,
    n_final_metabolites = length(final_selected_met),
    roc_obj = roc_obj,
    final_cv = final_cv,
    outer_predictions = outer_pred,
    selected_features_all = selected_features_all
  )
}

run_lasso_nested_cv_modelA <- function(data, outcome_bin_col, metabolite_cols,
                                       alpha = 1, n_outer = 10, n_inner = 10, seed = 123) {
  
  library(glmnet)
  library(pROC)
  library(caret)
  
  set.seed(seed)
  
  n_samples <- nrow(data)
  if (length(unique(data[[outcome_bin_col]])) != 2) {
    stop("Outcome must be binary with exactly 2 levels.")
  }
  
  # Metabolites only
  x_met_raw <- model.matrix(~ . - 1, data[, metabolite_cols, drop = FALSE])
  y_var <- as.numeric(as.factor(data[[outcome_bin_col]])) - 1
  
  # ---- Global scaling (once, outside folds) ----
  x_met <- scale(x_met_raw)
  
  # Outer folds
  outer_folds <- createFolds(y_var, k = n_outer, list = TRUE, returnTrain = FALSE)
  
  outer_pred <- rep(NA_real_, n_samples)
  selected_features_all <- vector("list", n_outer)
  
  for (i in seq_along(outer_folds)) {
    test_idx <- outer_folds[[i]]
    train_idx <- setdiff(seq_len(n_samples), test_idx)
    
    # Split (already scaled globally)
    x_train <- x_met[train_idx, ]
    x_test  <- x_met[test_idx,  ]
    
    # Inner CV
    inner_cv <- cv.glmnet(
      x = x_train,
      y = y_var[train_idx],
      family = "binomial",
      alpha = alpha,
      type.measure = "deviance",
      nfolds = n_inner
    )
    
    # Predict using lambda.1se
    pred_test <- predict(inner_cv, newx = x_test, s = "lambda.1se", type = "response")
    outer_pred[test_idx] <- as.vector(pred_test)
    
    # Collect selected metabolites at lambda.1se
    coef_1se <- coef(inner_cv, s = "lambda.1se")
    sel_all  <- rownames(coef_1se)[which(coef_1se != 0)]
    sel_met  <- setdiff(sel_all, "(Intercept)")
    selected_features_all[[i]] <- sel_met
  }
  
  # ROC and metrics
  roc_obj <- pROC::roc(
    response   = y_var,
    predictor  = outer_pred,
    levels     = c(0, 1),
    direction  = "<"
  )
  
  opt_thresh <- as.numeric(
    pROC::coords(roc_obj, "best", input = "threshold", best.method = "youden")["threshold"]
  )
  
  pred_class <- ifelse(outer_pred > opt_thresh, 1, 0)
  
  valid_idx <- complete.cases(data.frame(pred_class, y_var))
  conf_matrix <- caret::confusionMatrix(
    factor(pred_class[valid_idx], levels = c(0, 1)),
    factor(y_var[valid_idx], levels = c(0, 1)),
    positive = "1"
  )
  
  metrics <- data.frame(
    Metric = c("AUC", "Sensitivity", "Specificity", "Accuracy", "Threshold"),
    Value = c(
      as.numeric(pROC::auc(roc_obj)),
      as.numeric(conf_matrix$byClass["Sensitivity"]),
      as.numeric(conf_matrix$byClass["Specificity"]),
      as.numeric(conf_matrix$overall["Accuracy"]),
      opt_thresh
    )
  )
  
  # Final model on all data (same global scaling)
  final_cv <- cv.glmnet(
    x = x_met,
    y = y_var,
    family = "binomial",
    alpha = alpha,
    type.measure = "deviance",
    nfolds = n_inner
  )
  
  # Use lambda.1se everywhere (including final selection)
  final_coef <- coef(final_cv, s = "lambda.1se")
  
  final_selected_met <- setdiff(
    rownames(final_coef)[which(final_coef != 0)],
    "(Intercept)"
  )
  
  list(
    metrics = metrics,
    final_selected_metabolites = final_selected_met,
    n_final_metabolites = length(final_selected_met),
    roc_obj = roc_obj,
    final_cv = final_cv,
    outer_predictions = outer_pred,
    selected_features_all = selected_features_all
  )
}
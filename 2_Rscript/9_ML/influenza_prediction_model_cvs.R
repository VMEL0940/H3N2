library(reticulate)
use_condaenv("r-tf210", required = TRUE)

library(tensorflow)
tf$`__version__`         # 2.10.0
tf$config$list_physical_devices("GPU")

library(keras)
library(dplyr)
library(ggplot2)
library(caret)
library(e1071)
library(glmnet)
library(randomForest)
library(pROC) 
library(lme4)  
library(tidyr)
library(ggrepel)

k_clear_session()

path <- "~/H3N2/1_Data/9_ML"
df <- read.csv(paste0(path, "H3N2_preprocessed.csv")) %>%
  filter(vaccine_code != "Dar21") %>%
  mutate(NA_NS_Nonsyn = NA_Nonsyn * NS_Nonsyn)

get_train_test_split <- function(df, train_start, train_end, test_value) {
  vaccine_list <- unique(df$vaccine_code)
  train_codes  <- vaccine_list[which(vaccine_list == train_start):which(vaccine_list == train_end)]
  list(
    train       = df %>% filter(vaccine_code %in% train_codes),
    test        = df %>% filter(vaccine_code == test_value),
    train_codes = train_codes,
    test_code   = test_value
  )
}

get_cv_folds <- function(codes, window = 4) {
  folds <- list()
  n <- length(codes)
  for (start in seq_len(n - window)) {
    train_fold <- codes[start:(start + window - 1)]
    val_code   <- codes[start + window]
    folds[[length(folds) + 1]] <- list(train = train_fold, val = val_code)
  }
  folds
}

fit_model <- function(model, X, y, params) {
  if (model == "SVM") {
    svm(x = X, y = y, kernel = "radial",
        cost = params$C, gamma = params$gamma,
        probability = TRUE)
  } else if (model == "RF") {
    randomForest(x = X, y = y,
                 mtry = params$mtry,
                 ntree = params$ntree)
  } else {
    stop("Unknown model: ", model)
  }
}

predict_prob <- function(fit, X, model) {
  if (model == "SVM") {
    attr(predict(fit, X, probability = TRUE), "probabilities")[,2]
  } else if (model == "RF") {
    predict(fit, X, type = "prob")[,2]
  }
}

auc_metric <- tf$keras$metrics$AUC(name = "auc")
es_cb <- callback_early_stopping(
  monitor           = "val_auc", 
  mode              = "max",
  patience          = 5,
  restore_best_weights = TRUE
)

splits <- list(
  get_train_test_split(df, "Mos99", "Swtz13", "HK15"),
  get_train_test_split(df, "Fuj02", "HK15", "Kan17"),
  get_train_test_split(df, "Cal04", "Kan17", "HK19"),
  get_train_test_split(df, "Wis05", "HK19", "Cam20")
)

Syn_vars <- names(df)[grep("Syn", names(df))]
Nonsyn_vars <- setdiff(names(df)[grep("Nonsyn", names(df))], "NA_NS_Nonsyn")

analysis_schemes <- list(
  Original_var = c("HA_Nonsyn", "NA_Nonsyn", "NS_Nonsyn", "NA_NS_Nonsyn"),
  HA_only = c("HA_Nonsyn"),
  Nonsyn_var = Nonsyn_vars,
  Full_var = c(Syn_vars, Nonsyn_vars)
)

analysis_name <- list(
  Original_var = "Original",
  HA_only      = "HA Non-synonymous Distances",
  Nonsyn_var     = "All Non-synonymous Distances",
  Full_var     = "All Distances"
)

selected_params <- data.frame(
  Scheme      = character(), Split = integer(), Train_Range = character(),
  Test_Set    = character(), Model = character(), Param = character(), Value = numeric(),
  stringsAsFactors = FALSE
)

default_model_order <- c("GLMM", "Logistic", "Ridge", "Lasso", "SVM", "RF", "DNN")

col_models <- c("snow4", "tomato", "cornflowerblue", "orange", "mediumvioletred", "darkgreen", "midnightblue")
names(col_models) <- default_model_order

results_all <- list()
plots_all <- list()

results_long_merge <- NULL


library(reshape2)


df_long <- df %>% melt(., id.vars = setdiff(names(df), c(Syn_vars, Nonsyn_vars))) %>% 
  mutate(vaccine_code = factor(vaccine_code, levels = unique(as.character(vaccine_code))),
         Dom = factor(Dom, levels = c(0, 1))) %>% 
  mutate(variable = factor(variable, levels = c(Syn_vars, Nonsyn_vars)))

p_0_1_1 <- ggplot(df_long %>% filter(variable %in% Syn_vars), aes(x = vaccine_code, y = value, fill = Dom)) +
  geom_boxplot(position = position_dodge2()) +
  facet_wrap(facets = vars(variable), ncol = 1, scales = "free_y") +
  theme_bw()

ggsave(plot = p_0_1_1, filename = paste0("EDA_0-1_Syn_box_tot.png"),
       path = paste0(path, "Result/"), dpi = 300, width = 12, height = 8)

p_0_1_2 <- ggplot(df_long %>% filter(variable %in% Nonsyn_vars), aes(x = vaccine_code, y = value, fill = Dom)) +
  geom_boxplot(position = position_dodge2()) +
  facet_wrap(facets = vars(variable), ncol = 1, scales = "free_y") +
  theme_bw()

ggsave(plot = p_0_1_2, filename = paste0("EDA_0-2_Nonsyn_box_tot.png"),
       path = paste0(path, "Result/"), dpi = 300, width = 12, height = 8)


df_long_ch <- df_long %>% mutate(var_type = ifelse(variable %in% Syn_vars, "Synonymous", "Non-Synonymous")) %>% 
  mutate(var_type = factor(var_type, levels = c("Synonymous", "Non-Synonymous"))) %>% 
  mutate(variable = gsub("_Nonsyn", "", gsub("_Syn", "", as.character(variable)))) %>% 
  mutate(variable = factor(variable, levels = unique(as.character(variable))))

p_0_1_3 <- ggplot(df_long_ch,
                  aes(x = variable, y = value, fill = Dom)) +
  geom_boxplot(position = position_dodge2()) +
  facet_wrap(facets = vars(var_type), ncol = 1, scales = "free_y") +
  theme_bw()
p_0_1_3
ggsave(plot = p_0_1_3, filename = paste0("EDA_0-3_marginal_box_tot.png"),
       path = paste0(path, "Result/"), dpi = 300, width = 12, height = 8)

for (scheme in names(analysis_schemes)) {
  predictors <- analysis_schemes[[scheme]]
  
  res <- data.frame(
    Train_Range = c("Mos99-Swtz13", "Fuj02-HK15", "Cal04-Kan17", "Wis05-HK19"),
    Test_Set = c("HK15", "Kan17", "HK19", "Cam20"),
    Logistic_Train_AUROC = NA, Logistic_Test_AUROC = NA,
    Ridge_Train_AUROC    = NA, Ridge_Test_AUROC    = NA,
    Lasso_Train_AUROC    = NA, Lasso_Test_AUROC    = NA,
    SVM_Train_AUROC      = NA, SVM_Test_AUROC      = NA,
    RF_Train_AUROC       = NA, RF_Test_AUROC       = NA,
    DNN_Train_AUROC      = NA, DNN_Test_AUROC      = NA,
    GLMM_Train_AUROC     = NA, GLMM_Test_AUROC     = NA
  )
  for (i in seq_along(splits)) {
    split      <- splits[[i]]
    train_df   <- split$train %>% mutate(Dom = factor(Dom, levels = c("0","1")))
    test_df    <- split$test  %>% mutate(Dom = factor(Dom, levels = c("0","1")))
    train_codes <- split$train_codes
    test_code   <- split$test_code
    # train/test X, y
    X_full <- as.matrix(train_df %>% select(all_of(predictors)))
    y_full <- train_df$Dom
    X_test <- as.matrix(test_df %>% select(all_of(predictors)))
    y_test <- test_df$Dom
    
    # Nested CV for folds
    folds <- get_cv_folds(train_codes, window = 5)
    
    #### 1. Logistic Regression
    set.seed(777)
    logit_model <- glm(Dom ~ ., data = train_df %>% select(all_of(predictors), Dom), family = "binomial")
    logit_train_pred <- predict(logit_model, newdata = train_df, type = "response")
    logit_test_pred  <- predict(logit_model, newdata = test_df,  type = "response")
    logit_train_auc  <- auc(roc(y_full, logit_train_pred))
    logit_test_auc   <- auc(roc(y_test,  logit_test_pred))
    

    #### 2. Ridge & Lasso
    if (length(predictors) > 1) {
      # Ridge: alpha = 0
      ridge_full <- glmnet(X_full, y_full, alpha = 0, family = "binomial")
      lam_ridge <- ridge_full$lambda
      
      agg_auc_ridge <- numeric(length(lam_ridge))
      for (f in seq_along(folds)) {
        df_tr <- df %>% filter(vaccine_code %in% folds[[f]]$train) %>% mutate(Dom = factor(Dom, levels = c("0","1")))
        df_va <- df %>% filter(vaccine_code == folds[[f]]$val)   %>% mutate(Dom = factor(Dom, levels = c("0","1")))
        X_tr  <- as.matrix(df_tr %>% select(all_of(predictors)))
        y_tr  <- df_tr$Dom
        X_va  <- as.matrix(df_va %>% select(all_of(predictors)))
        y_va  <- df_va$Dom
        # glmnet fits all lambdas at once
        mod_tr <- glmnet(X_tr, y_tr, alpha = 0, family = "binomial", lambda = lam_ridge)
        preds_mat <- predict(mod_tr, X_va, type = "response")
        for (j in seq_along(lam_ridge)) {
          agg_auc_ridge[j] <- agg_auc_ridge[j] + auc(roc(y_va, preds_mat[, j]))
        }
      }
      mean_auc_ridge <- agg_auc_ridge / length(folds)
      best_idx_ridge <- which.max(mean_auc_ridge)
      best_lambda_ridge <- lam_ridge[best_idx_ridge]
      final_ridge <- glmnet(X_full, y_full, alpha = 0, family = "binomial", lambda = best_lambda_ridge)
      ridge_train_pred <- as.vector(predict(final_ridge, X_full, type = "response"))
      ridge_train_auc  <- auc(roc(y_full, ridge_train_pred))
      ridge_test_pred <- as.vector(predict(final_ridge, X_test, type = "response"))
      ridge_test_auc  <- auc(roc(y_test, ridge_test_pred))
      selected_params <- rbind(selected_params, data.frame(
        Scheme = scheme, Split = i,
        Train_Range = paste(train_codes[1], tail(train_codes,1), sep = "-"),
        Test_Set = test_code, Model = "Ridge",
        Param = "lambda", Value = best_lambda_ridge
      ))
      
      # Lasso: alpha = 1
      lasso_full <- glmnet(X_full, y_full, alpha = 1, family = "binomial")
      lam_lasso <- lasso_full$lambda
      agg_auc_lasso <- numeric(length(lam_lasso))
      for (f in seq_along(folds)) {
        df_tr <- df %>% filter(vaccine_code %in% folds[[f]]$train) %>% mutate(Dom = factor(Dom, levels = c("0","1")))
        df_va <- df %>% filter(vaccine_code == folds[[f]]$val)   %>% mutate(Dom = factor(Dom, levels = c("0","1")))
        X_tr  <- as.matrix(df_tr %>% select(all_of(predictors)))
        y_tr  <- df_tr$Dom
        X_va  <- as.matrix(df_va %>% select(all_of(predictors)))
        y_va  <- df_va$Dom
        mod_tr_l <- glmnet(X_tr, y_tr, alpha = 1, family = "binomial", lambda = lam_lasso)
        preds_mat_l <- predict(mod_tr_l, X_va, type = "response")
        for (j in seq_along(lam_lasso)) {
          agg_auc_lasso[j] <- agg_auc_lasso[j] + auc(roc(y_va, preds_mat_l[, j]))
        }
      }
      mean_auc_lasso <- agg_auc_lasso / length(folds)
      best_idx_lasso <- which.max(mean_auc_lasso)
      best_lambda_lasso <- lam_lasso[best_idx_lasso]
      final_lasso <- glmnet(X_full, y_full, alpha = 1, family = "binomial", lambda = best_lambda_lasso)
      lasso_train_pred <- as.vector(predict(final_lasso, X_full, type = "response"))
      lasso_train_auc  <- auc(roc(y_full, lasso_train_pred))
      lasso_test_pred <- as.vector(predict(final_lasso, X_test, type = "response"))
      lasso_test_auc  <- auc(roc(y_test, lasso_test_pred))
      selected_params <- rbind(selected_params, data.frame(
        Scheme = scheme, Split = i,
        Train_Range = paste(train_codes[1], tail(train_codes,1), sep = "-"),
        Test_Set = test_code, Model = "Lasso",
        Param = "lambda", Value = best_lambda_lasso
      ))
    }
    
    
    #### 3. SVM
    svm_grid <- expand.grid(
      C     = c(0.1, 1, 10, 100),
      gamma = c(1/ncol(X_full), 0.01, 0.1, 1)
    )
    
    
    agg_auc <- numeric(nrow(svm_grid))
    for (f in seq_along(folds)) {
      df_tr <- df %>% filter(vaccine_code %in% folds[[f]]$train) %>% mutate(Dom = factor(Dom, levels = c("0","1")))
      df_va <- df %>% filter(vaccine_code == folds[[f]]$val)   %>% mutate(Dom = factor(Dom, levels = c("0","1")))
      X_tr  <- as.matrix(df_tr %>% select(all_of(predictors)))
      y_tr  <- df_tr$Dom
      X_va  <- as.matrix(df_va %>% select(all_of(predictors)))
      y_va  <- df_va$Dom
      for (g in seq_len(nrow(svm_grid))) {
        mod <- fit_model("SVM", X_tr, y_tr, svm_grid[g, ])
        preds <- predict_prob(mod, X_va, "SVM")
        agg_auc[g] <- agg_auc[g] + auc(roc(y_va, preds))
      }
    }
    mean_auc <- agg_auc / length(folds)
    best_idx <- which.max(mean_auc)
    best_svm <- svm_grid[best_idx, ]
    final_svm <- fit_model("SVM", X_full, y_full, best_svm)
    svm_train_pred <- predict_prob(final_svm, X_full, "SVM")
    svm_train_auc  <- auc(roc(y_full, svm_train_pred))
    svm_test_pred <- predict_prob(final_svm, X_test, "SVM")
    svm_test_auc  <- auc(roc(y_test, svm_test_pred))
    selected_params <- rbind(selected_params, data.frame(
      Scheme = scheme, Split = i,
      Train_Range = paste(train_codes[1], tail(train_codes,1), sep="-"),
      Test_Set = test_code, Model="SVM",
      Param = c("C","gamma"), Value = c(best_svm$C, best_svm$gamma)
    ))
    
    #### 4. RF
    rf_grid <- expand.grid(
      mtry  = c(floor(sqrt(ncol(X_full))), ceiling(ncol(X_full)/3), ncol(X_full)),
      ntree = c(300, 500, 1000)
    )
    agg_auc_rf <- numeric(nrow(rf_grid))
    for (f in seq_along(folds)) {
      df_tr <- df %>% filter(vaccine_code %in% folds[[f]]$train) %>% mutate(Dom = factor(Dom, levels = c("0","1")))
      df_va <- df %>% filter(vaccine_code == folds[[f]]$val)   %>% mutate(Dom = factor(Dom, levels = c("0","1")))
      X_tr  <- as.matrix(df_tr %>% select(all_of(predictors)))
      y_tr  <- df_tr$Dom
      X_va  <- as.matrix(df_va %>% select(all_of(predictors)))
      y_va  <- df_va$Dom
      for (g in seq_len(nrow(rf_grid))) {
        mod <- fit_model("RF", X_tr, y_tr, rf_grid[g, ])
        preds <- predict_prob(mod, X_va, "RF")
        agg_auc_rf[g] <- agg_auc_rf[g] + auc(roc(y_va, preds))
      }
    }
    mean_auc_rf <- agg_auc_rf / length(folds)
    best_idx_rf <- which.max(mean_auc_rf)
    best_rf <- rf_grid[best_idx_rf, ]
    final_rf <- fit_model("RF", X_full, y_full, best_rf)
    
    df_rf_imp <- data.frame(final_rf$importance)
    
    rf_train_pred <- predict_prob(final_rf, X_full, "RF")
    rf_train_auc  <- auc(roc(y_full, rf_train_pred))
    rf_test_pred <- predict_prob(final_rf, X_test, "RF")
    rf_test_auc  <- auc(roc(y_test, rf_test_pred))
    selected_params <- rbind(selected_params, data.frame(
      Scheme=scheme, Split=i,
      Train_Range=paste(train_codes[1], tail(train_codes,1), sep="-"),
      Test_Set=test_code, Model="RF",
      Param=c("mtry","ntree"), Value=c(best_rf$mtry, best_rf$ntree)
    ))
    analysis_schemes
    
    #### 5. DNN
    n_predictors <- length(predictors)
    dnn_grid <- expand.grid(
      units      = c(n_predictors, n_predictors*2),  
      layers     = c(1, 2),                         
      dropout    = c(0.2, 0.5),                    
      lr         = c(1e-3, 5e-4),                 
      batch_size = c(32),                             
      epochs     = c(100),                           
      optimizer  = c("adam")                          
    )
    
    agg_auc_dnn <- numeric(nrow(dnn_grid))
    
    auc_metric <- tf$keras$metrics$AUC(name = "auc")
    es_cb <- callback_early_stopping(
      monitor              = "val_auc",
      mode                 = "max",
      patience             = 5,
      restore_best_weights = TRUE
    )
    
    for (f in seq_along(folds)) {
      df_tr <- df %>% filter(vaccine_code %in% folds[[f]]$train) %>% mutate(Dom = factor(Dom, levels = c("0","1")))
      df_va <- df %>% filter(vaccine_code ==   folds[[f]]$val)   %>% mutate(Dom = factor(Dom, levels = c("0","1")))
      X_tr  <- as.matrix(df_tr %>% select(all_of(predictors)))
      y_tr  <- as.numeric(as.character(df_tr$Dom))
      X_va  <- as.matrix(df_va %>% select(all_of(predictors)))
      y_va  <- as.numeric(as.character(df_va$Dom))
      
      for (g in seq_len(nrow(dnn_grid))) {
        params <- dnn_grid[g, ]
        model <- keras_model_sequential() %>%
          layer_dense(units = params$units, activation = "relu", input_shape = ncol(X_tr)) %>%
          layer_dropout(rate = params$dropout)
        if (params$layers == 2) {
          model <- model %>%
            layer_dense(units = params$units %/% 2, activation = "relu") %>%
            layer_dropout(rate = params$dropout)
        }
        model <- model %>%
          layer_dense(units = 1, activation = "sigmoid") %>%
          compile(
            optimizer = optimizer_adam(lr = params$lr),
            loss      = "binary_crossentropy",
            metrics   = list(auc_metric)
          )
        
        history <- model %>% fit(
          x               = X_tr, y = y_tr,
          validation_data = list(X_va, y_va),
          batch_size      = params$batch_size,
          epochs          = params$epochs,
          callbacks       = list(es_cb),
          verbose         = 0
        )
        
        best_val_auc <- max(history$metrics$val_auc)
        agg_auc_dnn[g] <- agg_auc_dnn[g] + best_val_auc
        
        k_clear_session()
      }
    }
    
    mean_auc_dnn <- agg_auc_dnn / length(folds)
    best_idx    <- which.max(mean_auc_dnn)
    best_dnn    <- dnn_grid[best_idx, ]

    model_final <- keras_model_sequential() %>%
      layer_dense(units = best_dnn$units, activation = "relu", input_shape = ncol(X_full)) %>%
      layer_dropout(rate = best_dnn$dropout)
    if (best_dnn$layers == 2) {
      model_final <- model_final %>%
        layer_dense(units = best_dnn$units %/% 2, activation = "relu") %>%
        layer_dropout(rate = best_dnn$dropout)
    }
    
    model_final <- model_final %>%
      layer_dense(units = 1, activation = "sigmoid") %>%
      compile(
        optimizer = optimizer_adam(lr = best_dnn$lr),
        loss      = "binary_crossentropy",
        metrics   = list(tf$keras$metrics$AUC(name = "val_auc"))
      )
    

    model_final %>% fit(
      x                = X_full,
      y                = as.numeric(as.character(y_full)),
      validation_split = 0.1,       
      batch_size       = best_dnn$batch_size,
      epochs           = best_dnn$epochs,
      callbacks        = list(es_cb),# EarlyStopping
      verbose          = 0
    )
    
    dnn_test_pred <- model_final %>% predict(X_test)
    dnn_train_pred<- model_final %>% predict(X_full)
    dnn_train_auc <- auc(roc(as.numeric(as.character(y_full)), dnn_train_pred))
    dnn_test_auc  <- auc(roc(as.numeric(as.character(y_test)), dnn_test_pred))
    
    selected_params <- rbind(
      selected_params,
      data.frame(
        Scheme      = scheme, Split = i,
        Train_Range = paste(train_codes[1], tail(train_codes,1), sep="-"),
        Test_Set    = test_code, Model = "DNN",
        Param       = c("units","layers","dropout","lr","batch_size","epochs", "optimizer"),
        Value       = c(unlist(best_dnn)[1:6], best_dnn[1,7] %>% as.character())
      )
    )
    
    #### 6. GLMM
    set.seed(777)
    form_str <- paste("Dom ~", paste(predictors, collapse=" + "), "+ (1 | vaccine_code)")
    glmm_model <- glmer(as.formula(form_str), data = train_df, family = "binomial")
    glmm_train_pred <- predict(glmm_model, newdata = train_df, type = "response")
    glmm_test_pred  <- predict(glmm_model, newdata = test_df,  type = "response", allow.new.levels = TRUE)
    glmm_train_auc  <- auc(roc(y_full, glmm_train_pred))
    glmm_test_auc   <- auc(roc(y_test,  glmm_test_pred))
    

    res[i, "Logistic_Train_AUROC"] <- ifelse(inherits(logit_model, "try-error"), NA, logit_train_auc)
    res[i, "Logistic_Test_AUROC"]  <- ifelse(inherits(logit_model, "try-error"), NA, logit_test_auc)
    
    res[i, "Ridge_Train_AUROC"] <- ifelse(inherits(final_ridge, "try-error"), NA, ridge_train_auc)
    res[i, "Ridge_Test_AUROC"]  <- ifelse(inherits(final_ridge, "try-error"), NA, ridge_test_auc)
    
    res[i, "Lasso_Train_AUROC"] <- ifelse(inherits(final_lasso, "try-error"), NA, lasso_train_auc)
    res[i, "Lasso_Test_AUROC"]  <- ifelse(inherits(final_lasso, "try-error"), NA, lasso_test_auc)
    
    res[i, "SVM_Train_AUROC"] <- ifelse(inherits(final_svm, "try-error"), NA, svm_train_auc)
    res[i, "SVM_Test_AUROC"]  <- ifelse(inherits(final_svm, "try-error"), NA, svm_test_auc)
    
    res[i, "RF_Train_AUROC"] <- ifelse(inherits(final_rf, "try-error"), NA, rf_train_auc)
    res[i, "RF_Test_AUROC"]  <- ifelse(inherits(final_rf, "try-error"), NA, rf_test_auc)
    
    # res 테이블에도 AUROC 추가
    res[i, "DNN_Train_AUROC"] <- dnn_train_auc
    res[i, "DNN_Test_AUROC"]  <- dnn_test_auc
    
    res[i, "GLMM_Train_AUROC"] <- ifelse(inherits(glmm_model, "try-error"), NA, glmm_train_auc)
    res[i, "GLMM_Test_AUROC"]  <- ifelse(inherits(glmm_model, "try-error"), NA, glmm_test_auc)
  }  # end for each split
  
  if(length(predictors) == 1){
    res <- res %>% select(-c(Ridge_Train_AUROC, Ridge_Test_AUROC, Lasso_Train_AUROC, Lasso_Test_AUROC))
    scheme_model_order <- c("GLMM", "Logistic", "SVM", "RF", "DNN")
  } else {
    scheme_model_order <- default_model_order
  }
  
  results_all[[scheme]] <- res
  
  results_long <- res %>%
    pivot_longer(cols = -c(Train_Range, Test_Set), names_to = "Model_Type", values_to = "AUROC") %>%
    mutate(
      Type = ifelse(grepl("Train", Model_Type), "Train", "Test"),
      Model = gsub("_Train_AUROC|_Test_AUROC", "", Model_Type)
    ) %>% 
    mutate(Test_Set = factor(Test_Set, levels = res$Test_Set))
  results_long$Model <- factor(results_long$Model, levels = scheme_model_order)
  results_long$Type <- factor(results_long$Type, levels = c("Train", "Test"))
  
  p <- ggplot(results_long, aes(x = Model, y = AUROC, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), alpha = 0.8) +
    scale_fill_manual(values = c("Train" = "cornflowerblue", "Test" = "tomato")) +
    scale_y_continuous(limits = c(0,1)) +
    geom_text_repel(aes(label = round(AUROC, 2)),
                    position = position_dodge(width = 0.9),
                    vjust = -0.3, size = 3, color = "black") +
    facet_grid(rows = vars(Test_Set)) +
    labs(title = paste("Model Performance Comparison (AUROC)"),
         subtitle = analysis_name[[scheme]],
         x = "Model", y = "AUROC") +
    theme_minimal() +
    theme(
      #axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank(),
      plot.background = element_rect(fill = "white", linetype = 0)
    )
  p
  plots_all[[scheme]] <- p
  ggsave(plot = p, filename = paste0("Model_performance_AUROC_", scheme, ".png"),
         path = paste0(path, "Result/"), dpi = 300, width = 8, height = 6)
  
  res$Test_Set
  
  res_ch <- res %>%
    pivot_longer(
      cols      = -c(Train_Range, Test_Set),
      names_to  = c("Model", "Type"),
      names_pattern = "(.*)_(Train|Test)_AUROC",
      values_to = "AUROC"
    ) %>%
    pivot_wider(
      id_cols     = c(Train_Range, Test_Set, Model),
      names_from  = Type,
      values_from = AUROC
    ) %>% 
    mutate(Model = factor(Model, levels = scheme_model_order),
           Test_Set = factor(Test_Set, levels = c("HK15", "Kan17", "HK19", "Cam20"))) 
  
  scheme_model_order_ch <- gsub("GLMM", "GLMM (Baseline)", scheme_model_order)
  

  res_fin <- res_ch %>% mutate(Model = ifelse(as.character(Model) == "GLMM", "GLMM (Baseline)", as.character(Model))) %>% 
    mutate(Model = factor(Model, levels = scheme_model_order_ch))

  write.csv(res_fin, paste0(path, "Result/Summarized_results_", scheme,".csv"), row.names = F)
  
  results_long_ch <- cbind(results_long, scheme)
  if(is.null(results_long_merge)){
    results_long_merge <- results_long_ch
  }else{
    results_long_merge <- rbind(results_long_merge, results_long_ch)
  }
}

write.csv(selected_params, file = paste0(path, "Result/selected_hyperparams_nestedCV.csv"), row.names = FALSE)
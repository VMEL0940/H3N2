library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)


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

default_model_order <- c("GLMM",  "Ridge", "Lasso", "SVM", "RF", "DNN")

col_models <- c("snow4",  "cornflowerblue", "orange", "mediumvioletred", "darkgreen", "midnightblue")
names(col_models) <- default_model_order

results_all <- list()
plots_all <- list()

results_long_merge <- NULL


library(reshape2)


df_long <- df %>% melt(., id.vars = setdiff(names(df), c(Syn_vars, Nonsyn_vars))) %>%
  mutate(vaccine_code = factor(vaccine_code, levels = unique(as.character(vaccine_code))),
         Dom = factor(Dom, levels = c(0, 1))) %>%
  mutate(variable = factor(variable, levels = c(Syn_vars, Nonsyn_vars)))

df_long %>% filter(variable == "HA_Nonsyn" & vaccine_code == "Cam20")

df_long_ch <- df_long %>% mutate(var_type = ifelse(variable %in% Syn_vars, "Synonymous", "Non-Synonymous")) %>%
  mutate(var_type = factor(var_type, levels = c("Synonymous", "Non-Synonymous"))) %>%
  mutate(variable = gsub("_Nonsyn", "", gsub("_Syn", "", as.character(variable)))) %>%
  mutate(variable = factor(variable, levels = unique(as.character(variable))))

################# Plot regeneration 
for (scheme in names(analysis_schemes)) {
  #scheme <- names(analysis_schemes)[1]
  
  predictors <- analysis_schemes[[scheme]]
  
  # HA_only (univariate)인 경우 Ridge와 Lasso 결과 컬럼 제거 및 모델 순서 재정의
  if(length(predictors) == 1){
    scheme_model_order <- c("GLMM", "SVM", "RF", "DNN")
  } else {
    scheme_model_order <- default_model_order
  }
  
  scheme_model_order_ch <- gsub("GLMM", "GLMM (Baseline)", scheme_model_order)
  
  
  res_fin <- read.csv(paste0(path, "Result/Summarized_results_", scheme,".csv")) %>% 
    filter(Model != "Logistic") %>%
    mutate(Model = factor(Model, levels = scheme_model_order_ch), 
           Test_Set = factor(Test_Set, levels = c("HK15", "Kan17", "HK19", "Cam20"))) 
  better_models <- res_fin %>% group_by(Test_Set) %>% filter(Test > Test[Model == "GLMM (Baseline)"])
  
  col_models_sub <- col_models[scheme_model_order] 
  
  names(col_models_sub)[names(col_models_sub) == "GLMM"] <- "GLMM (Baseline)"
  
  shape_models <- c(4, rep(16, length(scheme_model_order_ch)-1))
  names(shape_models) <- scheme_model_order_ch
  
  p_1 <- ggplot(res_fin, aes(x = Train, y = Test, col = Model, shape = Model)) +
    geom_point(alpha = 0.8) +
    scale_shape_manual(values = shape_models)+
    scale_color_manual(values = col_models_sub) +
    geom_abline(slope = 1, intercept = 0, col = "grey", alpha = 0.9, linetype = "dashed") +
    scale_x_continuous(limits = c(0,1)) +
    scale_y_continuous(limits = c(0,1)) +
    facet_wrap(facets = vars(Test_Set), nrow = 2) +
    labs(
      x = "Train AUROC", y = "Test AUROC") +
    theme_bw() +
    theme(
      #axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank(),
      plot.background = element_rect(fill = "white", linetype = 0)
    )
  p_1
  #plots_all[[scheme]] <- p
  ggsave(plot = p_1, filename = paste0("Fig_6_B_", scheme, ".svg"),
         path = paste0(path, "Result/"), dpi = 600, width = 6, height = 4.5)
}


df_res_tot <- NULL
for(scheme in names(analysis_schemes)){
  print(scheme)
  df_i <- read.csv(paste0(path, "Result/Summarized_results_",scheme,".csv"))
  df_i %>% group_by(Test_Set) %>% filter(Train > Train[Model == "GLMM (Baseline)"] & Test > Test[Model == "GLMM (Baseline)"]) %>% print()
  df_res_sub <- data.frame(df_i, analysis_scheme = scheme)
  if(is.null(df_res_tot)){
    df_res_tot <- df_res_sub
  }else{
    df_res_tot <- rbind(df_res_tot, df_res_sub)
  }
  
}

df_res_tot <- NULL
for(scheme in names(analysis_schemes)){
  print(scheme)
  df_i <- read.csv(paste0(path, "Result/Summarized_results_",scheme,".csv"))
  df_i %>% group_by(Test_Set) %>% filter(Train > Train[Model == "GLMM (Baseline)"] & Test > Test[Model == "GLMM (Baseline)"]) %>% print()
  df_res_sub <- data.frame(df_i, analysis_scheme = scheme)
  if(is.null(df_res_tot)){
    df_res_tot <- df_res_sub
  }else{
    df_res_tot <- rbind(df_res_tot, df_res_sub)
  }
  
}

for(ts_s in c("HK15", "Kan17", "HK19", "Cam20")){
  
  df_res_tot_sub <- df_res_tot %>% filter(Test_Set == ts_s) %>% 
    filter(Model !="Logistic") %>% 
    mutate(Model = factor(Model, levels = scheme_model_order_ch),
           analysis_scheme = factor(analysis_scheme, c("HA_only","Original_var","Nonsyn_var","Full_var")))
  
  
  df_res_tot_long <- df_res_tot %>% filter(Test_Set == ts_s) %>% 
    melt(., id.vars = setdiff(names(df_res_tot), c("Train", "Test"))) %>% 
    mutate(Model = factor(Model, levels = scheme_model_order_ch),
           analysis_scheme = factor(analysis_scheme, c("HA_only","Original_var","Nonsyn_var","Full_var")))
  
  
  df_res_dummy <- df_res_tot_long[1:4,]
  df_res_dummy$analysis_scheme <- "HA_only"
  df_res_dummy$Model <- rep(c("Ridge", "Lasso"), each = 2)
  df_res_dummy$variable <- rep(c("Train", "Test"), 2)
  df_res_dummy$value <- 0
  
  df_res_fin <- rbind(df_res_tot_long, df_res_dummy)
  
  # p_res <- ggplot(df_res_fin, aes(x = Model, y = value, fill = analysis_scheme)) +
  #   geom_col(position = position_dodge2()) +
  #   facet_wrap(vars(variable), nrow = 2) +
  #   theme_bw()
  # 
  # ggsave(plot = p_res, filename = paste0("res_summary_bar_",ts_s,".png"),
  #        path = paste0(path, "Result/"), dpi = 300, width = 12, height = 8)
  
  
  # 1) Overfit best: Train > Test 중에서 Test max
  best_overfit <- df_res_tot_sub %>%
    filter(Train > Test) %>%
    slice_max(Test, n = 1, with_ties = TRUE) %>%
    ungroup() %>%
    mutate(Criterion = "Criteria 1")
  
  # 2) Overall best: 전체 중에서 Test max
  best_overall <- df_res_tot_sub %>%
    slice_max(Test, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Criterion = "Criteria 2")
  
  # 3) 두 테이블 합치기
  best_both <- bind_rows(best_overfit, best_overall)
  
  # overlap_labels <-  df_res_tot_sub %>%
  #   # ① 반올림용 컬럼 추가 (소수점 3자리까지)
  #   mutate(
  #     Train_r = round(Train, 3),
  #     Test_r  = round(Test,  3)
  #   ) %>%
  #   # ② 반올림된 값으로 그룹핑
  #   group_by(Test_Set, Train_r, Test_r) %>%
  #   filter(n() > 1) %>%
  #   summarise(
  #     # 실제 레이블에는 원래 값(또는 둘 다 보여줄 수도) 사용
  #     Label = paste(Model, collapse = ", "),
  #     # 레이블 위치로는 반올림된 좌표를 씁니다
  #     Train = first(Train_r),
  #     Test  = first(Test_r),
  #     .groups = "drop"
  #   )
  #print(best_overfit)
  
  p_res_point <- ggplot(df_res_tot_sub, aes(x = Train, y = Test, col = Model, size = analysis_scheme)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = col_models_sub) +
    #scale_size_manual(name = "Predictor set") +
    labs(size = "Predictor set")+
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))+
    # geom_text_repel(
    #   data         = best_both,
    #   aes(x = Train, y = Test, label = paste0(format(Train, digits = 3), ", ", format(Test, digits = 3))),
    #   inherit.aes  = FALSE,
    #   arrow        = arrow(length = unit(0.02, "inches")),
    #   box.padding  = 0.5,
    #   point.padding= 0.3,
    #   segment.size = 0.6,
    #   segment.color= "grey40",
    #   min.segment.length = unit(0, "lines"),
    #   show.legend  = FALSE
    # ) +
    theme_bw()
  
  p_res_point
  
  ggsave(plot = p_res_point, filename = paste0("Fig_6_C_",ts_s,".svg"),
         path = paste0(path, "Result/"), dpi = 600, width = 8.5, height = 7)
  
  #write.csv(df_res_tot_sub %>% arrange(-Test, -Train) %>% select(-Train_Range), paste0(path, "Result/res_table_by_periods_",ts_s,".csv"),row.names = F)
}

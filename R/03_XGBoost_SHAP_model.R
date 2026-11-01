#XGBOOST-SHAP MODEL FRAMEWORK

#Load libraries
pacman::p_load(tidyverse, lubridate, xgboost, reticulate, corrplot, shapviz, dplyr,
               patchwork, knitr, kableExtra, ggplot2, viridis)


#Load daily methane dataset created using the script 02_script_to_create_comp_dataset.R
set.seed(42)

comp_data <- read_csv(
  "C:/Users/bibek/Desktop/Bibek_files/5_yrs_FCR_flux_analysis/comp_data.csv",
  show_col_types = FALSE
) |>
  mutate(Date = date) |>
  arrange(Date) |>
  mutate(
    mean_Chla_lag1 = dplyr::lag(mean_Chla, n = 1),
    DOY            = lubridate::yday(Date)
  )

comp_data <- comp_data %>% 
  mutate(daily_CH4 = daily_CH4 * 1386.3) #change scale to mg/m2/d from umol/m2/s
####

#Calculate water stability index (Relative Water Stability Index following Marcon et al., 2024)
# UNESCO freshwater density equation (kg/m3)
water_density <- function(T) {
  1000 * (1 -
            ((T + 288.9414) /
               (508929.2 * (T + 68.12963))) *
            (T - 3.9863)^2)
}

# Reference density difference
rho_4 <- water_density(4)
rho_5 <- water_density(5)

rho_reference <- rho_4 - rho_5

comp_data <- comp_data %>%
  mutate(
    rho_surface = water_density(temp_1.6m),
    rho_bottom = water_density(temp_9m),
    
    RWCS = (rho_bottom - rho_surface) / rho_reference,
    
    stratification = case_when(
      RWCS > 56.5 ~ "Stratified",
      RWCS < 16.3 ~ "Mixed",
      TRUE ~ "Partially stratified"
    )
  )


#Read in water chemistry dataset
carbon <- read_csv("https://pasta.lternet.edu/package/data/eml/edi/199/13/3f09a3d23b7b5dd32ed7d28e9bc1b081?key=yltMpS4UEIk12AvB9L7OL5uRiG0") |> 
  mutate(DateTime = as.POSIXct(strptime(DateTime, "%Y-%m-%d %H:%M", tz="EST"))) |>  
  filter(Reservoir == "FCR" & Site == 50 & Depth_m == 0.1) |> 
  filter(DateTime >= "2020-04-04") |>  
  mutate(DateTime = round_date(DateTime, "30 mins"),
         date = as.Date(DateTime)) |> 
  select(date, DC_mgL, TP_ugL)

comp_data <- comp_data |> 
  left_join(carbon, by = "date")

################################################################################
#Check autocorrelation and lags
AR_df <- comp_data |>
  dplyr::select(-date, -Date)

# remove NA values
chla_series <- na.omit(AR_df$mean_Chla)

# PACF plot
pacf(
  chla_series,
  lag.max = 30,   # adjust as needed
  main    = "",
  col     = "#2166AC",
  lwd     = 2
)


###########
#Check correlation
#prepare data 
cor_data <- comp_data %>%
  select(temp_1.6m, temp_5m, temp_9m, DOsat_surface, DOsat_9m, DOsat_5m, mean_fDOM, mean_Chla_lag1,
         WindSpeed_Average_m_s, AirTemp_C_Average, DOY, daily_CH4, sch_stab, thermo.depth, TP_ugL)

#correlation matrix
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs")

#cleaner plot
corrplot(
  cor_matrix,
  method   = "color",
  type     = "upper",
  order    = "hclust",
  
  # appearance
  col      = colorRampPalette(
    c("#3B4CC0", "white", "#B40426")
  )(80),
  
  # labels
  tl.col   = "black",
  tl.cex   = 0.6,
  tl.srt   = 35,
  
  # show correlation coefficients
  diag        = FALSE,
  addCoef.col = "black",
  number.cex  = 0.6,     # size of the coefficient text
  number.digits = 2,      # decimal places shown
  
  # remove visual clutter
  cl.cex   = 0.7,
  mar      = c(0, 0, 1, 0),
  
  title    = "Correlation Matrix"
)
############

#SHARED SETTINGS
xvars_clean <- c("WindSpeed_Average_m_s", "AirTemp_C_Average", "DOsat_5m",
                 "mean_Chla_lag1", "temp_5m", "DOY", "thermo.depth", "RWCS", "TP_ugL")


##################
#Creating virtual Python environment for reticulate to work in R
# install_python(version = "3.12")
# virtualenv_create("r-reticulate")
# use_virtualenv("r-reticulate", required = TRUE)
# 
# py_config()

#################

use_virtualenv("r-reticulate", required = TRUE)
optuna <- reticulate::import("optuna")
optuna$logging$set_verbosity(optuna$logging$WARNING)

print_progress <- function(current, total, best_rmse, width = 40) {
  pct    <- current / total
  filled <- round(pct * width)
  bar    <- paste0(strrep("=", filled), strrep("-", width - filled))
  cat(sprintf("\r  [%s] %3.0f%% | Trial %d/%d | Best RMSE: %.7f",
              bar, pct * 100, current, total, best_rmse))
  if (current == total) cat("\n")
  flush.console()
}

rmse_fn <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

run_optuna <- function(train_x, train_y, n_trials = 50L) {
  best_rmse    <<- Inf
  best_nrounds <<- 100L
  trial_count  <<- 0L
  
  objective <- function(trial) {
    trainMatrix_cv <- xgb.DMatrix(data = train_x, label = train_y)
    
    params <- list(
      booster          = "gbtree",
      objective        = "reg:squarederror",
      eval_metric      = "rmse",
      eta              = trial$suggest_float("eta",              0.01,  0.1),
      max_depth        = trial$suggest_int("max_depth",          1L,    7L),
      min_child_weight = trial$suggest_int("min_child_weight",   1L,   15L),
      subsample        = trial$suggest_float("subsample",        0.5,   0.9),
      colsample_bytree = trial$suggest_float("colsample_bytree", 0.4,   0.9),
      lambda           = trial$suggest_float("lambda",           1e-3, 15, log = TRUE),
      alpha            = trial$suggest_float("alpha",            1e-3, 20, log = TRUE),
      seed             = 42L
    )

    cv_result <- xgb.cv(
      params                = params,
      data                  = trainMatrix_cv,
      nrounds               = 100,
      nfold                 = 5,
      verbose               = 0,
      early_stopping_rounds = 30
    )
    
    trial_rmse <- min(cv_result$evaluation_log$test_rmse_mean)
    
    if (trial_rmse < best_rmse) {
      best_rmse    <<- trial_rmse
      best_nrounds <<- cv_result$best_iteration
    }
    
    trial_count <<- trial_count + 1L
    print_progress(trial_count, n_trials, best_rmse)
    return(trial_rmse)
  }
  
  sampler <- optuna$samplers$TPESampler(seed = 42L)
  study   <- optuna$create_study(direction = "minimize", sampler = sampler)
  study$optimize(objective, n_trials = n_trials)
  
  list(study = study, best_nrounds = best_nrounds)
}

fit_xgb <- function(train_x, train_y, test_x, test_y, best_params,
                    nrounds = 500L, verbose = 0) {
  trainMatrix <- xgb.DMatrix(data = train_x, label = train_y)
  testMatrix  <- xgb.DMatrix(data = test_x,  label = test_y)
  
  fit <- xgb.train(
    params = list(
      booster          = "gbtree",
      objective        = "reg:squarederror",
      eta              = best_params$eta,
      max_depth        = as.integer(best_params$max_depth),
      min_child_weight = as.integer(best_params$min_child_weight),
      subsample        = best_params$subsample,
      colsample_bytree = best_params$colsample_bytree,
      lambda           = best_params$lambda,
      alpha            = best_params$alpha,
      seed             = 42L
    ),
    data                  = trainMatrix,
    nrounds               = nrounds,
    evals                 = list(train = trainMatrix, test = testMatrix),
    early_stopping_rounds = 20L,
    verbose               = verbose
  )
  list(fit = fit, trainMatrix = trainMatrix, testMatrix = testMatrix)
}

eval_performance <- function(actual_train, pred_train,
                             actual_test,  pred_test,
                             model_label,
                             transform = c("none", "log1p")) {
  transform <- match.arg(transform)
  if (transform == "log1p") {
    pred_train <- expm1(pred_train)
    pred_test  <- expm1(pred_test)
  }
  
  tr_rmse <- rmse_fn(actual_train, pred_train)
  te_rmse <- rmse_fn(actual_test,  pred_test)
  
  ss_res_tr <- sum((actual_train - pred_train)^2,                       na.rm = TRUE)
  ss_tot_tr <- sum((actual_train - mean(actual_train, na.rm = TRUE))^2, na.rm = TRUE)
  r2_train  <- 1 - ss_res_tr / ss_tot_tr
  
  ss_res    <- sum((actual_test - pred_test)^2,                       na.rm = TRUE)
  ss_tot    <- sum((actual_test - mean(actual_test, na.rm = TRUE))^2, na.rm = TRUE)
  r2        <- 1 - ss_res / ss_tot
  
  cat(sprintf("\n── %s Performance ──────────────────\n", model_label))
  cat("Train RMSE:", round(tr_rmse, 6), "\n")
  cat("Test RMSE: ", round(te_rmse, 6), "\n")
  cat("Train R²:  ", round(r2_train, 3), "\n")
  cat("Test R²:   ", round(r2,       3), "\n")
  
  data.frame(
    Metric = c("Train R²", "Test R²", "Train RMSE", "Test RMSE"),
    Value  = c(round(r2_train, 3),
               round(r2,      3),
               round(tr_rmse, 6),
               round(te_rmse, 6))
  )
}

impute_from_medians <- function(df, medians) {
  df |>
    left_join(medians, by = "month", suffix = c("", "_med")) |>
    mutate(across(all_of(xvars_clean),
                  ~ coalesce(.x, get(paste0(cur_column(), "_med"))))) |>
    dplyr::select(-ends_with("_med"))
}

#XGBoost-SHAP driver analysis using non gap-filled raw data

cat("\n", strrep("═", 60), "\n")
cat("  MODEL A: Raw CH4 → SHAP Analysis\n")
cat(strrep("═", 60), "\n")

model_data_raw <- comp_data |>
  mutate(month = lubridate::month(date)) |>
  arrange(date) |>  # make sure it's sorted
  dplyr::select(all_of(c("daily_CH4", xvars_clean, "month"))) |>
  dplyr::filter(!is.na(daily_CH4))



cat("Model A training rows (raw CH4):", nrow(model_data_raw), "\n")

set.seed(42)
signed_log1p  <- function(x) sign(x) * log1p(abs(x))
signed_expm1  <- function(x) sign(x) * expm1(abs(x))

train_idx_A  <- sample(nrow(model_data_raw), 0.75 * nrow(model_data_raw))
trainData_A  <- model_data_raw[train_idx_A, ]
testData_A   <- model_data_raw[-train_idx_A, ]

train_med_A  <- trainData_A |>
  group_by(month) |>
  summarise(across(all_of(xvars_clean), ~ median(.x, na.rm = TRUE)),
            .groups = "drop")

trainData_A  <- impute_from_medians(trainData_A, train_med_A)
testData_A   <- impute_from_medians(testData_A,  train_med_A)

train_x_A    <- as.matrix(trainData_A[, xvars_clean])
train_y_A    <- (trainData_A$daily_CH4)
test_x_A     <- as.matrix(testData_A[,  xvars_clean])
test_y_A     <- (testData_A$daily_CH4)

cat("\nTuning Model A...\n")
tuneA        <- run_optuna(train_x_A, train_y_A, n_trials = 50L)
fitA_out     <- fit_xgb(train_x_A, train_y_A, test_x_A, test_y_A,
                        tuneA$study$best_params, verbose = 0)
fit_raw      <- fitA_out$fit

perf_A <- eval_performance(
  actual_train = trainData_A$daily_CH4,
  pred_train   = (predict(fit_raw, fitA_out$trainMatrix)),
  actual_test  = testData_A$daily_CH4,
  pred_test    = (predict(fit_raw, fitA_out$testMatrix)),
  model_label  = "Model A — Raw CH4 (SHAP)",
  transform    = "none"                  
)

#SHAP model
set.seed(42)
X_explain_A <- model_data_raw[sample(nrow(model_data_raw),
                                     min(415, nrow(model_data_raw))),
                              xvars_clean]

shp_A <- shapviz(fit_raw,
                 X_pred = data.matrix(X_explain_A),
                 X      = X_explain_A)

#rename variables for figures
shap_labels <- c(
  "DOsat_5m"                  = "Saturated DO (5 m)",
  "AirTemp_C_Average"         = "Air temperature",
  "temp_5m"                   = "Temperature (5 m)",
  "DOsat_9m"                  = "Saturated DO (9 m)",
  "temp_9m"                   = "Temperature (9 m)",
  "WindSpeed_Average_m_s"     = "Wind speed",
  "mean_Chla_lag1"            = "Chla lag1"
)


p_imp_A <- sv_importance(shp_A, show_numbers = TRUE) + scale_y_discrete(
  labels = shap_labels) +
  labs(title = "") +
  theme_classic()

p_bee_A <- sv_importance(shp_A, kind = "beeswarm") + scale_y_discrete(
  labels = shap_labels) +
  labs(title = "") +
  scale_color_gradient(low = "blue", high = "red",
                       na.value = NA, name = "Feature value") +
  labs(title = "") + theme_classic()

print(p_imp_A)
print(p_bee_A)

#Dependence plots
p1 <- sv_dependence(shp_A, v = "DOsat_5m", color_var = NULL) +
  labs(x = "% DO saturation (5 m)", y = "SHAP value") 
p2 <- sv_dependence(shp_A, v = "AirTemp_C_Average", color_var = NULL)
p3 <- sv_dependence(shp_A, v = "temp_5m", color_var = NULL) +
  labs(x = "Temperature (5 m)", y = "SHAP value")
p4 <- sv_dependence(shp_A, v = "WindSpeed_Average_m_s", color_var = NULL)

p1 + p2 + p3 + p4

# Figure 3
p1 <- p1 + theme_classic()
p3 <- p3 + theme_classic()

shap_plots <- 
  (p_imp_A + p_bee_A) /
  p1 /
  p3 +
  plot_layout(
    heights = c(1, 1.3, 1.3)
  ) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 14
      )
    )
  )

ggsave(
  "Figure_3_SHAP.png",
  shap_plots,
  width = 11,
  height = 12,
  units = "in",
  dpi = 600
)

#Interaction plot
sv_dependence(shp_A, v = "temp_5m",   color_var = "DOsat_5m")
sv_dependence(shp_A, v = "temp_9m",   color_var = "DOsat_9m")

# ── Model A performance table ─────────────────────────────────────────────────
perf_A |>
  knitr::kable(
    align   = c("l", "r"),
    caption = "Model A performance — raw daily CH4 (driver analysis)"
  ) |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE) |>
  row_spec(0, bold = TRUE)

###




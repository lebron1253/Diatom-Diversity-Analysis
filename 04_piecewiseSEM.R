# 04_piecewiseSEM.R
# Piecewise structural equation models for HFP–biodiversity pathways
# Packages: piecewiseSEM, glmmTMB, vegan, e1071, dplyr, tidyr

library(dplyr)
library(tidyr)
library(piecewiseSEM)
library(glmmTMB)
library(vegan)
library(e1071)

data_raw        <- read.table("metadata.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
data_raw$Basin  <- as.factor(data_raw$Basin)
data_raw$Type   <- as.factor(data_raw$Type)

# Water quality PCA
vars_water <- c("TN", "TP", "EC", "pH", "DO", "COD")
valid_rows <- complete.cases(data_raw[, vars_water])
pca_water  <- prcomp(data_raw[valid_rows, vars_water], scale.=TRUE)
data_raw$Water_PC1         <- NA
data_raw$Water_PC2_Organic <- NA
data_raw$Water_PC1[valid_rows]         <- pca_water$x[, 1]
data_raw$Water_PC2_Organic[valid_rows] <- pca_water$x[, 2]

vars_pred <- c("Elevation", "CatchArea", "MAT", "MAP", "Qmon",
               "Water_PC1", "Water_PC2_Organic", "HFP")
vars_div  <- c("TD", "FD", "PD", "LCBDtax", "LCBDfunc", "LCBDphylo")
all_vars  <- c(vars_pred, vars_div)

# Log transformation + Z-score
process_data <- data_raw
vars_to_check_skew <- c("CatchArea", "Elevation", "Qmon", "TD", "TD", "FD", "PD", "LCBDtax", "LCBDfunc", "LCBDphylo")
for (v in vars_to_check_skew) {
  if (!v %in% names(process_data)) next
  x <- process_data[[v]]
  sk <- tryCatch(skewness(x, na.rm=TRUE), error=function(e) NA)
  if (!is.na(sk) && abs(sk) > 1) {
    mn <- min(x, na.rm=TRUE)
    if (mn > 0) {
      process_data[[v]] <- log10(x)
    } else {
      process_data[[v]] <- log10(x + abs(mn) + 1)
    }
  }
}
process_data <- process_data[complete.cases(process_data[, all_vars]), ]
for (v in all_vars) {
  if (!v %in% names(process_data)) next
  process_data[[v]] <- as.numeric(scale(process_data[[v]], center = TRUE, scale = TRUE)[,1])
}

# Component GLMM
fit_glmm <- function(response, predictors, data) {
  form <- as.formula(paste(response, "~",
                           paste(predictors, collapse=" + "), "+ (1|Basin)"))
  glmmTMB(form, data=data, family=gaussian())
}
# piecewiseSEM
run_sem <- function(data, div_metric) {
  data <- as.data.frame(data)
  m_mat   <- fit_glmm("MAT",               c("Elevation"), data)
  m_map   <- fit_glmm("MAP",               c("Elevation"), data)
  m_hydro <- fit_glmm("Qmon",              c("CatchArea", "MAP", "HFP", "Elevation"), data)
  m_wpc1  <- fit_glmm("Water_PC1",         c("Elevation", "CatchArea", "MAT", "MAP", "Qmon", "HFP"), data)
  m_wpc2  <- fit_glmm("Water_PC2_Organic", c("Elevation", "CatchArea", "MAT", "MAP", "Qmon", "HFP"), data)
  m_div   <- fit_glmm(div_metric,          c("Elevation", "CatchArea", "MAT", "MAP", "Qmon",
                                             "Water_PC1", "Water_PC2_Organic", "HFP",
                                             "HFP:MAT", "HFP:CatchArea"), data)
  sem_obj <- as.psem(list(m_mat, m_map, m_hydro, m_wpc1, m_wpc2, m_div))
  sem_obj <- update(sem_obj, Water_PC1 %~~% Water_PC2_Organic)
  return(sem_obj)
}

all_sem_results <- list()
for (div in vars_div) {
  nat <- try(run_sem(process_data, div), silent=TRUE)
  data_large <- process_data[process_data$Type == "Large", ]
  lrg <- if (nrow(data_large) > 30) try(run_sem(data_large, div), silent=TRUE) else NULL
  data_small <- process_data[process_data$Type == "Small", ]
  sml <- if (nrow(data_small) > 30) try(run_sem(data_small, div), silent=TRUE) else NULL
  all_sem_results[[div]] <- list(National=nat, Large=lrg, Small=sml)
}

extract_coefs <- function(sem_obj, div, scale) {
  if (!inherits(sem_obj, "psem")) return(NULL)
  summ <- tryCatch(summary(sem_obj, standardize="scale", .progressBar=FALSE), error=function(e) NULL)
  if (is.null(summ)) return(NULL)
  coefs <- as.data.frame(summ$coefficients)
  coefs$Diversity_Metric <- div
  coefs$Scale            <- scale
  coefs
}
get_est <- function(d, resp, pred, coef_col) {
  v <- d[[coef_col]][d$Response == resp & d$Predictor == pred]
  if (length(v) == 0 || all(is.na(v))) 0 else v[1]
}

calc_indirect <- function(d, div, coef_col) {
  a1  <- get_est(d, "Qmon",              "HFP",  coef_col)
  a2  <- get_est(d, "Water_PC1",         "HFP",  coef_col)
  a3  <- get_est(d, "Water_PC2_Organic", "HFP",  coef_col)
  b12 <- get_est(d, "Water_PC1",         "Qmon", coef_col)
  b13 <- get_est(d, "Water_PC2_Organic", "Qmon", coef_col)
  c1  <- get_est(d, div, "Qmon",              coef_col)
  c2  <- get_est(d, div, "Water_PC1",         coef_col)
  c3  <- get_est(d, div, "Water_PC2_Organic", coef_col)
  data.frame(
    Direct_HFP                 = get_est(d, div, "HFP", coef_col),
    Qmon_mediated              = a1 * c1,
    Nutrients_ions_mediated    = (a2 + a1 * b12) * c2,
    Organic_pollution_mediated = (a3 + a1 * b13) * c3
  )
}
effects_results <- list()
for (div in vars_div) {
  for (sc in c("National", "Large", "Small")) {
    sem_obj <- all_sem_results[[div]][[sc]]
    if (inherits(sem_obj, "try-error") || is.null(sem_obj)) next
    coefs    <- extract_coefs(sem_obj, div, sc)
    if (is.null(coefs)) next
    coef_col <- if ("Std.Estimate" %in% names(coefs)) "Std.Estimate" else "Estimate"
    coefs[[coef_col]] <- as.numeric(coefs[[coef_col]])
    eff <- calc_indirect(coefs, div, coef_col)
    eff$Total_indirect <- eff$Qmon_mediated + eff$Nutrients_ions_mediated + eff$Organic_pollution_mediated
    eff$Total_effect   <- eff$Direct_HFP + eff$Total_indirect
    eff$Diversity      <- div
    eff$Scale          <- sc
    effects_results[[paste(div, sc)]] <- eff
  }
}
effects_df <- bind_rows(effects_results)
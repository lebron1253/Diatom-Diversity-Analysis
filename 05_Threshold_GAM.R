# 05_Threshold_GAM.R
# Threshold detection and 3D-GAM spatial grid prediction
# Packages: mgcv, segmented, chngpt, dplyr, tidyr

library(mgcv)
library(segmented)
library(chngpt)
library(dplyr)
library(tidyr)

data_raw <- read.table("metadata.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
responses <- c("TD", "FD", "PD", "LCBDtax", "LCBDfunc", "LCBDphylo")

df_clean <- data_raw %>%
  drop_na(all_of(c(responses, "HFP", "MAT", "CatchArea"))) %>%
  mutate(log_CatchArea = log10(CatchArea + 1))

# Utility functions
calc_aicc <- function(model, n) {
  k <- tryCatch(attr(logLik(model), "df"), error=function(e) NA)
  if (is.na(k)) k <- tryCatch(length(coef(model)) + 1, error=function(e) NA)
  if (is.na(k) || n <= k + 1) return(Inf)
  AIC(model) + (2 * k * (k + 1)) / (n - k - 1)
}

calc_r2 <- function(model, df) {
  preds <- as.numeric(predict(model))
  1 - sum((df$y - preds)^2, na.rm=TRUE) / sum((df$y - mean(df$y, na.rm=TRUE))^2, na.rm=TRUE)
}

# Threshold model fitting
fit_threshold_models <- function(x, y, delta_cutoff=-2) {
  df <- data.frame(x=x, y=y) %>% drop_na()
  n  <- nrow(df)
  results <- list()
  
  mod_lm <- lm(y ~ x, data=df)
  results$LM <- list(aicc=calc_aicc(mod_lm, n), bic=BIC(mod_lm),
                     r2=summary(mod_lm)$r.squared, threshold=NA, model=mod_lm)
  
  mod_gam <- tryCatch(gam(y ~ s(x, k=5, bs="cr"), data=df, method="REML"), error=function(e) NULL)
  if (!is.null(mod_gam)) {
    x_seq <- seq(min(df$x), max(df$x), length.out=500)
    pred  <- predict(mod_gam, newdata=data.frame(x=x_seq))
    d2    <- diff(diff(pred) / diff(x_seq)) / diff(x_seq[-1])
    results$GAM <- list(aicc=calc_aicc(mod_gam, n), bic=BIC(mod_gam),
                        r2=summary(mod_gam)$r.sq,
                        threshold=x_seq[which.max(abs(d2)) + 1], model=mod_gam)
  }
  
  mod_seg <- tryCatch(segmented(mod_lm, seg.Z=~x, psi=median(df$x)), error=function(e) NULL)
  if (!is.null(mod_seg))
    results$Segmented <- list(aicc=calc_aicc(mod_seg, n), bic=BIC(mod_seg),
                              r2=summary(mod_seg)$adj.r.squared,
                              threshold=mod_seg$psi[1,"Est."], model=mod_seg)
  
  for (type in c("step", "stegmented", "hinge", "upperhinge")) {
    mod <- tryCatch(chngptm(y~1, ~x, data=df, family="gaussian", type=type), error=function(e) NULL)
    if (!is.null(mod)) {
      key <- switch(type, step="Step", stegmented="StepSegmented",
                    hinge="Hinge", upperhinge="UpperHinge")
      results[[key]] <- list(aicc=calc_aicc(mod, n), bic=BIC(mod),
                             r2=calc_r2(mod, df), threshold=mod$chngpt, model=mod)
    }
  }
  
  aicc_vec    <- sapply(results, function(r) r$aicc)
  lm_aicc     <- aicc_vec["LM"]
  best_name   <- names(which.min(aicc_vec))
  best        <- results[[best_name]]
  delta_vs_lm <- best$aicc - lm_aicc
  
  threshold_supported <- !is.na(best$threshold) &&
    best_name != "LM" && delta_vs_lm <= delta_cutoff
  
  list(best_model       = best_name,
       threshold        = if (threshold_supported) best$threshold else NA,
       threshold_supported = threshold_supported,
       delta_aicc_vs_lm = delta_vs_lm,
       model            = best$model)
}

# National threshold detection
global_thresh <- list()

for (y_var in responses) {
  res <- fit_threshold_models(df_clean$HFP, df_clean[[y_var]])
  global_thresh[[y_var]] <- res$threshold
}

# 3D-GAM spatial grid prediction
grid_res     <- 100
pred_seq_res <- 100

mat_seq      <- seq(min(df_clean$MAT,           na.rm=TRUE), max(df_clean$MAT,           na.rm=TRUE), length.out=grid_res)
catcharea_seq <- seq(min(df_clean$log_CatchArea, na.rm=TRUE), max(df_clean$log_CatchArea, na.rm=TRUE), length.out=grid_res)
hfp_seq      <- seq(min(df_clean$HFP,           na.rm=TRUE), max(df_clean$HFP,           na.rm=TRUE), length.out=pred_seq_res)

grid_base       <- expand.grid(MAT=mat_seq, log_CatchArea=catcharea_seq)
all_grids_results <- data.frame()

for (var in responses) {
  form    <- as.formula(paste(var, "~ te(HFP, MAT, log_CatchArea, k=c(5,5,5), bs=c('cr','cr','cr'))"))
  mod_3d  <- gam(form, data=df_clean, method="REML")
  thresholds <- numeric(nrow(grid_base))
  
  for (i in seq_len(nrow(grid_base))) {
    newdata <- data.frame(HFP=hfp_seq, MAT=grid_base$MAT[i], log_CatchArea=grid_base$log_CatchArea[i])
    pred_y  <- predict(mod_3d, newdata=newdata)
    res_i   <- tryCatch(fit_threshold_models(hfp_seq, pred_y), error=function(e) list(threshold=NA))
    thresholds[i] <- res_i$threshold
  }
  
  grid_var          <- grid_base
  grid_var$HFPsite  <- thresholds
  grid_var$Index    <- var
  all_grids_results <- rbind(all_grids_results, grid_var)
}

# Context-dependent threshold summaries
national_thresholds <- c(
  "TD"=17.53, "FD"=16.67, "PD"=7.62,
  "LCBD_Tax"=9.92, "LCBD_Func"=11.97, "LCBD_Phylo"=7.62)

df_plot <- all_grids_results %>% filter(!is.na(HFPsite))

mat_breaks  <- quantile(df_plot$MAT,           probs=c(0, 0.33, 0.66, 1), na.rm=TRUE)
area_breaks <- quantile(df_plot$log_CatchArea,  probs=c(0, 0.33, 0.66, 1), na.rm=TRUE)

df_risk <- df_plot %>%
  filter(Index %in% names(national_thresholds)) %>%
  mutate(
    National_Val   = national_thresholds[Index],
    MAT_Gradient   = case_when(
      MAT <= mat_breaks[2]                           ~ sprintf("Low MAT (<%.1f°C)",      mat_breaks[2]),
      MAT > mat_breaks[2] & MAT <= mat_breaks[3]    ~ sprintf("Mid MAT (%.1f-%.1f°C)",  mat_breaks[2], mat_breaks[3]),
      MAT > mat_breaks[3]                            ~ sprintf("High MAT (>%.1f°C)",     mat_breaks[3])),
    Area_Gradient  = case_when(
      log_CatchArea <= area_breaks[2]                          ~ sprintf("Small CA (<%.1f)",       area_breaks[2]),
      log_CatchArea > area_breaks[2] & log_CatchArea <= area_breaks[3] ~ sprintf("Mid CA (%.1f-%.1f)", area_breaks[2], area_breaks[3]),
      log_CatchArea > area_breaks[3]                           ~ sprintf("Large CA (>%.1f)",       area_breaks[3])))
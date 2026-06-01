# 02_BRT.R
# Boosted regression trees for six diatom diversity facets across four spatial scales
# Packages: dismo, gbm, dplyr

library(dismo)
library(gbm)
library(dplyr)

RESPONSE_VARS <- c("TD", "FD", "PD", "LCBDtax", "LCBDfunc", "LCBDphylo")

PREDICTOR_VARS <- c(
  "Elevation", "Gradient", "SPI", "Qmon", "Qann", "CatchArea", "Order", "Out_dist",
  "MAT", "Tsea", "Ttrend", "MAP", "Psea", "Ptrend",
  "HFP", "Urban", "Cropland",
  "TN", "TP", "EC", "pH", "DO", "COD")

VAR_CATEGORIES <- list(
  Hydrography    = c("Elevation", "Gradient", "SPI", "Qmon", "Qann", "CatchArea", "Order", "Out_dist"),
  Climate        = c("MAT", "Tsea", "Ttrend", "MAP", "Psea", "Ptrend"),
  Human_pressure = c("HFP", "Urban", "Cropland"),
  Water_quality  = c("TN", "TP", "EC", "pH", "DO", "COD"))

run_brt <- function(data, response, predictors) {
  df <- na.omit(data[, c(response, predictors)])
  set.seed(123)
  model <- tryCatch({
    gbm.step(
      data            = df,
      gbm.x           = which(names(df) %in% predictors),
      gbm.y           = which(names(df) == response),
      family          = "gaussian",
      tree.complexity = 3,
      learning.rate   = 0.001,
      bag.fraction    = 0.5,
      max.trees       = 10000, 
      plot.main       = FALSE,
      silent          = TRUE)
  }, error = function(e) { NULL })
  return(model)
}

get_importance <- function(model) {
  var_imp <- summary(model, plotit = FALSE) %>%
    rename(Variable = var, Importance = rel.inf) %>%
    arrange(desc(Importance))
  
  cat_imp <- data.frame(
    Category   = names(VAR_CATEGORIES),
    Importance = sapply(VAR_CATEGORIES, function(vars)
      sum(var_imp$Importance[var_imp$Variable %in% vars]))) %>%
    arrange(desc(Importance))
  list(var_imp = var_imp, cat_imp = cat_imp)
}

run_scale <- function(data, label) {
  var_summary <- data.frame(Variable = PREDICTOR_VARS)
  for (resp in RESPONSE_VARS) {
    model <- run_brt(data, resp, PREDICTOR_VARS)
    if (is.null(model)) {
      var_summary <- mutate(var_summary, !!resp := NA)
      next}
    imp <- get_importance(model)
    var_summary <- left_join(var_summary,
                             imp$var_imp %>% rename(!!resp := Importance),
                             by = "Variable")
    cat_summary <- left_join(cat_summary,
                             imp$cat_imp %>% rename(!!resp := Importance),
                             by = "Category")}
}

data <- read.table("metadata.txt", header = TRUE, sep = "\t")

# National
run_scale(data, "national")
# Large/small rivers
run_scale(data[data$Order >= 5, ], "large_river")
run_scale(data[data$Order <= 4, ], "small_river")
# Basin level
for (b in unique(data$basin)) {
  run_scale(data[data$basin == b, ], paste0("basin_", b))}
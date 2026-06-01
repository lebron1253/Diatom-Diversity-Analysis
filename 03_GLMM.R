# 03_GLMM.R
# GLMMs for diversity facets, CWM traits, and basin-level HFP slope analysis
# Packages: glmmTMB, MuMIn, FD, e1071, car, dplyr, broom

library(glmmTMB)
library(MuMIn)
library(FD)
library(e1071)
library(car)
library(dplyr)
library(broom)

source("01_Diversity_Calculation.R")

df_raw     <- read.table("metadata.txt",       header=TRUE, sep="\t", stringsAsFactors=FALSE)
comm       <- read.table("abundace.txt",              header=TRUE, sep="\t", row.names=1, check.names=FALSE)
traits     <- read.table("traits.txt",          header=TRUE, sep="\t", row.names=1, check.names=FALSE)
alpha_div  <- alpha_result
lcbd_basin <- res_basin

df_raw$Basin <- as.factor(df_raw$Basin)

# CWM calculation
common_sp     <- intersect(colnames(comm), rownames(traits))
cwm           <- functcomp(x=traits[common_sp,], a=as.matrix(comm[,common_sp]), CWM.type="all")
colnames(cwm) <- paste0("CWM_", colnames(cwm))
cwm$Site      <- rownames(cwm)
df_raw        <- left_join(df_raw, cwm, by="Site")

# Variable lists
response_vars <- c("TD", "FD", "PD", "LCBD_Tax", "LCBD_Func", "LCBD_Phylo")
env_vars      <- c("Elevation", "CatchArea", "MAT", "HFP", "EC", "TP", "COD")
trait_vars    <- grep("^CWM_", names(df_raw), value=TRUE)
all_vars      <- c(response_vars, env_vars, trait_vars)

# Log transformation
log_transform <- function(data, vars) {
  for (v in vars) {
    if (!v %in% names(data)) next
    sk <- skewness(data[[v]], na.rm=TRUE)
    if (!is.na(sk) && abs(sk) > 1) {
      x  <- data[[v]]
      mn <- min(x, na.rm=TRUE)
      data[[v]] <- if (mn > 0) log10(x) else log10(x + abs(mn) + 1)
    }
  }
  data
}

df_proc  <- log_transform(df_raw, all_vars)
df_final <- df_proc
df_final[all_vars[all_vars %in% names(df_final)]] <-
  scale(df_final[all_vars[all_vars %in% names(df_final)]])

# ── VIF check ──────────────────────────────────────────────────────────────
print(vif(lm(TD ~ Elevation + CatchArea + MAT + HFP + EC + TP + COD, data=df_final)))

# ── Interaction classification ─────────────────────────────────────────────
classify_interaction <- function(b1, b2, bi, p1, p2, pi, n1, n2, alpha=0.05) {
  if (any(is.na(c(b1, b2, bi, p1, p2, pi)))) return(NA)
  
  s1 <- p1 < alpha
  s2 <- p2 < alpha
  si <- pi < alpha

  if (!si) {
    if (s1 & s2)  return("Additive")
    if (s1 & !s2) return(paste0("Dominance (", n1, ")"))
    if (!s1 & s2) return(paste0("Dominance (", n2, ")"))
    return("None")
  }
 
  if (!s1 | !s2) {
    return("Unclassified") 
  }
  
  # Calculate Reversal
  exp_eff <- b1 + b2
  obs_eff <- b1 + b2 + bi
  
  if (exp_eff != 0 && sign(obs_eff) != sign(exp_eff)) {
    return("Reversal")
  }
  
  # Remaining significant interactions based on coefficient signs
  if (sign(b1) == sign(b2)) {
    if (sign(bi) == sign(b1)) {
      return("Synergistic")
    } else {
      return("Antagonistic")
    }
  } else {
    # Opposing: main effects have opposite signs
    if (sign(bi) == sign(b1)) {
      return(paste0("Opposing (", n1, ")"))
    } else if (sign(bi) == sign(b2)) {
      return(paste0("Opposing (", n2, ")"))
    } else {
      return("Unclassified")
    }
  }
}

# Diversity GLMM
run_glmm(df_final,                          response_vars, "diversity_national")
run_glmm(df_final[df_final$Type=="Large",], response_vars, "diversity_large")
run_glmm(df_final[df_final$Type=="Small",], response_vars, "diversity_small")

# CWM trait GLMM
run_glmm(df_final,                          trait_vars, "CWM_national")
run_glmm(df_final[df_final$Type=="Large",], trait_vars, "CWM_large")
run_glmm(df_final[df_final$Type=="Small",], trait_vars, "CWM_small")

# Basin-level HFP slope analysis
df_basin <- df_raw %>%
  left_join(alpha_div, by="Site") %>%
  left_join(lcbd_basin[, c("Site","LCBDtax","LCBDfunc","LCBDphylo")], by="Site") %>%
  rename(LCBD_Tax=LCBDtax, LCBD_Func=LCBDfunc, LCBD_Phylo=LCBDphylo)

df_basin <- log_transform(df_basin, c(response_vars, env_vars))
df_basin[c(response_vars, env_vars)] <- scale(df_basin[c(response_vars, env_vars)])

slope_df <- lapply(response_vars, function(resp) {
  df_basin %>%
    group_by(Basin) %>%
    do(tidy(lm(as.formula(paste(resp,
                                "~ HFP + Elevation + EC + TP + COD + CatchArea + MAT")), data=.))) %>%
    filter(term == "HFP") %>%
    mutate(Response = resp)
}) %>% bind_rows()

basin_means <- df_raw %>%
  group_by(Basin) %>%
  summarise(Mean_MAT       = mean(MAT,       na.rm=TRUE),
            Mean_CatchArea = mean(CatchArea,  na.rm=TRUE))

slope_df  <- left_join(slope_df, basin_means, by="Basin")

slope_reg <- lapply(response_vars, function(resp) {
  sub <- slope_df[slope_df$Response == resp, ]
  r_mat  <- summary(lm(estimate ~ Mean_MAT,      data=sub))
  r_area <- summary(lm(estimate ~ Mean_CatchArea, data=sub))
  data.frame(Response  = resp,
             MAT_slope  = r_mat$coefficients[2,1],  MAT_p  = r_mat$coefficients[2,4],
             Area_slope = r_area$coefficients[2,1],  Area_p = r_area$coefficients[2,4])
}) %>% bind_rows()
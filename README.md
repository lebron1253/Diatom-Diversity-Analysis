# River hydrography modulates anthropogenic impacts on benthic diatom diversity in China

This repository contains the R scripts and analytical pipelines required to reproduce the statistical analyses and figures presented in the manuscript.

# All data files are included in this repository. 
Please download the datasets and place them directly into the root working directory before executing the scripts.

# Operating System: Windows 10/11, macOS, or Linux.
Software: R (version 4.5.0 or higher recommended) and RStudio.
Dependencies: The scripts rely on the following R packages. 
You can install them using standard `install.packages()` or `BiocManager::install()` commands:
  `hillR`, `vegan`, `ape`, `cluster`, `adespatial`, `SYNCSA`, `dismo`, `gbm`, `dplyr`, `tidyr`, `glmmTMB`, `MuMIn`, `FD`, `e1071`, `car`, `broom`, `piecewiseSEM`, `mgcv`, `segmented`, `chngpt`, `ggplot2`, `patchwork`, `ggpubr`.

# Scripts should be executed in the following sequential order. Ensure your working directory is set to the folder containing the scripts and the downloaded Figshare data.

`01_Diversity_Calculation.R` — Calculates taxonomic, functional, and phylogenetic diversity (Hill q = 1) and local contributions to beta diversity (LCBD).
`02_BRT.R` — Constructs boosted regression trees (BRTs) to evaluate the relative importance of environmental drivers.
`03_GLMM.R` — Fits Generalized Linear Mixed Models (GLMMs) to examine HFP–biodiversity relationships and conducts basin-level slope analysis.
`04_SEM.R` — Constructs piecewise structural equation models (piecewise SEMs) to evaluate indirect ecological pathways mediated by hydrological and hydrochemical factors.
`05_Threshold_GAM.R` — Performs threshold detection (using segmented and step models) and 3D-GAM spatial grid prediction. *(Note: The spatial grid prediction involves fitting thousands of models. Depending on your CPU, this script may take 1-3 hours to complete).

For any questions regarding the code or data, please contact Yu Ma at yma@mail.bnu.edu.cn

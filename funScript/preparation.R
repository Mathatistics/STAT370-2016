## Loading Packages
req_pkgs <- c('pls', 'MASS', 'plsVarSel', 'ggplot2', 'reshape2', 'plyr', 'dplyr', 'data.table', 'purrr')
invisible(lapply(req_pkgs, require, character.only = TRUE, warn.conflicts = FALSE, quietly = TRUE))
rm(req_pkgs)

## Load Datasets
load("Data/MALDITOF_Milk.RData")
load("Data/GeneExpr_Prostate.RData")
load("Data/NIR_Raman_PUFA.RData")
GeneExpr_Prostate <- GeneExpr_Prostate[GeneExpr_Prostate$train,]

## Soursing Relevant Scripts
source('funScript/basicFun.R')
source('funScript/expressions.R')

## Setting up Training and Test Set
MALDITOF_Milk$train <- ifelse(as.numeric(rownames(MALDITOF_Milk)) <= 120, TRUE, FALSE)
GeneExpr_Prostate <- splitMe(GeneExpr_Prostate, 0.7)
NIR_Raman_PUFA <- splitMe(NIR_Raman_PUFA, 0.7)

## Initial Preparation
eval(prepInit)

## Model Fitting ::: takes long time if not loaded from previous
load_or_run('Data/FittedModels.Rdata', list(model.fit))

if (!file.exists("Data/FittedModels.Rdata")) {
  
} else {
  load("Data/FittedModels.Rdata")
}
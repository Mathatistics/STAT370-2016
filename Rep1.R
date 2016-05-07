## ----frontMatter, child='FrontMatters.Rnw'-------------------------------
rm(list = ls())

## ----loadPackage, echo = FALSE, message=FALSE, warning=FALSE-------------
req.pkgs <- c('pls', 'MASS', 'plsVarSel', 'data.table', 'ggplot2', 'reshape2', 'plyr', 'dplyr', 'xtable', 'cowplot', 'nlme')
invisible(lapply(req.pkgs, require, character.only = TRUE))
theme_set(theme_bw(10))

## ----knitrSetup, include=FALSE-------------------------------------------
knitr::opts_chunk$set(echo = FALSE, comment = NA, size = 'small', fig.height = 4, fig.pos = 'H')

## ----dataLoad------------------------------------------------------------
load('Data/MALDITOF_Milk.RData')
load('Data/GeneExpr_Prostate.RData')
load('Data/NIR_Raman_PUFA.RData')
GeneExpr_Prostate <- GeneExpr_Prostate[GeneExpr_Prostate$train,]
GeneExpr_Prostate <- within(GeneExpr_Prostate, {
  colnames(GeneExpr) <- 1:ncol(GeneExpr)
  rownames(GeneExpr) <- 1:nrow(GeneExpr)
})
data <- list(milk = MALDITOF_Milk, nir = NIR_Raman_PUFA, gene = GeneExpr_Prostate)

source('funScript/basicFun.R')
source("funScript/plda.R")

## ----modelling, child='Modelling.Rnw'------------------------------------

## ----DataSplit-----------------------------------------------------------
data$gene <- splitMe(data$gene, 0.7, seed = 3000)
data$raman <- data$nir <- splitMe(data$nir, 0.7, seed = 3000)
data$milk$train <- ifelse(1:nrow(data$milk) %between% c(61, 180), TRUE, FALSE)
xy.info <- list(
  c(x = 'MS', y = "milk"), c(x = 'NIR', y = 'PUFA'), c(x = 'GeneExpr', y = 'tumor'), c(x = 'Raman', y = 'PUFA')
)
names(xy.info) <- names(data)

## ----Model_Fitting-------------------------------------------------------

pls <- list()
## Initilizing Parallel
pls.options(parallel = 6)

## PLS model for Milk dataset
pls$milk <- lapply(colnames(data$milk$milk), function(rsp){
    plsr(milk[, rsp] ~ MS, data = data$milk[data$milk$train, ], 
           validation = "CV", segments = 30, segment.type = "consecutive", 
           jackknife = TRUE)
})
names(pls$milk) <- colnames(data$milk$milk)

## PLS model for NIR dataset
pls$nir <- lapply(colnames(data$nir$PUFA), function(rsp){
    plsr(PUFA[, rsp] ~ NIR, data = data$nir[data$nir$train, ], 
         validation = "CV", jackknife = TRUE)

})

## PLS model for Raman dataset
pls$raman <- lapply(colnames(data$raman$PUFA), function(rsp){
    plsr(PUFA[, rsp] ~ Raman, data = data$raman[data$raman$train, ], 
         validation = "CV", jackknife = TRUE)
})
names(pls$nir) <- names(pls$raman) <- colnames(data$nir$PUFA)

gene <- list()

gene$plda <- plda(data$gene$GeneExpr[data$gene$train,],
                  cat.resp = data$gene$tumor[data$gene$train],
                  fn = 'plsr', fitComp = 25,
                  split = 10)

gene$r2 <- classR2(gene$plda, 
        newX = data$gene$GeneExpr[!data$gene$train,],
        newY = data$gene$tumor[!data$gene$train], ncomp = 25)

# Variables Selection --------------------------------------------------------
opt.comp <- unique(gene$r2$max.correct.prop$comp)

## VIP Method (vip > 1)
gene$subset$vip <- VIP(gene$plda, opt.comp = opt.comp)
gene$subvaridx$vip <- as.numeric(which(gene$subset$vip > 1))
gene$submodel$vip <- plda(data$gene$GeneExpr[data$gene$train, gene$subvaridx$vip],
                  cat.resp = data$gene$tumor[data$gene$train],
                  fn = 'plsr', fitComp = opt.comp,
                  split = 10)
gene$sub_r2$vip <- classR2(gene$submodel$vip, 
        newX = data$gene$GeneExpr[!data$gene$train, gene$subvaridx$vip],
        newY = data$gene$tumor[!data$gene$train], ncomp = opt.comp)

## MCUV-PLS Method
gene$subset$mcuv <- mcuve_pls(data$gene$tumor[data$gene$train],
                         data$gene$GeneExpr[data$gene$train,],
                         ncomp = opt.comp)
gene$subvaridx$mcuv <- gene$subset$mcuv[[1]]
gene$submodel$mcuv <- plda(data$gene$GeneExpr[data$gene$train, gene$subvaridx$mcuv],
                  cat.resp = data$gene$tumor[data$gene$train],
                  fn = 'plsr', fitComp = opt.comp,
                  split = 10)
gene$sub_r2$mcuv <- classR2(gene$submodel$mcuv, 
        newX = data$gene$GeneExpr[!data$gene$train, gene$subvaridx$mcuv],
        newY = data$gene$tumor[!data$gene$train], ncomp = opt.comp)

## Shaving Method
gene$subset$shaving <- shaving(data$gene$tumor[data$gene$train],
                         data$gene$GeneExpr[data$gene$train,],
                         ncomp = opt.comp, method = "sMC")
gene$subvaridx$shaving <- gene$subset$shaving$variables[
  which.min(gene$subset$shaving$error)][[1]]
gene$submodel$shaving <- plda(data$gene$GeneExpr[data$gene$train, gene$subvaridx$shaving],
                  cat.resp = data$gene$tumor[data$gene$train],
                  fn = 'plsr', fitComp = opt.comp,
                  split = 10)
gene$sub_r2$shaving <- classR2(gene$submodel$shaving, 
        newX = data$gene$GeneExpr[!data$gene$train, gene$subvaridx$shaving],
        newY = data$gene$tumor[!data$gene$train], ncomp = opt.comp)


## Truncated Method
gene$submodel$truncation <- plda(data$gene$GeneExpr[data$gene$train,],
                  cat.resp = data$gene$tumor[data$gene$train],
                  fn = 'truncation', fitComp = 25,
                  split = 10, truncation = "Lenth", trunc.width = 0.95)

gene$sub_r2$truncation <- classR2(gene$submodel$truncation, 
                                  data$gene$GeneExpr[!data$gene$train, ], 
                                  newY = data$gene$tumor[!data$gene$train], 
                                  ncomp = opt.comp)


## ----Validation----------------------------------------------------------
r2 <- lapply(names(pls), function(mdls){
    getValidation(pls[[mdls]], newdata = data[[mdls]][!data[[mdls]]$train, ])
})
names(r2) <- names(pls)  

## ----ValidataionPlots, fig.height=4, results='hide', fig.cap="Maximum variation explained on response in Milk, NIR and Raman dataset and maximium correct classification proportion for gene Expression dataset. Each line represent a model for which the R2 predicted is calculated with number of components chosen by cross-validation. The number at the end of each line represent the components used."----
plt.list <- lapply(names(r2), function(x){
  plot(r2[[x]], full = FALSE) + ggtitle(toupper(x)) + theme(title = element_text(size = 7))
})
plt.list$gene <- plot(gene$r2, full = F) + ggtitle("Gene Expression") + theme(title = element_text(size = 7))
plot_grid(plotlist = plt.list, align = 'v', nrow = 2)


## ----classification, child='Classification.Rnw'--------------------------



## ----variableSelection, child='VariableSelection.Rnw'--------------------

## ----Subset_Models-------------------------------------------------------

  ## Initilize List
subset <- subvaridx <- subdata <- submodel <- sub_r2 <- list()

## Filter with VIP
subset$vip <- lapply(names(pls), function(dta, vald = r2){
  mdls <- lapply(names(pls[[dta]]), function(mdl){
    opt.comp <- vald[[dta]][[2]][response == mdl & estimate == "CV", comp]
    vip <- VIP(pls[[dta]][[mdl]], opt.comp = opt.comp)
    vip.idx <- as.numeric(which(vip > 1))
    subdata <- within(data[[dta]], {assign(xy.info[[dta]][['x']], data[[dta]][, xy.info[[dta]]['x']][, vip.idx])})
    subMdl <- update(pls[[dta]][[mdl]], . ~ ., data = subdata[subdata$train, ])
    
    ## Output (Export to respective List)
    subvaridx$vip[[dta]][[mdl]] <<- vip.idx
    subdata$vip[[dta]][[mdl]] <<- subdata
    submodel$vip[[dta]][[mdl]] <<- subMdl
    
    return(vip)
  })
  names(mdls) <- names(pls[[dta]])
  sub_r2$vip[[dta]] <<- getValidation(submodel$vip[[dta]], subdata$vip[[dta]])
  return(mdls)
})
names(subset$vip) <- names(pls)

## Fitting mcuv_pls
subset$mcuv <- lapply(names(pls), function(dta, vald = r2){
  mdls <- lapply(names(pls[[dta]]), function(mdl){
    opt.comp <- vald[[dta]][[2]][response == mdl & estimate == "CV", comp]
    mcuv.test <- mcuve_pls(data[[dta]][[xy.info[[dta]]["y"]]][, mdl], data[[dta]][[xy.info[[dta]]["x"]]])[[1]]
    mcuv.idx <- mcuv.test
    subdata <- within(data[[dta]], {assign(xy.info[[dta]][['x']], data[[dta]][, xy.info[[dta]]['x']][, mcuv.idx])})
    subMdl <- update(pls[[dta]][[mdl]], . ~ ., data = subdata[subdata$train, ])
    
    ## Output (Export to respective List)
    subvaridx$mcuv[[dta]][[mdl]] <<- mcuv.idx
    subdata$mcuv[[dta]][[mdl]] <<- subdata
    submodel$mcuv[[dta]][[mdl]] <<- subMdl
    
    return(mcuv.test)
  })
  names(mdls) <- names(pls[[dta]])
  sub_r2$mcuv[[dta]] <<- getValidation(submodel$mcuv[[dta]], subdata$mcuv[[dta]])
  return(mdls)
})
names(subset$mcuv) <- names(pls)

## Fitting shaving_pls
subset$shaving <- lapply(names(pls), function(dta, vald = r2){
  mdls <- lapply(names(pls[[dta]]), function(mdl){
    opt.comp <- vald[[dta]][[2]][response == mdl & estimate == "CV", comp]
    shaving.test <- shaving(data[[dta]][[xy.info[[dta]]["y"]]][, mdl], data[[dta]][[xy.info[[dta]]["x"]]])
    shaving.idx <- shaving.test$variables[which.min(shaving.test$error)][[1]]
    
    subdata <- within(data[[dta]], {assign(xy.info[[dta]][['x']], data[[dta]][, xy.info[[dta]]['x']][, shaving.idx])})
    subMdl <- update(pls[[dta]][[mdl]], . ~ ., data = subdata[subdata$train, ])
    
    ## Output (Export to respective List)
    subvaridx$shaving[[dta]][[mdl]] <<- shaving.idx
    subdata$shaving[[dta]][[mdl]] <<- subdata
    submodel$shaving[[dta]][[mdl]] <<- subMdl
    
    return(shaving.test)
  })
  names(mdls) <- names(pls[[dta]])
  sub_r2$shaving[[dta]] <<- getValidation(submodel$shaving[[dta]], subdata$shaving[[dta]])
  return(mdls)
})
names(subset$shaving) <- names(pls)

## Fitting truncation_pls
invisible(
  lapply(names(pls), function(dta, vald = r2){
    lapply(names(pls[[dta]]), function(rsp){
      opt.comp <- vald[[dta]][[2]][response == rsp & estimate == "CV", comp]
      mdl.call <- as.list(getCall(pls[[dta]][[rsp]]))
      mdl.call[[1]] <- as.name("truncation")
      mdl.call$truncation = "Lenth"
      mdl.call$trunc.width = 0.95
      mdl.call$ncomp = opt.comp
      submodel$truncation[[dta]][[rsp]] <<- eval(as.call(mdl.call))
    })
    sub_r2$truncation[[dta]] <<- getValidation(submodel$truncation[[dta]], data[[dta]][!data[[dta]]$train, ])
  })
)


## ----truncPlot, fig.height=3, fig.cap='Histogram of Loading Weights for full and truncated model for (left) MALDI-TOF milk data with response cow (right) Gene Expression data. Corresponding inset shows the maximum R-sq predicted for train, CV and test.', message=FALSE, warning=FALSE----
trunc.cow <- truncPlot(pls$milk$cow, submodel$truncation$milk$cow, bins = 1000, newdata = data$milk[!data$milk$train, ], inset = T, invert = T)
trunc.gene <- truncPlot(gene$plda, gene$submodel$truncation, bins = 1000, inset = T, invert = F)
grid_arrange_shared_legend(trunc.cow, trunc.gene, n.col = 2)


## ----AnalysisDiscussion, child='AnalysisDiscussion.Rnw'------------------

## ----TransferGene--------------------------------------------------------
r2$gene <- gene$r2
for (x in names(sub_r2)) {
  sub_r2[[x]]$gene <- gene$sub_r2[[x]]
}

## ----MetaDataset---------------------------------------------------------
full.test.r2 <- setnames(rbindlist(lapply(r2, function(rsp){
      rsp[[2]][estimate == "test"]
  }), idcol = TRUE), '.id', "data")
full.test.r2 <- cbind(method = "complete", full.test.r2)

sub.test.r2 <- (setnames(rbindlist(lapply(sub_r2, function(mthd) {
    setnames(rbindlist(lapply(mthd, function(rsp){
        rsp[[2]][estimate == "test"]
    }), idcol = TRUE), '.id', "data")
}), idcol = TRUE), '.id', 'method'))

test.r2 <- setnames(rbindlist(list(full = full.test.r2, subset = sub.test.r2), idcol = TRUE), '.id', 'type')
test.r2[, c("method", "data") := list(as.factor(method), as.factor(data))]
rm(full.test.r2, sub.test.r2)

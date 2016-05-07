## ----frontMatter, child='FrontMatters.Rnw'-------------------------------

## ----loadPackage, echo = FALSE, message=FALSE, warning=FALSE-------------
req.pkgs <- c('pls', 'MASS', 'plsVarSel', 'data.table', 'ggplot2', 'reshape2', 'plyr', 'dplyr', 'xtable', 'cowplot', 'nlme', 'mixlm')
invisible(lapply(req.pkgs, require, character.only = TRUE))
theme_set(theme_bw(10))

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


## ----dataset_desc, results='asis'----------------------------------------
dataDesc <- data.table(
  Datasets = paste0('\\texttt{', c('MALDITOF\\_Milk', 'NIR\\_Raman\\_PUFA', 'GeneExpr\\_Prostate'),'}'),
  Types = c('Spectrometric', 'Spectroscopic', 'Microarray'),
  `Response Variables` = paste0('\\texttt{', c(paste('cow', 'goat', 'ewe', sep = ', '), 
                                              paste('PUFA\\_total', 'PUFA\\_fat', sep = ', '), 
                                              'tumor'), '}'),
  `Predictor Variables` = c(
    "Mass Spectra (\\texttt{MS})", "\\texttt{NIR} and \\texttt{Raman}", "Gene Expression (\\texttt{GeneExpr})"
  )
)

print.xtable(xtable(dataDesc, align = 'llllX', caption = "Overview of Datasets", 
                    label = 'tbl:data-n-nature'), include.rownames = FALSE, 
             tabular.environment = "tabularx", width = '0.9\\textwidth',
             sanitize.text.function = function(x) {x}, size = 'small', table.placement = 'H')

## ----samplePlot, fig.height=3.8, fig.cap="Sample observations from a) MALDI-TOF (top-left) b) NIR (top-right) c) Raman (bottom-left) and d) GeneExpression (bottom-right)", fig.pos='H'----
plts <- list(
  getSeries(MALDITOF_Milk$MS, nsamp = 5, title = "MALDI-TOF milk dataset"),
  getSeries(NIR_Raman_PUFA$NIR, nsamp = 5, title = "NIR dataset"),
  getSeries(NIR_Raman_PUFA$Raman, nsamp = 5, title = "Raman dataset"),
  getSeries(GeneExpr_Prostate$GeneExpr, nsamp = 5, title = "Gene Expression dataset")
)
do.call(gridExtra::grid.arrange, plts)

## ----varsel_desc, results='asis', eval = F-------------------------------
varSelDesc <- data.table(
  Methods = c('Jack Knife', '\\texttt{mcuv}', 'Shaving (\\texttt{sMC})','Truncation'),
  Types = c('Filter', 'Wrapper', 'Wrapper', 'Embeded'),
  `Description` = c(
    'Jack Knifing is used for variable selection where variables having p-value less than 0.05 are only selected.',
    'Predictor variables having lower ``importance" than the artificial noise added to the predictors set are eliminated before model fitting. The process is repeated until a criterion is attained.',
    'Variables are sorted with respect to some importance measure using \\texttt{sMC} from which least informative variables are eliminated using a threshold value. A model is fitted again with remaining variables and the procedure is repeated until maximum model performance is achieved.',
    'During the model fitting, based on confidence interval for modeling, loading weights are truncated around their median before returning to the NIPALS algorithm again.'
  )
)

varSelDesc.xtbl <- xtable(varSelDesc, align = 'lllX', caption = "Variable Selection Methods (\\cite{mehmood2012review})", label = 'tbl:varsel-method')
print.xtable(varSelDesc.xtbl, include.rownames = FALSE, 
             tabular.environment = "tabularx", width = '0.9\\textwidth',
             sanitize.text.function = function(x) {x}, size = 'small', 
             table.placement = 'H', hline.after = c(-1, 0, 0:nrow(varSelDesc)))

## ----modelling, child='Modelling.Rnw'------------------------------------

## ----PCA, fig.height=2---------------------------------------------------
pca <- list()
pca$milk <- prcomp(data$milk$MS)
pca$raman <- prcomp(data$nir$Raman)
pca$NIR <- prcomp(data$nir$NIR)
pca$gene <- prcomp(data$gene$GeneExpr)

## ----screePlot, eval = FALSE---------------------------------------------
pcaPlots <- lapply(names(pca), function(x){
  plotPCA(pca[[x]], npcs = 20, title = x, base.size = 8)
})
do.call(gridExtra::grid.arrange, `$<-`(pcaPlots, ncol, 4))

## ----scorePlots, fig.height=4.5, fig.cap='Principal components of (top-left) MALDI-TOF dataset with colored proportion of cow, goat and ewe milk, (top-right) Gene Expression dataset with presence and absense of tumer on the sample, (bottom) NIR and Raman dataset with opacity factor with percentage of (left) PUFA in total sample and (right) PUFA in total fat content', fig.show='H'----
pca.scoreplot <- pca.score <- list()
pca.score$Milk <- melt(data.table(pca$milk$x[, 1:2], unclass(MALDITOF_Milk$milk)), 1:2)
pca.score$Gene <- melt(data.table(pca$gene$x[, 1:2], tumor = as.factor(GeneExpr_Prostate$tumor)), 1:2)
pca.score$NIR <- melt(data.table(pca$NIR$x[, 1:2], unclass(NIR_Raman_PUFA$PUFA)), 1:2)[, .(PC1, PC2, value = value/max(value)), by = variable]
pca.score$Raman <- melt(data.table(pca$raman$x[, 1:2], unclass(NIR_Raman_PUFA$PUFA)), 1:2)[, .(PC1, PC2, value = value/max(value)), by = variable]

pca.scoreplot <- lapply(names(pca.score), function(i){
  if (i == "Milk") {
    mapping <- aes(PC1, PC2, fill = variable, alpha = value)
    my_point <- geom_point(shape = 21)
  } else if (i == "Gene") {
    mapping <- aes(PC1, PC2, fill = as.factor(value))
    my_point <- geom_point(shape = 21)
  } else {
    mapping <- aes(PC1, PC2, alpha = value)
    my_point <- geom_point(shape = 21, fill = "gray")
  }
  plt <- ggplot(pca.score[[i]], mapping) + 
    my_point +
    scale_alpha_continuous(breaks = NULL) +
    theme(legend.position = "top", legend.title = element_blank()) +
    ggtitle(paste('Scoreplot for', i, 'dataset')) +
    theme(title = element_text(size = rel(0.9)))
  if (i == "NIR" | i == "Raman") {
    plt <- plt + facet_grid(.~variable)
  }
  return(plt)
})

do.call(gridExtra::grid.arrange, `$<-`(`$<-`(pca.scoreplot, 'ncol', 2), 'heights', c(3,2)))

## ----DataSplit-----------------------------------------------------------
data$gene <- splitMe(data$gene, 0.7)
data$raman <- data$nir <- splitMe(data$nir, 0.7)
data$milk$train <- ifelse(1:nrow(data$milk) <= 120, TRUE, FALSE)
xy.info <- list(
  c(x = 'MS', y = "milk"), c(x = 'NIR', y = 'PUFA'), c(x = 'GeneExpr', y = 'tumor'), c(x = 'Raman', y = 'PUFA')
)
names(xy.info) <- names(data)

## ----Model_Fitting-------------------------------------------------------
if (file.exists('robjects/fittedModels.rdata')) {
  load('robjects/fittedModels.rdata')
} else {
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
}

if (file.exists("robjects/gene.rdata")) {
  load("robjects/gene.rdata")
} else {
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
}


## ----Validation----------------------------------------------------------
if (file.exists("robjects/fittedModels.rdata")) {
  load("robjects/fittedModels.rdata")
} else {
  r2 <- lapply(names(pls), function(mdls){
      getValidation(pls[[mdls]], newdata = data[[mdls]][!data[[mdls]]$train, ])
  })
  names(r2) <- names(pls)  
}

## ----ValidationPlots, fig.height=4.5, results='hide', fig.cap="Maximum variation explained on response in Milk, NIR and Raman dataset and maximium correct classification proportion for gene Expression dataset. Each line represent a model for which the R2 predicted is calculated with number of components chosen by cross-validation (end of each line).", fig.pos='H'----
plt.list <- lapply(names(r2), function(x){
  plot(r2[[x]], full = FALSE) + ggtitle(toupper(x)) + theme(title = element_text(size = 7))
})
plot_grid(plotlist = plt.list, align = 'v', nrow = 2)


## ----classification, child='Classification.Rnw'--------------------------

## ----classify------------------------------------------------------------
gene$min.comp <- gene$r2$max.correct.prop[estimate == "CV", comp]
gene$test.pred <- getPredicted(gene$plda, newdata = data$gene$GeneExpr[!data$gene$train,], 
             ncomp = gene$min.comp, newY = data$gene$tumor[!data$gene$train])

## ----classifyPlot, message=FALSE-----------------------------------------
plt1 <- ggplot(data.table(gene$test.pred$test.pred$x, gene$test.pred$test.pred$class), aes(LD1, fill = V2, color = V2)) + 
  geom_density(alpha = 0.5) + scale_color_discrete(l = 40, name = "Test:") + 
  scale_fill_discrete(name = "Test:") + theme_bw(base_size = 18) +
  theme(legend.position = "top") + labs(x = 'LD1', y = 'Density')
plt2 <- ggplot(data.table(gene$test.pred$train.pred$x, gene$test.pred$train.pred$class), aes(LD1, fill = V2, color = V2)) + 
  geom_density(alpha = 0.5) + scale_color_discrete(l = 40, name = "Train:") + 
  scale_fill_discrete(name = "Train:") + theme_bw(base_size = 18) +
  theme(legend.position = "top") + labs(x = 'LD1', y = 'Density')
plt <- plot_grid(plt1, plt2, ncol = 1)
ggsave(file = "figure/classifPlot.pdf", plot = plt, width = 5, height = 7)


## ----variableSelection, child='VariableSelection.Rnw'--------------------

## ----Subset_Models-------------------------------------------------------
if (file.exists('robjects/SubModels.Rdata')) {
  load('robjects/SubModels.Rdata')
} else {
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
}

## ----VIPcowCompare, fig.height=2.3, fig.width='\\textwidth', fig.cap='Comparison of Models with complete set of variables and variables obtained from VIP filter'----
plt1 <- plotSeries(submodel$vip$milk$cow, pls$milk$cow, direction = 1, alpha = 0.75, ann.size = 3,
                     modelName = c("VIP", "Fullmodel")) + theme(legend.title = element_blank())
plt2 <- compareValidation(r2$milk, sub_r2$vip$milk, "cow", direction = 1) +
  theme(legend.title = element_blank())
grid_arrange_shared_legend(plt1, plt2, n.col = 2, width = c(4, 2))

## ----ShavingMCUV, fig.height=2.5, fig.cap="(left) Variables selected on MALDI-TOF Milk Spectra data where goat is a response. (right) Maximum $R^2$ Predicted from model with variables selected from Fullmodel, MCUV, Shaving method."----
plt1 <- compareValidation(r2$milk, 
                          subModel_r2df = list(sub_r2$mcuv$milk, sub_r2$shaving$milk), 
                          rsp = 'goat', direction = 1,
                          mdl.names = c('MCUV', 'Shaving', 'Full Model'))
plt2 <- plotSeries(submodel$shaving$milk$goat, pls$milk$goat, modelName = c("Shaving", "Fullmodel"), alpha = 0.85, ann.size = 2.5) + 
  geom_point(data = melt(as.table(loadings(submodel$mcuv$milk$goat)[, 1])), 
             aes(x = Var1, y = value, color = "MCUV"), size = 2)

grid_arrange_shared_legend(plt2, plt1, n.col = 2, width = c(4, 2))

## ----truncPlot, fig.height=3, fig.cap='Histogram of Loading Weights for full and truncated model for (left) MALDI-TOF milk data with response cow (right) Gene Expression data. Corresponding inset shows the maximum R-sq predicted for train, CV and test.', message=FALSE, warning=FALSE, fig.pos="H"----
trunc.cow <- truncPlot(pls$milk$cow, submodel$truncation$milk$cow, bins = 1000, 
                       newdata = data$milk[!data$milk$train, ], inset = T, invert = T)
trunc.gene <- truncPlot(gene$plda, gene$submodel$truncation, bins = 1000, inset = T, invert = F) + 
  coord_cartesian(xlim = c(-0.1, 0.35), ylim = c(0, 2500))
grid_arrange_shared_legend(trunc.cow, trunc.gene, n.col = 2, width = c(5, 5))


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

## ----NestedMixedModel----------------------------------------------------
test.r2 <- setnames(rbindlist(
  lapply(seq_along(dir("Reps")), function(idx){
  rd <- dir("Reps")[idx]
  load(file = file.path("Reps", rd), envir = .GlobalEnv)
  return(test.r2)
}), idcol = TRUE
), ".id", "replicates")
setnames(test.r2, "value", "rsq_pred")
test.r2[, c('method', 'data', 'response') := .(as.factor(method), as.factor(data), as.factor(response))]

## ----testTblPrint, results='asis', eval = FALSE--------------------------
test.xtbl <- xtable(test.r2[sample(nrow(test.r2), 5), .(replicates, method, data, response, model, rsq_pred)], align = "lrllXrr", 
                    caption = "Five random sample observation from dataset with $R^2$ predicted for test",
                    label = "tbl:testdf")
print(test.xtbl, table.placement = "H", caption.placement = "top", tabular.environment = "tabularx", width = '0.8\\textwidth', include.rownames = F)

## ----testPlot, fig.height=2.8, fig.cap="Boxplot for each data and method combination subdivided into their respective responses. Model from MCUVE in NIR and Raman data have drastically poor performance than other."----
ggplot(test.r2, aes(data, rsq_pred, group = response, fill = response)) + 
  geom_boxplot() + facet_wrap(~method, nrow = 1) + 
  theme(legend.position = "top", 
        legend.title = element_blank()) + 
  guides(fill = guide_legend(nrow = 1)) + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "R-sq Predicted")

## ----modelFitting, warning=FALSE, message=FALSE, results='hide'----------
test.r2[, log_rsq := log1p(rsq_pred * (1 - rsq_pred))]
mdl1 <- lm(log_rsq ~ method * r(data) + response %in% r(data), data = test.r2)
mdl2 <- lm(log_rsq ~ method * r(data) + response %in% r(data), data = test.r2, subset  = !(method == "mcuv" & data == "nir"))
mdl3 <- lm(log_rsq ~ method * r(data) + response %in% r(data), data = test.r2, subset  = method != "mcuv")

## ----anovaTbl, warning=FALSE, message=FALSE, results='hide'--------------
mdl <- list(Model1 = mdl1, Model2 = mdl2, Model3 = mdl3)
aovdf <- setnames(rbindlist(lapply(mdl, function(mdl){
  aov.tbl <- as.data.table(data.frame(Anova(mdl, type = "III")$anova), 
                         keep.rownames = T)
  setnames(aov.tbl, "Pr..F.", "p.value")
  return(na.omit(aov.tbl))
}), idcol = TRUE), ".id", "Model")
aovdf[, significance := ifelse(p.value <= 0.05, TRUE, FALSE)]

mdlsumry <- setnames(as.data.table(t(sapply(mdl, function(mdl) {
    list(sigma = summary(mdl)$sigma, rsq = summary(mdl)$r.squared, adj.rsq = summary(mdl)$adj.r.squared)
})), keep.rownames = T), "rn", "Model")

## ----rvfPlot, fig.height=2.5, warning=FALSE, message=FALSE, results='hide', fig.cap="Residuals vs Fitted Plot", fig.width='0.8\\textwidth', fig.pos='H'----
op <- par(mfrow = c(1, 3))
invisible(lapply(names(mdl), function(mdl.name) {
  plot(mdl[[mdl.name]], 1, main = mdl.name)
}))
par(op)

## ----anovaPlot, fig.height=3, fig.cap="Plot from Anova table. Number on top of each bars are p-value", fig.pos='H'----
ggplot(aovdf, aes(rn, F.value, fill = significance)) + 
  geom_bar(stat = "identity", na.rm = TRUE, position = "dodge") + 
  facet_grid(.~Model) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") + 
  geom_text(aes(label = round(p.value, 3)), nudge_y = 10, 
            family = "mono", size = rel(3), na.rm = T) + 
  scale_fill_brewer(palette = "Set1", direction = -1) +
  labs(x = NULL)


## ----codeUsed, child='CodeUsed.Rnw'--------------------------------------




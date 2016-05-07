# PCA model Fitting
pca$milk <- prcomp(data$milk$MS)
## PLS model for MALDI-TOF milk Dataset
pls$milk <- lapply(colnames(data$milk$milk), function(rsp) {
  plsr(milk[, rsp] ~ MS, data = data$milk[data$milk$train, ], validation = "CV", segments = 30, segment.type = "consecutive", jackknife = TRUE)
})
names(pls$milk) <- colnames(data$milk$milk)
## PLS with integrated LDA for GeneExpression
gene$plda <- plda(data$gene$GeneExpr[data$gene$train, ], cat.resp = data$gene$tumor[data$gene$train], 
                  fn = "plsr", fitComp = 25, split = 10)
## Getting Validation (R-sq predicted)
r2 <- lapply(names(pls), function(mdls) {
  getValidation(pls[[mdls]], newdata = data[[mdls]][!data[[mdls]]$train, ])
})
names(r2) <- names(pls)
# Variable Selection
subset <- subvaridx <- subdata <- submodel <- sub_r2 <- list()
## Filter with VIP
subset$vip <- lapply(names(pls), function(dta, vald = r2) {
  mdls <- lapply(names(pls[[dta]]), function(mdl) {
    opt.comp <- vald[[dta]][[2]][response == mdl & estimate == "CV", comp]
    vip <- VIP(pls[[dta]][[mdl]], opt.comp = opt.comp)
    vip.idx <- as.numeric(which(vip > 1))
    subdata <- within(data[[dta]], {assign(xy.info[[dta]][["x"]], data[[dta]][, xy.info[[dta]]["x"]][, vip.idx])})
    subMdl <- update(pls[[dta]][[mdl]], . ~ ., data = subdata[subdata$train, ])
    subvaridx$vip[[dta]][[mdl]] <<- vip.idx
    subdata$vip[[dta]][[mdl]] <<- subdata
    submodel$vip[[dta]][[mdl]] <<- subMdl
    return(vip)})
  names(mdls) <- names(pls[[dta]])
  sub_r2$vip[[dta]] <<- getValidation(submodel$vip[[dta]], subdata$vip[[dta]])
  return(mdls)
})
names(subset$vip) <- names(pls)
# Model Fitting for model comparison with data as random factor and response nested on it
mdl <- lm(value ~ method * r(data) + response %in% r(data), data = test.r2)
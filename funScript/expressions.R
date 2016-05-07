## Prepare Data.Table
prepInit <- quote({
  dataInfo <- {rbind(
    expand.grid(dataName = "MALDITOF_Milk", x = 'MS', y = 'milk', 
                cv.seg = 30, seg.type = 'consecutive', jack.knife = TRUE,
                scale = FALSE, stringsAsFactors = F),
    expand.grid(dataName = "GeneExpr_Prostate", x = 'GeneExpr', y = 'tumor', 
                cv.seg = 10, seg.type = 'random', jack.knife = TRUE,
                scale = FALSE, stringsAsFactors = F),
    expand.grid(dataName = "NIR_Raman_PUFA", x = c('NIR', 'Raman'), y = 'PUFA', 
                cv.seg = 10, seg.type = 'random', jack.knife = TRUE,
                scale = FALSE, stringsAsFactors = F)
  )}
  modelInfo <- {expand.grid(
    model = c("pcr", "plsr"),
    dataName = c('MALDITOF_Milk', 
                 'GeneExpr_Prostate', 
                 'NIR_Raman_PUFA'),
    stringsAsFactors = F
  )}
  
  info <- data.table(merge(dataInfo, modelInfo, by = "dataName"))
  info[, formula := paste0(y, '~', x)]
  info[, train := lapply(dataName, function(x) get(x)$train)]
  rm(dataInfo, modelInfo)
})

## Principal Component Analysis
pca.fit <- quote(
  info[, pca.fit := pmap(.(dataName, scale), function(db, sc){
    prcomp(get(db), center = sc)
  })]
)

## Model Fitting (PCR and PLS)
model.fit <- quote(
  info[, model.fit := 
         pmap(.(formula, model, dataName, cv.seg, seg.type, jack.knife, scale, train), 
              function(f, mdl, db, cv, stype, jk, sc, trn){
                match.fun(mdl)(as.formula(f), 
                               data = get(db)[trn,],
                               scale = sc,
                               validation = "CV",
                               segments = cv,
                               segment.type = stype,
                               jackknife = jk)
              })
       ]
)
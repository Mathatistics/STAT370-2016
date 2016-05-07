## Used Functions ---------------

grid_arrange_shared_legend <- function(..., n.col = 1, width = NULL) {
  plots <- list(...)
  legend <- gtable::gtable_filter(ggplotGrob(
    plots[[1]] + theme(legend.position = "top")), "guide-box")
  lheight <- sum(legend$height)
  
  plt.lst <- lapply(plots, function(x) x + theme(legend.position = "none"))
  plt.lst$widths <- width
  plt.tbl <- do.call(gridExtra::arrangeGrob, plt.lst)
  
  gridExtra::grid.arrange(legend, plt.tbl, ncol = 1,
    heights = grid::unit.c(lheight, grid::unit(1, "npc") - lheight))
}

plda <- function(x, cat.resp, fn = 'plsr', split = 10, 
                 fitComp = min(dim(x)), prior = F, ldafn = 'lda', ...){
  fn <- match.fun(fn)
  
  train.err.mat <- matrix(ncol = split, nrow = fitComp, 
                          dimnames = list(Comp = 1:fitComp,
                                          Split = 1:split))
  test.err.mat <- matrix(ncol = split, nrow = fitComp, 
                         dimnames = list(Comp = 1:fitComp,
                                         Split = 1:split))

  y.mdl.mtrx <- model.matrix(~cat.resp - 1)
  # y.mdl.mtrx <- apply(y.mdl.mtrx, 2, as.factor)
  myData <- data.frame(y = I(y.mdl.mtrx), x = I(x))
  cvSplits <- cvsegments(nrow(myData), k = split, type = 'random')
  
  ## Start cv-loop
  for (k in seq_along(cvSplits)) {
    splt <- cvSplits[[k]] 
    ## Fitting model on whole dataset
    
    if (identical(fn, match.fun("truncation"))) {
      mdl <- fn(y ~ x, data = myData[-splt, ], ncomp = fitComp, 
                truncation = "Lenth", trunc.width = 0.95)
    } else {
      mdl <- fn(y ~ x, data = myData[-splt, ], ncomp = fitComp, ...)
    }
    
    x.scrs <- mdl$scores
    ## Predict test Scores
    x.tst.scrs <- predict(mdl, newdata = myData[splt, ], 
                          type = 'scores')
    for (cmp in 1:fitComp) {
      ## Predict X-Scores to use as predictor in LDA
      ldaData.train <- data.frame(y = cat.resp[-splt], 
                                  x = I(x.scrs[, 1:cmp, drop = F]))
      
      ldaData.test <- data.frame(y = cat.resp[splt], 
                                 x = I(x.tst.scrs[, 1:cmp, drop = F]))

      if (!prior) {
        prior.prob <- as.numeric(1/table(factor(ldaData.train$y)))
        prior.prob <- prior.prob/sum(prior.prob)
      } else {
        prior.prob <- prior
      }

      ## LDA fit
      if (ldafn != 'lda') {
        lda.fit <- klaR::rda(y ~ x, data = ldaData.train, gamma = 0, lambda = 1)
      } else {
        lda.fit <- MASS::lda(y ~ x, data = ldaData.train)
      }
   
      ## LDA Prediction
      lda.trn.pred <- predict(lda.fit)
      lda.tst.pred <- predict(lda.fit, newdata = ldaData.test)
      
      ## Confusion Table
      trn.conf.tbl <- table(lda.trn.pred$class, ldaData.train$y)
      tst.conf.tbl <- table(lda.tst.pred$class, ldaData.test$y)
      
      ## Prediction Error
      tst.pred.err <- 1 - sum(diag(tst.conf.tbl)) / sum(tst.conf.tbl)
      trn.pred.err <- 1 - sum(diag(trn.conf.tbl)) / sum(trn.conf.tbl)
      
      ## Saving Results
      test.err.mat[cmp, k] <- tst.pred.err
      train.err.mat[cmp, k] <- trn.pred.err
    }
  }
  
  ## We have got a full fitted mvr model
  ## Need minimum component required to get minimum test error
  if (identical(fn, match.fun("truncation"))) {
    mdl <- fn(y ~ x, data = myData, ncomp = fitComp,
              truncation = "Lenth", trunc.width = 0.95)
  } else {
    mdl <- fn(y ~ x, data = myData, ncomp = fitComp, ...)
  }
  
  mdl$cat.resp <- cat.resp
  mdl$test.err <- test.err.mat
  mdl$train.err <- train.err.mat
  mdl$lda.fun <- ldafn
  class(mdl) <- append(class(mdl), 'plda')
  
  return(mdl)
}

## Classification method for cv
classify <- function(x, ncomp, ldafn) UseMethod('classify', x)
classify.plda <- function(x, ncomp){
  if (!any(class(x) %in% 'plda'))
    stop('Only works for objects with class plda')
  
  ldaData <- data.frame(x = I(x$scores[, 1:ncomp]), 
                        y = x$cat.resp)
  
  if (x$lda.fun != 'lda') {
    lda.fit <- klaR::rda(y ~ x, data = ldaData, gamma = 0, lambda = 1)
  } else {
    lda.fit <- MASS::lda(y ~ x, data = ldaData)
  }
  
  
  lda.fit$cat.resp <- ldaData$y
  
  class(lda.fit) <- append(class(lda.fit), 'cplda')
  return(lda.fit)
}

getPredicted <- function(x, newdata, ncomp, newY) UseMethod('getPredicted', x)
getPredicted.plda <- function(x, ncomp, newdata = NULL, newY = NULL){
  if (!any(class(x) %in% 'plda'))
    stop('Only works for objects with class plda')
  ret <- list()
  
  ldaModel <- classify(x, ncomp = ncomp)

  if (!is.null(newdata)) {
    ldaNewData <- predict(x, 
                          newdata = data.frame(x = I(as.matrix(newdata))), 
                          type = 'scores')[, 1:ncomp]
    
    ret$test.pred <- predict(ldaModel, newdata = data.frame(x = I(ldaNewData)))
  }
  
  ret$train.pred <- predict(ldaModel)
  
  ## Errors
  if (!is.null(newY)) {
    ret$test.err <- table(Original = newY, Predicted = ret$test.pred$class)
  }
  
  ret$train.err <- table(Original = x$cat.resp, Predicted = ret$train.pred$class)
  ret$ncomp <- ncomp
  class(ret) <- 'plda.pred'

  return(ret)
}

msc <- function(x) UseMethod('msc', x)
msc.plda.pred <- function(x){
  if (!any(class(x) %in% 'plda.pred'))
    stop('Only works for objects with class plda.pred')
  
  require(data.table)
  
  err.lst <- list()
  if (exists('test.err', x)) {
    test.err.tbl <- x$test.err
    test.msc <- 1 - sum(diag(test.err.tbl))/sum(test.err.tbl)
    err.lst$test <- test.err.tbl
  }
  train.err.tbl <- x$train.err
  train.msc <- 1 - sum(diag(train.err.tbl))/sum(train.err.tbl)
  err.lst$train <- train.err.tbl
  
  # err.lst <- list(test = test.err.tbl, train = train.err.tbl)
  
  conf.dt <- data.table(reshape2::melt(err.lst))

  conf.dt <- conf.dt[, Correct := as.character(Original) == as.character(Predicted), 
                     by = L1]
  conf.dt <- conf.dt[, .(value = sum(value)), by = .(Original, Correct, L1)]

  if (exists('test.err', x))
    err.lst$test.msc <- test.msc
  
  err.lst$train.msc <- train.msc
  err.lst$conf.dt <- conf.dt
  setkeyv(err.lst$conf.dt, c('Original', 'Correct', 'L1'))
  err.lst$ncomp <- x$ncomp
  class(err.lst) <- 'plda.err'
  
  return(err.lst)
}

plot.plda.err <- function(x){
  dt <- x$conf.dt
  require(ggplot2)
  plt <- ggplot(dt, aes(Original, value)) + 
    geom_bar(stat = 'identity',
             position = 'fill',
             aes(fill = Correct)) +
    theme_bw() +
    theme(legend.position = 'top',
          axis.title = element_blank()) +
    guides(fill = guide_legend(title = 'Correct Predictions',
                               title.position = 'top', 
                               title.hjust = 0.5)) +
    coord_flip() +
    facet_grid(.~L1, as.table = TRUE) + 
    geom_text(aes(ymax = value/sum(value), 
                  label = value,
                  ymin = 0), 
              position = 'fill',  
              hjust = 1, size = 3.5)
  return(plt)
  }

getError <- function(x, type, plot) UseMethod('getError', x)
getError.plda <- function(x, type = c('test', 'train', 'both'), plot = FALSE){
  ret <- list()
  if (any(type %in% c('test', 'both'))) {
    ret$cvtest <- x$test.err
  }
  if (any(type %in% c('train', 'both'))) {
    ret$train <- x$train.err
  }
  
  if (plot) {
    require(data.table)
    require(ggplot2)
    require(reshape2)
    dt <- data.table(melt(ret))[, .(value = mean(value)), by = .(Comp, L1)]
    plt <- ggplot(dt, aes(Comp, value))
    if (type == 'both') {
      plt <- plt + geom_line(aes(color = L1))
    } else {
      plt <- plt + geom_line()
    }
    # plt <- plt + geom_point(shape = 21, fill = 'gray', size = 0.75)
    plt <- plt + theme_bw()
    plt <- plt + theme(legend.title = element_blank(),
                       legend.position = 'top')
    plt <- plt + labs(x = 'Components', y = 'Misclassification Error')
    plt <- plt + annotate('text', x = Inf, y = Inf, 
                          label = paste('Comp:', dt[which.min(value), Comp],
                                         '\nError:', dt[which.min(value), round(value, 2)]),
                          vjust = 1.5, hjust = 1.5)
    ret$plot <- plt
  }
  require(data.table)
  ret$err.dt <- data.table(melt(ret[-3]))
  ret$err.dt <- ret$err.dt[, .(value = mean(value)), by = .(L1, Comp)]
  setkeyv(ret$err.dt, c('L1', 'Comp'))
  return(invisible(ret))
}

plotScore <- function(model, title, cat.var = train.1$y, point.size = 1.8, comps = 1:2){
  scores <- data.table(model$scores[, 1:3])[, event := cat.var]
  e <- parent.env(environment())
  names(scores) <- gsub('[[:space:]]', '', names(scores))
  which.comps <- names(scores)[comps]
  # browser()
  plt <- ggplot(scores, aes_string(which.comps[1],
                            which.comps[2], 
                            color = 'event'),
                environment = e) + 
    geom_point(size = point.size) +
    theme_bw() +
    theme(legend.position = 'top',
          legend.direction = 'horizontal') +
    guides(color = guide_legend(nrow = 3, title.position = 'top',
                                title.hjust = 0.5)) +
    labs(x = paste(which.comps[1], ' (', round(explvar(model)[comps[1]], 2), 
                   '%)', sep = ''),
         y = paste(which.comps[2], ' (', round(explvar(model)[comps[2]], 2), 
                   '%)', sep = '')) +
    ggtitle(title)
}

mscPlot <- function(mscdt, type){
  # Misclassification Plot
  plt <- ggplot(mscdt[L1 == type], aes(Original, value)) + 
    geom_bar(stat = 'identity',
             position = 'fill',
             aes(fill = Correct)) +
    theme_bw() +
    theme(legend.position = 'top',
          axis.title = element_blank(),
          text = element_text(size = 8)) +
    guides(fill = guide_legend(title = 'Correct Predictions',
                               title.position = 'top', 
                               title.hjust = 0.5)) +
    coord_flip()
  if ("variable" %in% names(mscdt)) {
    plt <- plt + facet_grid(. ~ variable, as.table = TRUE, drop = T)
  }
     plt <- plt + geom_text(aes(ymax = value/sum(value), 
                  label = value,
                  ymin = 0), 
              position = 'fill',  
              hjust = 1, size = 2)
}
splitMe <- function(df, prob, seed = 1000){
  set.seed(seed)
  df$train <- as.logical(rbinom(nrow(df), size = 1, prob = prob))
  return(df)
}

load_or_run <- function(file, expr.list) {
  if (!file.exists(file)) {
    invisible(lapply(expr.list, function(expr){
      eval(expr, envir = .GlobalEnv)
    }))
  } else {
    load(file, envir = .GlobalEnv)
  }
}

plotPCA <- function(pca.obj, label = F, npcs = 10, title = NULL, base.size = 10){
  require(ggplot2)
  comps <- factor(colnames(pca.obj$x), levels = colnames(pca.obj$x))
  comps <- as.numeric(gsub('[[:alpha:][:blank:]]', '', comps))
  sdev.df <- data.frame(comps = comps,
                        sdev = pca.obj$sdev,
                        explvar = pls::explvar(pca.obj))
  sdev.df <- sdev.df[1:npcs, ]
  plt <- sdev.df %>% 
    ggplot(aes(comps, sdev, group = 1)) +
    geom_line() + geom_point(size = 0.5)
  if (label) {
    plt <- plt + geom_text(aes(label = round(explvar, 2)), 
                           nudge_y = 0.5, nudge_x = 0.15)
  }
  plt <- plt + theme_bw(base_size = base.size) +
    labs(x = "Components", y = "Standard Deviation") +
    ggtitle(title)
  return(plt)
}

getValidation <- function(model, newdata, type = 'R2') {
  if ("mvr" %in% class(model)) model <- list(model)
  if (class(newdata) == "list") {
    if (length(newdata) != length(model))
      stop("Model and New Data should have same length or 
           New Data should be of length 1")
  }
  name_seq <- if (is.null(names(model))) seq_along(model) else names(model)
  
  vald.list <- llply(name_seq, function(x) {
    if (class(newdata) == "list" & length(newdata) > 1) 
      newdata <- newdata[[x]]
      newdata <- newdata[!newdata$train, ]
    
    rmsep <- match.fun(type)(model[[x]],  estimate = "all", 
                   newdata = newdata)$val
    df <- data.table(melt(rmsep))[, -2, with = FALSE]
    df[, comp := as.integer(ifelse(grepl('^\\(Int', model), 0, 
                                   gsub('[[:alpha:]]', '', model)))]
    return(df)
  })
  names(vald.list) <- name_seq
  vald <- rbindlist(vald.list, idcol = TRUE)
  
  names(vald)[1] <- "response"
  opt.comp <- vald[estimate == "CV", 
                   .(comp = ifelse(type == "R2", comp[which.max(value)], 
                                   comp[which.min(value)])), 
                   by = response]
  opt.vald <- merge(opt.comp, vald, by = c('comp', 'response'))
  
  ret <- list(vald, opt.vald, type)
  names(ret) <- c(tolower(type), 
                  paste0(ifelse(type == "R2", 'max.', 'min.'), tolower(type)),
                  'type')
  class(ret) <- 'validation'
  return(ret)
}

plot.validation <- function(valdObj, full = TRUE, ncomp = 20) {
  if (full) {
    qplot(comp, value, color = estimate,
          data = valdObj[[1]][comp <= ncomp], geom = c('line'), 
          xlab = 'Components', ylab = valdObj[[3]]) + 
      geom_point(size = 0.5) +
      facet_wrap(~response, scales = 'free_y') +
      theme(legend.position = 'top', legend.title = element_blank())
  } else {
    ggplot(valdObj[[2]], aes(estimate, value, color = response, group = response)) + 
      geom_line() + geom_point() + 
      geom_text(data = valdObj[[2]][estimate == "test"], show.legend = FALSE,
                aes(x = "test", y = value, label = comp), nudge_x = 0.25) +
      scale_color_brewer(palette = 'Set1') +
      theme(axis.title.x = element_blank(), 
            legend.position = if (length(unique(valdObj[[2]]$response)) > 1) 
              "top" else "none",
            legend.title = element_blank()) +
      ylab(paste(ifelse(valdObj[[3]] %in% c("R2", "Classification Prop."), 
                        "Max.", "Min."), valdObj[[3]]))
  }
}

prepData <- function(model, fn, ncomp) {
  require(data.table)
  require(reshape2)
  
  if (fn == 'coef') {
    dta <- drop(match.fun(fn)(model, ncomp = ncomp))
    if (any(class(model) == "plda")) dta <- dta[,1]
    dta <- data.table(melt(dta), keep.rownames = TRUE)
  } else {
    dta <- data.table(match.fun(fn)(model)[, ncomp, drop = FALSE], 
                      keep.rownames = TRUE)
  }
  dta <- dta[, idx := .I]
  dta <- dta[, rn := as.numeric(gsub('[[:alpha:]]', '', rn))]
  setnames(dta, names(dta), c('Variables', fn, 'idx'))
  setkey(dta, "Variables")
  return(dta)
}

plotSeries <- function(fullModel, subModel = NULL, 
                       type = "loading.weights", invert = FALSE, 
                       ncomp = 3, line = TRUE, na.rm = FALSE,
                       modelName = c("Full-Model", "Sub-Model"),
                       direction = 1, alpha = 0.75, title = FALSE,
                       ann.size = 3) {
  
  require(pls)
  require(data.table)
  require(ggplot2)
  
  if (!type %in% c("loadings", "loading.weights", "coef")) stop("Undefined Type")
  full_dta <- prepData(fullModel, type, ncomp)
  if (is.null(subModel)) {
    plt <- ggplot(full_dta, aes_string('Variables', type))
    plt <- plt + geom_line()
  }
  if (!is.null(subModel)) {
    sub_dta <- prepData(subModel, type, ncomp)
    if (invert) sub_dta[[type]] <- -sub_dta[[type]]
    unique_keys <- unique(c(full_dta$Variables, sub_dta$Variables))
    
    dta <- sub_dta[full_dta[J(unique_keys)]][,-c(3,5), with = F]
    setnames(dta, names(dta)[-1], c(modelName[2],  modelName[1]))
    dta <- melt(dta, 1)
    setnames(dta, c("variable", "value"), c("Model", type))
    
    if (na.rm) dta <- na.omit(dta)
    
    anotTbl <- na.omit(dta)[, .N, by = Model]
    anotText <- paste(apply(anotTbl, 1, paste, collapse = ':'), collapse = '\n')
    
    
    plt <- ggplot(dta, aes_string('Variables', type))
    if (line) {
      plt <- plt + geom_line(aes(color = Model), na.rm = TRUE, alpha = alpha)
    } else {
      plt <- plt + geom_point(aes(color = Model), na.rm = TRUE, size = 0.5, alpha = alpha)
    }
    plt <- plt + scale_color_brewer(palette = 'Set1', direction = direction)
    plt <- plt + theme(legend.position = "top")
    if (title) plt <- plt + 
      ggtitle(paste('Plot for', type, 'of', paste(modelName, collapse = " and ")))
    plt <- plt + annotate(geom = 'text', Inf, Inf, 
                          label = paste(anotText), hjust = 1.2, vjust = 1.5,
                          size = ann.size)
  }
  return(plt)
}

getSeries <- function(dMatrix, nsamp, xlab = NULL, ylab = NULL, title = NULL){
  df <- data.table(melt(unclass(dMatrix)
                        [sample(1:nrow(dMatrix), nsamp), ,drop = F]))
  ggplot(df, aes(Var2, value, color = as.factor(Var1))) + 
    geom_line(alpha = 0.75) + 
    theme(legend.position = "none") + 
    scale_color_brewer(palette = 'Spectral') + 
    labs(x = xlab, y = ylab) +
    ggtitle(title)
}

truncPlot <- function(fullModel, truncModel, ncomp = 1, 
                      bins = ncol(fullModel$model[[2]]), 
                      newdata = NULL, invert = FALSE, 
                      inset = FALSE, limit = NULL) {
  fullLwt <- melt(as.table(fullModel$loading.weights[, ncomp, drop = F]))
  truncLwt <- melt(as.table(truncModel$loading.weights[, ncomp, drop = F]))
  
  binObj <- qplot(value, data = fullLwt, bins = bins)
  
  if (is.null(limit)) limit <- layer_scales(binObj)$y$range$range[2] * 1.1

  if (invert) truncLwt <- within(truncLwt, value <- -value)
  dta <- setnames(rbindlist(list(
    fullModel = fullLwt, 
    truncModel = truncLwt
  ), idcol = TRUE), 
  c('.id', 'Var1', 'Var2', 'value'), 
  c('model.type', 'x.var', 'comp', 'loading.weights'))
  
  if (inset) {
    if ("plda" %in% class(truncModel)) {
      mdl.name <- as.character(substitute(truncModel)[[2]][[2]])
      plt.inset <- compareValidation(get(mdl.name)[['r2']], get(mdl.name)[['sub_r2']][['truncation']], rsp = 'tumor')
    } else {
      mdl.name <- as.character(substitute(truncModel)[[2]][[3]])
      rsp.name <- as.character(substitute(truncModel)[[3]])
      plt.inset <- compareValidation(r2[[mdl.name]], sub_r2$truncation[[mdl.name]], rsp = rsp.name)
    }
    plt.inset <- plt.inset + theme(text = element_text(size = rel(3)), axis.title = element_blank())
    plt.inset <- ggplotGrob(plt.inset)[-c(2:3),]
  }
  plt <- ggplot(dta, aes(loading.weights, fill = model.type)) + 
    geom_histogram(bins = bins, position = "identity") + 
    facet_grid(.~comp, scales = 'free', space = 'free') + 
    labs(x = 'Loadings Weights', y = 'Frequency') +
    theme(legend.title = element_blank(), legend.position = "top") +
    scale_fill_brewer(palette = 'Set1') +
    coord_cartesian(ylim = c(0, limit))
  
  if (inset) {
  xlims <- ggplot2::layer_scales(binObj)[[1]]$range$range
  ylims <- ggplot2::layer_scales(binObj)[[2]]$range$range
    plt <- plt + annotation_custom(
      plt.inset,
      xmin = xlims[2] * 0.2, 
      xmax = xlims[2],
      ymin = ylims[2] * 0.5, 
      ymax = ylims[2])
  }
  return(plt)
}

classR2 <- function(plda.mdl, newX, newY, ncomp) {
  vald <- rbindlist(lapply(1:ncomp, function(comp){
    pred <- getPredicted(plda.mdl, newdata = newX, ncomp = comp, newY = newY)
    r2 <- data.table(
      response = "tumor",
      comp = comp,
      train = 1 - msc(pred)$train.msc,
      CV = apply(1 - plda.mdl$test.err, 1, mean)[comp],
      test = 1 - msc(pred)$test.msc
    )
    ret <- melt(data.table(r2), 1:2, variable.name = 'estimate')
    ret[, model := paste(comp, "comps")]
    ret <- ret[, .(response, comp, estimate, model, value = as.numeric(value))]
    return(ret)
  }))
  opt.comp <- vald[estimate == "CV", 
                   .(comp = comp[which.max(value)]), 
                   by = response]
  opt.vald <- merge(opt.comp, vald, by = c('comp', 'response'))
  
  ret <- list(vald, opt.vald, type = "Classification Proportion")
  names(ret) <- c("correct.prop", "max.correct.prop", 'type')
  class(ret) <- 'validation'
  return(ret)
}

compareValidation <- function(fullModel_r2df, subModel_r2df, rsp, mdl.names = NULL, direction = 1) {
  if (!rsp %in% unique(fullModel_r2df[[2]]$response)) stop("Input Correct Response")
  if ('validation' %in% class(subModel_r2df)) subModel_r2df <- list(subModel_r2df)
  if (is.null(mdl.names)) 
    mdl.names <- c(paste("SubModel", 1:length(subModel_r2df), sep = "_"), "FullModel")
  
  df.list <- subModel_r2df
  df.list$fullModel <- fullModel_r2df
  df.list <- lapply(df.list, function(x) x[[2]])
  names(df.list) <- mdl.names
  dt <- setnames(rbindlist(df.list, idcol = TRUE), '.id', 'Model')
  
  plt <- ggplot(dt[response == rsp], aes(estimate, value, color = Model, group = Model)) + 
    geom_point() + geom_line() +
    theme(legend.position = "top", legend.title = element_blank()) +
    labs(x = NULL, y = "Maximum R-sq Predicted") +
    scale_color_brewer(palette = 'Set1', direction = direction)
  return(plt)
}
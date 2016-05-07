pca <- function(X, scale = TRUE, center = TRUE, ncomp = min(dim(X))){
  ## Scaling the X object
  X <- scale(X, center = center, scale = scale)
  
  ## Get names
  nvar <- paste("Comp", 1:min(ncol(X), nrow(X)), sep = ".")
  comp.names <- nvar[1:ncomp]
  colNames <- colnames(X)
  rowNames <- rownames(X)
  
  ## Singular Value Decomposition
  udv <- svd(X, nu = ncomp, nv = ncomp)
  
  
  ## Get Scores
  scores <- udv$u %*% diag(udv$d[1:ncomp])
  rownames(scores) <- rowNames
  colnames(scores) <- comp.names
  
  ## Get Loadings
  loadings <- udv$v
  rownames(loadings) <- colNames
  colnames(loadings) <- comp.names
  
  
  ## Compute Std.Dev
  sdev <- udv$d / (nrow(X) - 1)
  names(sdev) <- nvar
  
  ## Explained Variation
  explvar <- sdev ^ 2 / sum(sdev ^ 2) * 100
  names(explvar) <- comp.names
  
  ret <- list(
    scores = scores,
    loadings = loadings,
    sdev = sdev,
    explvar = explvar,
    comps = nvar,
    ncomps = comp.names,
    data = X
  )
  class(ret) <- c("pca", "mvr")
  return(ret)
}
model.matrix.pca <- function(object)object$data

## ----summaryMethod, echo=FALSE-------------------------------------------
summary.pca <- function(pca.obj){
  ret <- rbind(pca.obj$sdev, pca.obj$explvar, cumsum(pca.obj$explvar))
  rownames(ret) <- c("Std.Dev", "Explained Variation", "Cummulative Exp.Var")
  colnames(ret) <- pca.obj$comps
  return(round(ret, 3))
}

## printMethod ------------------------------------------------------------
print.pca <- function(pca.obj){
  print(summary(pca.obj))
}

## ----plotMethod----------------------------------------------------------
plot.pca <- function(pca.obj, label = F){
  require(ggplot2)
  comps <- factor(pca.obj$comps, levels = pca.obj$comps)
  sdev.df <- data.frame(comps = comps,
                        sdev = pca.obj$sdev,
                        explvar = pca.obj$explvar)
  plt <- sdev.df %>% 
    ggplot(aes(comps, sdev, group = 1)) +
    geom_line() + geom_point()
  if (label) {
    plt <- plt + geom_text(aes(label = round(explvar, 2)), 
                           nudge_y = 0.5, nudge_x = 0.15)
  }
  plt <- plt + theme_bw() +
    labs(x = "Components", y = "Standard Deviation") +
    ggtitle("Scree Plot for\nPrincipal Component Analysis")
  return(plt)
}

## ----scorePlotMethod, echo = FALSE---------------------------------------
scorePlot <- function(x, ...) {
  UseMethod("scorePlot", x)
}
scorePlot.pca <- function(pca.obj, comps = 1:2, useNames = FALSE){
  scores <- as.data.frame(pca.obj$scores[ , comps])
  vars <- names(scores)
  explvar <- pca.obj$explvar[comps]
  labs <- paste(vars, "\n(", round(explvar, 2), "%)", sep = "")
  if (length(comps) > 2) {
    require(GGally)
    
    diag.list <- lapply(seq_along(labs), function(x){
      ggplot(expand.grid(1, 1)) + 
        geom_text(label = labs[x], x = 0.5, y = 0.5, 
                  size = 6) +
        theme(panel.background = element_rect(fill = NA))
    })
    
    plt <- ggpairs(scores, 
                   upper = list(continuous = "points", combo = "dot"),
                   lower = list(continuous = "points", combo = "dot"), 
                   columnLabels = rep("", length(labs))) + theme_bw()
    
    for (x in seq_along(labs)) {
      plt <- putPlot(plt, diag.list[[x]], x, x)
    }
    
    for (col in 1:plt$ncol) {
      for (row in 1:plt$nrow) {
        if (col == row) {
          plt[row, col] <- plt[row, col]
        } else {
          if (useNames) {
            plt[row, col] <- plt[row, col] +
              geom_text(label = rownames(scores), size = 3)
            plt[row, col]$layers <- plt[row, col]$layers[-1]
          }
          plt[row, col] <- plt[row, col] + 
            geom_vline(xintercept = 0, col = "blue", linetype = 2) +
            geom_hline(yintercept = 0, col = "blue", linetype = 2)
        }
      }
    }
  } else {
    require(ggplot2)
    labs <- gsub("\n", " ", labs)
    plt <- scores %>% 
      ggplot(aes_string(vars[1], vars[2]))
    if (useNames) {
      plt <- plt + geom_text(label = rownames(scores), size = 3)
    } else {
      plt <- plt + geom_point()
    }
    plt <- plt + theme_bw() +
      labs(x = labs[1], y = labs[2]) +
      geom_vline(xintercept = 0, col = "blue", linetype = 2) +
      geom_hline(yintercept = 0, col = "blue", linetype = 2) +
      ggtitle("Score Plot of\nPrincipal Component Analysis")
  }
  return(plt)
}

## ----loadingPlotMethod, echo = FALSE-------------------------------------
loadingPlot <- function(x, ...) {
  UseMethod("loadingPlot", x)
}
loadingPlot.pca <- function(pca.obj, comps = 1:2, scatter = FALSE, useNames = FALSE){
  loads <- as.data.frame(pca.obj$loadings[ , comps])
  comps <- names(loads)
  vars <- rownames(loads)
  explvar <- pca.obj$explvar[comps]
  labs <- paste(comps, "\n(", round(explvar, 2), "%)", sep = "")
  
  if (!scatter) {
    require(ggplot2)
    plt <- loads %>% 
      dplyr::add_rownames("variables") %>% 
      tidyr::gather("Comp", "value", -variables) %>% 
      ggplot(aes(variables, value, group = Comp)) +
      geom_line(aes(linetype = Comp)) +
      theme_bw() + theme(legend.position = "top") +
      labs(x = "Variables", y = "Loadings") +
      ggtitle("Loading Plot for\nPrincipal Component Analysis")
  } else {
    if (length(comps) > 2) {
      require(GGally)
      
      diag.list <- lapply(seq_along(labs), function(x){
        ggplot(expand.grid(1, 1)) + 
          geom_text(label = labs[x], x = 0.5, y = 0.5, 
                    size = 6) +
          theme(panel.background = element_rect(fill = NA))
      })
      
      plt <- ggpairs(loads, 
                     upper = list(continuous = "points", combo = "dot"),
                     lower = list(continuous = "points", combo = "dot"), 
                     columnLabels = rep("", length(labs))) + theme_bw()
      
      for (x in seq_along(labs)) {
        plt <- putPlot(plt, diag.list[[x]], x, x)
      }
      
      for (col in 1:plt$ncol) {
        for (row in 1:plt$nrow) {
          if (col == row) {
            plt[row, col] <- plt[row, col]
          } else {
            if (useNames) {
              plt[row, col] <- plt[row, col] +
                geom_text(label = vars, size = 3)
              plt[row, col]$layers <- plt[row, col]$layers[-1]
            }
            plt[row, col] <- plt[row, col] + 
              geom_vline(xintercept = 0, col = "blue", linetype = 2) +
              geom_hline(yintercept = 0, col = "blue", linetype = 2)
          }
        }
      }
    } else {
      require(ggplot2)
      labs <- gsub("\n", " ", labs)
      
      plt <- loads %>% 
        ggplot(aes_string(comps[1], comps[2]))
      if (useNames) {
        plt <- plt + geom_text(label = vars, size = 3)
      } else {
        plt <- plt + geom_point()
      }
      plt <- plt + theme_bw() +
        labs(x = labs[1], y = labs[2]) +
        geom_vline(xintercept = 0, col = "blue", linetype = 2) +
        geom_hline(yintercept = 0, col = "blue", linetype = 2) +
        ggtitle("Score Plot of\nPrincipal Component Analysis")
    }
  }
  return(plt)
}
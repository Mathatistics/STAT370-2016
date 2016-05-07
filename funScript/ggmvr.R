ggscore <- function(x, ...) {
  UseMethod("ggscore", x)
}

ggscore.mvr <- function(mvr.obj, comps = 1:2, useNames = FALSE){
  scores <- as.data.frame(mvr.obj$scores[ , comps])
  names(scores) <- gsub("[[:space:]]", "", names(scores))
  vars <- names(scores)
  explvar <- explvar(mvr.obj)[comps]
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
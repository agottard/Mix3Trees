#' Plot a tree of mixed effects 3Trees models.
#'
#' @param mod3T The output of fit3trees.
#' @param tree Tree indicator: Must be 1, 2, or 3.
#' @return A plot for one of the tree components
#'
#' @examples
#' data(mydat)
#' X<-mydat[,1:5]
#' mod3T <- fit3Trees(Y=mydat$Y, X=X, gr=mydat$gr, covLin=names(X), 
#'                  covT1=names(X[,1:3]), covT2=names(X[,4:5]), 
#'                  covT3=colnames(X),
#'                  niter = 50, re_form = "(1|gr)",
#'                  selective.inference = FALSE)
#'                  
#' plot3T(mod3T, tree = 1)
#' plot3T(mod3T, tree = 2)
#' plot3T(mod3T, tree = 3)
#'
#' @import rpart rpart.plot
#' @export

plot3T <- function(mod3T, tree = 1){
  requireNamespace("rpart", quietly = TRUE)
  requireNamespace("rpart.plot", quietly = TRUE)

  ## Function for getting the tree  
  get_tree <- function(obj, tree) {
    switch(tree,
           `1` = obj$Tree1,
           `2` = obj$Tree2,
           `3` = obj$Tree3,
           stop("Invalid tree number. Must be 1, 2, or 3.")
    )
  }
  # Get the tree
  myTree <- get_tree(mod3T, tree)
  # Get the estimates
  nome <- paste0("^", "Tree", tree)
  myrow <- grep(nome, rownames(summary(mod3T$mod.final)$coefficients))  
  myval <- round(summary(mod3T$mod.final)$coefficients[myrow,1] ,4)
  ncoef <- length(myrow)
  nleaf <- length(myTree$frame$yval[myTree$frame$var == "<leaf>"])
  if (nleaf > ncoef) myval <- c(rep(0,(nleaf-ncoef)),myval)
  
  # Replace values
  
  myTree$frame$yval[myTree$frame$var == "<leaf>"]<- myval 
  
  # Plot the rpart tree 
  rpart.plot(myTree, roundint=FALSE, type=0, extra = 1)
}


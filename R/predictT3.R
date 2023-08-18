#' Predicts the outcome using random intercept and slope mixed effects 3Trees models.
#'
#' @param mod3T The output of fit3trees.
#' @param newdata New data to which to apply the prediction.
#' @param re Logical indicating whether to include random effects in the prediction. Default is TRUE.
#' @return Predicted outcomes based on the fitted model.
#'
#' @import lme4
#' @import treeClust
#' @importFrom treeClust rpart.predict.leaves
#' @export

predictT3 <- function(mod3T, newdata, re=TRUE){
  require(lme4)
  requireNamespace("treeClust", quietly = TRUE)

  n <- nrow(newdata)
   
  ## recover the trees
  Tree1 <- nomi1 <- Tree2 <- nomi2 <- Tree3 <- nomi3 <- NULL
  if (!is.null(mod3T$Tree1)) {
    Tree1 <- treeClust::rpart.predict.leaves(mod3T$Tree1, newdata, type = "where")
    nomi1 <- .mylist.rules.rpart(mod3T$Tree1)
  }
  if (!is.null(mod3T$Tree2)) {
    Tree2 <- treeClust::rpart.predict.leaves(mod3T$Tree2, newdata, type = "where")
    nomi2 <- .mylist.rules.rpart(mod3T$Tree2)
  }
  if (!is.null(mod3T$Tree3)) {
    Tree3 <- treeClust::rpart.predict.leaves(mod3T$Tree3, newdata, type = "where")
    nomi3 <- .mylist.rules.rpart(mod3T$Tree3)
  }
  
  mydata2 <- as.data.frame(newdata)
  if(!is.null(Tree1))mydata2$Tree1 <- Tree1
  if(!is.null(Tree2))mydata2$Tree2 <- Tree2
  if(!is.null(Tree3))mydata2$Tree3 <- Tree3
  if (sum(Tree1)>n){mydata2$Tree1 <- factor(mydata2$Tree1, labels=nomi1)}
  if (sum(Tree2)>n){mydata2$Tree2 <- factor(mydata2$Tree2, labels=nomi2)}
  if (sum(Tree3)>n){mydata2$Tree3 <- factor(mydata2$Tree3, labels=nomi3)}
  
  # Predict using the 3Trees model
  reform <- if (re) NULL else NA
  Yhat <- predict(mod3T$mod.final, newdata = mydata2, re.form = reform, random.only = FALSE, type = "response")
  
  return(Yhat)
}

#' Calculate MSE, RMSE, and MAE for 3Trees models.
#'
#' @param mod3T The output of fit3trees.
#' @param newdata New data to evaluate the predictions.
#' @param re Logical indicating whether to include random effects in the prediction. Default is TRUE.
#' @param save Logical indicating whether to return the metrics as a list. Default is FALSE.
#' @return If save is TRUE, a list of MSE, RMSE, and MAE. Otherwise, the metrics are printed.
#'
#' @importFrom treeClust rpart.predict.leaves
#' @export

MSE3T <- function(mod3T, newdata, re = TRUE, save = FALSE) {
  
  if (!inherits(mod3T, "list") || !all(c("mod.final", "Tree1", "Tree2", "Tree3") %in% names(mod3T))) {
    stop("Invalid input: mod3T must be the output of fit3trees.")
  }
  
  if (!inherits(newdata, "data.frame")) {
    stop("Invalid input: newdata must be a data frame.")
  }
  
  if (!is.logical(re)) {
    stop("Invalid input: re must be a logical value.")
  }
  
  if (!is.logical(save)) {
    stop("Invalid input: save must be a logical value.")
  }
  
  if (!("Y" %in% names(newdata))) {
    stop("Invalid input: newdata must contain the outcome variable 'Y'.")
  }
  
  if (!is.numeric(newdata$Y)) {
    stop("Invalid input: 'Y' in newdata must be a numeric vector.")
  }
  
  yhat <- predictT3(mod3T, newdata, re = re)
  mse <- mean((newdata$Y - yhat)^2)
  mae <- mean(abs(newdata$Y - yhat))
  if (!save) {
    cat("MSE =", mse, "\n")
    cat("RMSE =", sqrt(mse), "\n")
    cat("MAE =", mae, "\n")
  } else {
    return(list(MSE = mse, RMSE = sqrt(mse), MAE = mae))
  }
}
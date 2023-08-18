#' Fit random intercept and slope mixed effects 3Trees models
#'
#' It can use a predictive or a post selection inference approach
#'
#' @import lme4
#' @import rpart
#' @param Y A vector of the outcome variable (continuous).
#' @param X A data frame or matrix containing the all the predictors, with labelled columns.
#' @param gr A vector or matrix representing the grouping variable(s). If selective.inference=TRUE, gr must be a numeric vector.
#' @param covLin A vector of predictor names to be included in the linear part. Default is NULL.
#' 
#' @param covT1 A vector of predictor names to be included in the first tree. Default is NULL.
#' @param covT2 A vector of predictor names to be included in the second tree. Default is NULL.
#' @param covT3 A vector of predictor names to be included in the third tree. Default is NULL.
#' @param re_form The formula for the random effect, as in \code{\link{lme4::lmer}}. Default is "(1|gr)".
#' @param selective.inference TRUE if post-selection inference is required; default is FALSE.
#' @param cp A vector of length 3 containing the rpart complexity parameters cp for each tree.
#' @param maxd A vector of length 3 containing the rpart maxdepth parameter for each tree.
#' @param nmin The minbucket setting for \code{\link{rpart::rpart}}. Assumed equal for all the trees.
#' @param prec The tolerance for the iterative tree search.
#' @param niter The maximum number of iterations.
#' @param seed The seed to set for selective inference data splitting.
#' @param extended When TRUE, it provides the final model matrix as output.
#'
#' @return A fitted mixed effects model object, containing the following elements:
#' \item{mse.final}{The mean squared error of the final model.}
#' \item{mse.random}{The mean squared error of the model with random effects.}
#' \item{mod.final}{The fitted mixed effects model object.}
#' \item{Tree1}{The \code{\link{rpart::rpart}} tree object for Tree1.}
#' \item{Tree2}{The \code{\link{rpart::rpart}} tree object for Tree2.}
#' \item{Tree3}{The \code{\link{rpart::rpart}} tree object for Tree3.}
#' \item{modmat}{The model matrix with tree predictions (if extended is TRUE).}
#'
#' @examples
#' 
#' data(mydat)
#' X<-mydat[,1:5]
#' mod3T <- fit3Trees(Y=mydat$Y, X=X, gr=mydat$gr, covLin=names(X), 
#'                  covT1=names(X[,1:3]), covT2=names(X[,4:5]), 
#'                  covT3=colnames(X),
#'                  niter = 50, re_form = "(1|gr)",
#'                  selective.inference = FALSE)
#' @seealso
#' \code{\link{lme4::lmer}}, \code{\link{rpart::rpart}}
#' 
#' @export

fit3Trees <- function(Y, X, gr,  covLin=NULL, covT1=NULL, covT2=NULL, covT3=NULL,
                                   niter = 50, re_form = "1|gr",
                                   selective.inference = FALSE,
                                   cp = c(0.0001, 0.0001, 0.0001),
                                   maxd = c(3,3,3), nmin=10,
                                   prec=1e-4, extended= FALSE, seed = 111){
  # Input validation
  if (NCOL(X) == 0) {
    stop("Invalid input: X must be a non-empty matrix.")
  }
  
  if (!is.numeric(Y) || length(Y) != NROW(X)) {
    stop("Invalid input: Y must be a numeric vector  with the same length as the same length as the number of rows in X.")
  }
  
  if (NROW(gr) != NROW(X)) {
    stop("Invalid input: gr must have the same number of rows as X.")
  }
  
  if (!is.null(covLin) & any(!(covLin %in% names(X)))) {
    stop("Invalid input: names of the variables covLin must be X.")
  }
  
  if (!is.null(covT1) & any(!(covT1 %in% names(X)))) {
    stop("Invalid input: names of the variables covT1 must be in X.")
  }
  
  if (!is.null(covT2) & any(!(covT2 %in% names(X)))) {
    stop("Invalid input: names of the variables covT2 must be in X.")
  }
  
  if (!is.null(covT3) & any(!(covT3 %in% names(X)))) {
    stop("Invalid input: names of the variables covT3 must be in X.")
  }
  
  if (niter <= 0) {
    stop("Invalid input: niter must be a positive integer value.")
  }
  
  if (length(cp)  != 3 || !all(cp >= 0 & cp <= 1) ) {
    stop("Invalid input: cp must be a vector of length three, of values between 0 and 1")
  }
  
  if (length(maxd)  != 3 ||   !all(maxd>0)  ) {
    stop("Invalid input: maxd must be a vector of length three, of values greater than 0")
  }
  
  
  ## call estimating function  
if(selective.inference){
  myris <- .PS3Tree(Y=Y, X=X, gr=gr,  covLin=covLin, covT1=covT1, covT2=covT2, covT3=covT3,
                                niter = niter, re_form = re_form,
                                cp = cp,
                                maxd = maxd, nmin=nmin,
                                prec=prec, extended= extended, seed = seed)
} else {
  myris <- .ML3Trees(Y=Y, X=X, gr=gr,  covLin=covLin, covT1=covT1, covT2=covT2, covT3=covT3,
                    niter = niter, re_form = re_form,
                    cp = cp,
                    maxd = maxd, nmin=nmin,
                    prec=prec, extended= extended)
}
  return(myris)
  }

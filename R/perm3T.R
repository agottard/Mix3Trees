#' Evaluate random intercept mixed effects 3Trees models when permuting the order of the three trees.
#'
#' It splits the data into train (trainprop  * n) and test (the rest).
#'
#' @import lme4
#' @import rpart
#' @param Y A vector of the outcome variable (continuous).
#' @param X A data frame or matrix containing the all the predictors, with labelled columns.
#' @param gr A vector or matrix representing the grouping variable(s). If selective.inference=TRUE, gr must be a numeric vector.
#' @param covLin A vector of predictor names to be included in the linear part. Default is NULL.
#' @param covT1 A vector of predictor names to be included in the first tree. Default is NULL.
#' @param covT2 A vector of predictor names to be included in the second tree. Default is NULL.
#' @param covT3 A vector of predictor names to be included in the third tree. Default is NULL.
#' @param re_form The formula for the random effect, as in lme4. Default is "(1|gr)".
#' @param cp The rpart complexity parameters cp, assumed equal all trees.
#' @param maxd The rpart maxdepth parameter, assumed equal all trees.
#' @param nmin The minbucket setting for rpart, assumed equal all trees.
#' @param prec The tolerance for the iterative tree search.
#' @param niter The maximum number of iterations.
#' @param seed The seed to set for train-and-test data splitting.
#' @param trainprop Proportion of units in the train set.
#' @return A matrix containing the MSE values for different combinations of permuted trees.
#'
#' @examples
#' data(mydat)
#' X<-mydat[,1:5]
#' mse_matrix <- perm3T(Y=mydat$Y, X, gr=mydat$gr, covLin = NULL, covT1 = names(X[,1:3]), 
#'                      covT2 = names(X[,4:5]), covT3 = colnames(X),
#'                      niter = 50, re_form = "(1|gr)", cp = 0.0001, maxd = 3, nmin = 10,
#'                      prec = 1e-4, seed = 111, trainprop = 2/3)
#'  mse_matrix                    
#'
#' @importFrom combinat permn
#' @importFrom treeClust rpart.predict.leaves
#'
#' @export

perm3T <- function(Y, X, gr,  covLin=NULL, covT1=NULL, covT2=NULL, covT3=NULL, 
                   niter = 50, re_form = "(1|gr)",
                   cp = 0.0001, maxd = 3, nmin = 10,
                   prec = 1e-4, seed = 111,
                   trainprop=2/3){
  
  if (!requireNamespace("combinat", quietly = TRUE)){
    install.packages("combinat")}
 
   # Input validation
  if (trainprop < 0 || trainprop > 1) {
    stop("trainprop must be between 0 and 1")
  }

# Define train and test data
  ### Split sample by units in each group
  mygr <- unique(gr)
  J <- length(mygr)
  n <- length(Y)
  id <- 1:n
  set.seed(seed)
  idtr <- unlist(lapply(1:J, function(j){
   idj <-  id[gr==mygr[j]]
   nj <- ceiling(length(idj)*trainprop)
  return(sample(idj, nj, replace=FALSE))
   }))
  
  idte <- setdiff(id, idtr)
  
  ### Prepare data
  Ytr <- Y[idtr]; Yte <- Y[idte]
  Xtr <- X[idtr,]; Xte <- X[idte,]
  grtr <- gr[idtr]; grte <- gr[idte]
  
  
# Define permutations  and fit the models
  all_perm <- do.call(rbind, combinat::permn(c("covT1", "covT2", "covT3")))

  Mods <- lapply(1:NROW(all_perm), function(i){
    fit3Trees(Y=Ytr, X=Xtr, gr=grtr,  covLin=covLin, 
              covT1=get(all_perm[i,1]), covT2=get(all_perm[i,2]), 
              covT3=get(all_perm[i,3]),niter = niter, re_form = re_form,
              selective.inference = FALSE,
              cp = rep(cp,3), maxd = rep(maxd, 3), nmin = nmin,
              prec = 1e-4, extended = FALSE
              )}) 
# Models comparison

pnames<- apply(all_perm, 1, paste, collapse = "-")
mydat <- data.frame(Y, X, gr)

test <- data.frame(Y, X, gr)[idte,]


ris <- lapply(1:NROW(all_perm), function(i){
           c(Mods[[i]]$mse.random, Mods[[i]]$mse.final, 
           MSE3T(Mods[[i]], newdata = test, re=TRUE, save=TRUE)[[1]], 
           MSE3T(Mods[[i]], newdata = test, re=FALSE, save=TRUE)[[1]])
   })

MSEmat <- matrix(unlist(ris), nrow = length(pnames), byrow = TRUE)
rownames(MSEmat) <- pnames
colnames(MSEmat) <- c("MSE.tr_re", "MSE.tr_no_re", "MSE.te_re", "MSE.te_no_re")
return(MSEmat)
}


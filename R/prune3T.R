#' Prune a random intercept and slope mixed effects 3Trees model
#'
#' @param mod3T The output of fit3trees estimated setting extended=TRUE
#' @param cutT1 A list of final nodes from Tree 1 to merge for pruning
#' @param cutT2 A list of final nodes from Tree 2 to merge for pruning
#' @param cutT3 A list of final nodes from Tree 3 to merge for pruning
#' @param cutLin A list of variables from covLin to omit
#' @param re_form The formula for the random effect, as in lme4, to be used to fit the pruned model. Default is "(1|gr)".
#' 
#' @return A fitted mixed effects model object.
#' 
#' @examples
#' data(mydat)
#' X<-mydat[,1:5]
#' mod3T <- fit3Trees(Y=mydat$Y, X=X, gr=mydat$gr, covLin=names(X), 
#'                  covT1=names(X[,1:3]), covT2=names(X[,4:5]), 
#'                  covT3=colnames(X),
#'                  niter = 50, re_form = "(1|gr)",
#'                  selective.inference = FALSE, extended=TRUE)
#'                  
#' mod3T.pruned <- prune3T(mod3T, cutT1=list(c(1,2)), cutT2=NULL, cutT3=NULL, cutLin=NULL)
#' summary(mod3T.pruned)
#' @import rpart lme4 rockchalk stringr PTXQC
#' @export

prune3T <- function(mod3T, cutT1=NULL, cutT2=NULL, cutT3=NULL, cutLin=NULL, re_form = "1|gr"){
  require(lme4)
  requireNamespace("rpart", quietly = TRUE)
  requireNamespace("rockchalk", quietly = TRUE)
  requireNamespace("stringr", quietly = TRUE)
  requireNamespace("PTXQC", quietly = TRUE)
  
  if (is.null(mod3T$modmat)){stop("Invalid input: mod3T has to be fitted with extended=TRUE to use prune3T.")}
  
  if (!is.null(cutT1) & !is.list(cutT1)) {
    stop("Invalid input: cutT1 has to be or a list or NULL.")
  }
  
  if (!is.null(cutT2) & !is.list(cutT2)) {
    stop("Invalid input: cutT2 has to be or a list or NULL.")
  }
  
  if (!is.null(cutT3) & !is.list(cutT3)) {
    stop("Invalid input: cutT3 has to be or a list or NULL.")
  }
  
  if (!is.null(cutLin) & !is.list(cutLin)) {
    stop("Invalid input: cutLin has to be or a list or NULL.")
  }

  
##### Nodes merging    
  mymod2 <- mod3T
  if(!is.null(cutT1)){
    for (i in 1: length(cutT1)){
      T11<-levels(mod3T$modmat$Tree1)[grepl( paste0(":",cutT1[[i]][1]) , levels(mod3T$modmat$Tree1))]
      T12<- levels(mod3T$modmat$Tree1)[grepl( paste0(":",cutT1[[i]][2]) , levels(mod3T$modmat$Tree1))]
      namenew<-PTXQC::LCSn(c(T11, T12), min_LCS_length = 0)
      namenew<- paste0(":", cutT1[[i]][1],"+", cutT1[[i]][2], namenew)
      namenew<-stringr::word(namenew,1,-4)
      mymod2$modmat$Tree1<-rockchalk::combineLevels(mod3T$modmat$Tree1,levs = c(T11, T12), newLabel = c(namenew) )
    }}

  if(!is.null(cutT2)){
    for (i in 1: length(cutT2)){
      T21<-levels(mod3T$modmat$Tree2)[grepl( paste0(":",cutT2[[i]][1]) , levels(mod3T$modmat$Tree2))]
      T22<- levels(mod3T$modmat$Tree2)[grepl( paste0(":",cutT2[[i]][2]) , levels(mod3T$modmat$Tree2))]
      namenew<-PTXQC::LCSn(c(T21, T22), min_LCS_length = 0)
      namenew<- paste0(":", cutT2[[i]][1],"+", cutT2[[i]][2], namenew)
      namenew<-stringr::word(namenew,1,-4)
      mymod2$modmat$Tree2<-rockchalk::combineLevels(mod3T$modmat$Tree2,levs = c(T21, T22), newLabel = c(namenew) )
    }}

  if(!is.null(cutT3)){
    for (i in 1: length(cutT3)){
      T31<-levels(mod3T$modmat$Tree3)[grepl( paste0(":",cutT3[[i]][1]) , levels(mod3T$modmat$Tree3))]
      T32<- levels(mod3T$modmat$Tree3)[grepl( paste0(":",cutT3[[i]][2]) , levels(mod3T$modmat$Tree3))]
      namenew<-PTXQC::LCSn(c(T31, T32), min_LCS_length = 0)
      namenew<- paste0(":", cutT3[[i]][1],"+", cutT3[[i]][2], namenew)
      namenew<-stringr::word(namenew,1,-4)
      mymod2$modmat$Tree3<-rockchalk::combineLevels(mod3T$modmat$Tree3,levs = c(T31, T32), newLabel = c(namenew) )
    }}

  #### Remove variables from covLin if specified
  newmodmat <- names(mod3T$modmat)
  if(!is.null(cutLin)) newmodmat <- setdiff(newmodmat, cutLin)
  
  ### Fit reduced model
  cova <- paste0(newmodmat[-c(which(newmodmat=="Y"),which(newmodmat=="gr"))], collapse="+")
  myform3 <- as.formula(paste0("Y ~ ", cova, " + (", re_form, ")"))
  mod.rid <- lmer(myform3, REML=FALSE, data=mymod2$modmat)
  
  return(mod.rid)
}

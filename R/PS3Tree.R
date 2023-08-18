.PS3Tree <- function(Y, X, gr,  covLin=NULL, covT1=NULL, covT2=NULL, covT3=NULL,
                     niter = 50, re_form = "(1|gr)",
                     cp = c(0.0001, 0.0001, 0.0001),
                     maxd = c(3,3,3), nmin=10,
                     prec=1e-4, extended= FALSE, seed = 111){
 
   # Input validation
  
   if (!is.vector(gr)) {
    stop("Invalid input: At the moment, the Selective inference procedure is avalable only for 
         two-level mixed effect model. gr has to be a vector")
  }
  
  if (length(gr) != nrow(X)) {
    stop("Invalid input: gr must have the same length as the number of rows in X.")
  }
  
  if (!grepl("^\\(.*\\)$", re_form)) {
    re_form <- paste0("(", re_form, ")")
  }
  
  
### Setting data splitting
  mydata <- as.data.frame(cbind(Y,X, gr))
  mydata <- mydata[order(mydata$gr),]
  p <- ncol(X)
  n <- nrow(X)
  
  ### Split sample by group
  mygr <- unique(gr)
  J <- length(mygr)
  set.seed(seed)
 
  id.J <- sort(sample(mygr, J2, replace=FALSE))
  id.J2 <- setdiff(mygr, id.J)
  
  ### Prepare data
  mydata1 <- mydata[mydata$gr %in%id.J,]   ## selection set
  mydata2 <- mydata[mydata$gr %in%id.J2,]  ## estimation set
  n1 <- dim(mydata1)[1]
  n2 <- dim(mydata1)[2]
  
  sample.id<-function(nj, J, prop.nj=1/3, seed=111){
    set.seed(seed)
    nt <- floor(nj*prop.nj)
    idplace <- sort(sample(nj, nt, replace = FALSE))
    idtest<-idplace+ rep(0:(J-1)*nj, each=nt)
    return(idtest)
  }
  

  if(!is.null(covT1)) covT1 <- paste0(covT1, collapse="+")
  if(!is.null(covT2)) covT2 <- paste0(covT2, collapse="+")
  if(!is.null(covT3)) covT3 <- paste0(covT3, collapse="+")
  if(!is.null(covLin)) covLin  <- paste0(covLin , collapse="+")
  if(is.null(covLin)) covLin  <- "1"

  myform1 <- as.formula(paste0("Y_pres ~ ", covLin, " + ", re_form))
  if(!is.null(covT1)) myform.t.1 <- as.formula(paste0("Y_pres~", covT1))
  if(!is.null(covT2)) myform.t.2 <- as.formula(paste0("Y_pres~", covT2))
  if(!is.null(covT3)) myform.t.3 <- as.formula(paste0("Y_pres~", covT3))

  #### Selection step
  library(lme4)
  cp1=cp[1]; cp2=cp[2]; cp3=cp[3]
  maxd1=maxd[1]; maxd2=maxd[2]; maxd3=maxd[3]
  
  #initialization 
  YhatT1 <- YhatT2 <- YhatT3 <- mean(mydata1$Y)/3
  YhatL <- predict(lm(as.formula(paste0("Y ~ ", covLin)), data=mydata1))
  
  msetrain<- mse.best<- 10000000
  t <-0 
  d <-1
  
  options(warn = -1)
  
  #### Start iterative procedure
  while (d != 0){
    msetrain_old <- msetrain
    ##### Compute residuals
    Y_residuals <- mydata1$Y - YhatL - YhatT1 - YhatT2 - YhatT3
    
    ### Linear part
    Y_pres <- Y_residuals + YhatL
    mod1 <- suppressWarnings(lmer(myform1, REML=FALSE, data=mydata1, verbose = 0))
    YhatL <- predict(mod1, re.form = NULL, random.only=FALSE, type="response")

    ### fit tree T1
    Y_pres <- Y_pres - YhatL + YhatT1
    YhatT1 <- rep(0,n1); treefit.1<- NULL
    if(!is.null(covT1))treefit.1<-rpart::rpart(formula = myform.t.1, maxdepth = maxd1, cp=cp1, data=mydata1, maxsurrogate = 0, minbucket=nmin)
    if(!is.null(covT1))YhatT1 <- predict(treefit.1)
    

    ### fit tree T2
    Y_pres <- Y_pres -  YhatT1 + YhatT2
    YhatT2 <- rep(0,n1); treefit.2=NULL
    if(!is.null(covT2))treefit.2<-rpart::rpart(formula = myform.t.2, maxdepth = maxd2, cp=cp2, data=mydata1, maxsurrogate = 0, minbucket=nmin)
    if(!is.null(covT2))YhatT2 <- predict(treefit.2)

    ### fit tree T3
    Y_pres <- Y_pres -  YhatT2 + YhatT3
    YhatT3 <- rep(0,n1); treefit.3=NULL
    if(!is.null(covT3))treefit.3<-rpart::rpart(formula = myform.t.3, maxdepth = maxd3, cp=cp3, data=mydata1, maxsurrogate = 0, minbucket=nmin)
    if(!is.null(covT3))YhatT3 <- predict(treefit.3)

    pred.final <- YhatL + YhatT1 + YhatT2 + YhatT3
    msetrain<- mean((mydata1$Y-pred.final)^2)
    
    t=t+1
    cat(paste0(rep('=', t), collapse = ''))
    d <- (abs(msetrain_old - msetrain)> prec)*(t < niter)
    #update
    if (msetrain<mse.best) {mse.best<- msetrain; treefit.best.1 <- treefit.1; treefit.best.2 <- treefit.2; treefit.best.3 <- treefit.3}
  }#### End iterative procedure
  
  options(warn = 0)

  #########################
  # Compute model on selection set
  Tree1 <- Tree2 <- Tree3<- NULL
  nomi.1 <- nomi.2 <- nomi.2 <- NULL
  if(!is.null(treefit.best.1))nomi.1 <- .mylist.rules.rpart(treefit.best.1)
  if(!is.null(treefit.best.1))Tree1 <- treeClust::rpart.predict.leaves(treefit.best.1, mydata1, type = "where")
  if(!is.null(treefit.best.2))nomi.2 <- .mylist.rules.rpart(treefit.best.2)
  if(!is.null(treefit.best.2))Tree2 <- treeClust::rpart.predict.leaves(treefit.best.2, mydata1, type = "where")
  if(!is.null(treefit.best.3))nomi.3 <- .mylist.rules.rpart(treefit.best.3)
  if(!is.null(treefit.best.3))Tree3 <- treeClust::rpart.predict.leaves(treefit.best.3, mydata1, type = "where")
  
  mydata1s <- mydata1
  if(!is.null(Tree1))mydata1s <- cbind.data.frame(mydata1s,Tree1)
  if(!is.null(Tree2))mydata1s <- cbind.data.frame(mydata1s,Tree2)
  if(!is.null(Tree3))mydata1s <- cbind.data.frame(mydata1s,Tree3)
  if (sum(Tree1)>n){mydata1s$Tree1 <- factor(mydata1s$Tree1, labels=nomi.1)}
  if (sum(Tree2)>n){mydata1s$Tree2 <- factor(mydata1s$Tree2, labels=nomi.2)}
  if (sum(Tree3)>n){mydata1s$Tree3 <- factor(mydata1s$Tree3, labels=nomi.3)}

  myform3 <- paste0("Y ~ ", covLin)
  if(covLin == "1")myform3 <- "Y ~ "
  if(!is.null(Tree1))myform3 <- paste0(myform3, "+Tree1")
  if(!is.null(Tree2))myform3 <- paste0(myform3, "+Tree2")
  if(!is.null(Tree3))myform3 <- paste0(myform3, "+Tree3")
  myform3 <- as.formula(paste0(myform3, "+", re_form))
  
  mod.sel <- lmer(myform3, REML=FALSE, data=mydata1s)

  #########################
  #  set data for estimation step

  if(!is.null(Tree1))Tree1 <- treeClust::rpart.predict.leaves(treefit.best.1, mydata2, type = "where")
  if(!is.null(Tree1))nomi.1 <- .mylist.rules.rpart(treefit.best.1)
  if(!is.null(Tree2))Tree2 <- treeClust::rpart.predict.leaves(treefit.best.2, mydata2, type = "where")
  if(!is.null(Tree2))nomi.2 <- .mylist.rules.rpart(treefit.best.2)
  if(!is.null(Tree3))Tree3 <- treeClust::rpart.predict.leaves(treefit.best.3, mydata2, type = "where")
  if(!is.null(Tree3))nomi.3 <- .mylist.rules.rpart(treefit.best.3)

  mydata3 <- mydata2
  if(!is.null(Tree1))mydata3 <- cbind.data.frame(mydata3,Tree1)
  if(!is.null(Tree2))mydata3 <- cbind.data.frame(mydata3,Tree2)
  if(!is.null(Tree3))mydata3 <- cbind.data.frame(mydata3,Tree3)
  
  if (sum(Tree1)>n){mydata3$Tree1 <- factor(mydata3$Tree1, labels=nomi.1)}
  if (sum(Tree2)>n){mydata3$Tree2 <- factor(mydata3$Tree2, labels=nomi.2)}
  if (sum(Tree3)>n){mydata3$Tree3 <- factor(mydata3$Tree3, labels=nomi.3)}

  myform3 <- paste0("Y ~ ", covLin)
  if(covLin == "1")myform3 <- "Y ~ "
  if(!is.null(Tree1))myform3 <- paste0(myform3, "+Tree1")
  if(!is.null(Tree2))myform3 <- paste0(myform3, "+Tree2")
  if(!is.null(Tree3))myform3 <- paste0(myform3, "+Tree3")
  myform3 <- as.formula(paste0(myform3, "+", re_form))

  mod.final <- lmer(myform3, REML=FALSE, data=mydata3)
  
  mse.final<- mean((mydata3$Y - predict(mod.final, re.form = ~0, random.only=FALSE, type="response"))^2)
  mse.random  <- mean((mydata3$Y - predict(mod.final, re.form = NULL, random.only=FALSE, type="response"))^2)
  ifelse(extended==TRUE,
         return(list(tree=list(Tree1=treefit.best.1, Tree2=treefit.best.2, Tree3=treefit.best.3), mse.final=mse.final, mod.selected = mod.sel, mse.random=mse.random, mod.final=mod.final, modmat = mydata3, modmat.sel = mydata1s, id_sets=list(sel=id.J, est=id.J2))),
         return(list(tree=list(Tree1=treefit.best.1, Tree2=treefit.best.2, Tree3=treefit.best.3), mse.final=mse.final,mse.random=mse.random, mod.final=mod.final, id_sets=list(sel=id.J, est=id.J2))))
}

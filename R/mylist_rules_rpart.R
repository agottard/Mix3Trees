# This is a modified version of the function originally 
# proposed by Bharat Warule and can be found at 
# https://github.com/Bwarule/Data-Mining/blob/master/list.rules.rpart.R. 
# I would like to acknowledge and give credit to the original author.

.mylist.rules.rpart <- function(model){
  if (!inherits(model, "rpart")) stop("Not a legitimate rpart tree")
  frm     <- model$frame
  names   <- row.names(frm)
  nodenames<- c()
  for (i in 1:nrow(frm)){
    if (frm[i,1] == "<leaf>"){
      a<- sprintf(":%s", names[i])
      pth <- rpart::path.rpart(model, nodes=as.numeric(names[i]), print.it=FALSE, pretty = 1)
      b<- capture.output(cat(sprintf("  %s ", unlist(pth)[-1]), sep="&"))
      nodenames[i] <- paste0(a,b)
    }
  }
  return(as.vector(na.omit(nodenames)))
}

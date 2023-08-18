#' Simulated data for using 3Trees models.
#'
#' This is an example dataset that contains simulated data.
#' It is used for demonstration purposes in the package examples.
#'
#' @format A data frame with 5000 rows and 7 columns.
#'   
#'   Columns 1-3: Unit level covariates X1,X2,X3, with X2 and X3 irrelevant.
#'   Columns 4-5: Group level covariates Z1,Z2, with Z1 having a linear effect.
#'   Column 6: Y response variable
#'   Column 7: gr group indicator
#'   
#' @source This dataset was created for the package examples.
#'
#' @examples
#' data(mydat)
#' summary(mydat)
#'
#' @seealso
#' \code{\link{fit3Trees}}, \code{\link{perm3T}}, \code{\link{prune.3T}},
#' \code{\link{MSE3T}}
#'
#' @name mydat
#' @docType data
#' @usage data(mydat)



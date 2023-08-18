
This `R` package provides functions to estimate *3Tree* models, as described in the paper

Gottard, A., Vannucci, G., Grilli, L., & Rampichini, C. (2023). Mixed-effect models with trees. Advances in Data Analysis and Classification, 17(2), 431-461.

## Installation


To install the package, run 
```
install.packages("devtools")
devtools::install_github("agottard/Mix3Trees")
library(Mix3Trees)
```

## Example

This is a basic example which shows you how to use the main function:

```{r example}
library(Mix3Trees)

### I will use a simulated dataset as illustration
data(mydat)

names(mydat)

X<-mydat[,1:5]
mod1 <- fit3Trees(Y=mydat$Y, X=X, gr=mydat$gr,
                  covLin=names(X), covT1=names(X[,1:3]),
                  covT2=names(X[,4:5]), covT3=colnames(X),
                  niter = 50, re_form = "(1|gr)",
                  selective.inference = FALSE)

summary(mod1$mod.final)

```






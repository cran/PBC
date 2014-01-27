###########
## tests ##
###########

require(PBC)

## Set the underlying graphical structure for the PBC model
g <- graph.formula(X1-X2, X2-X3, X3-X4, X4-X5, simplify = FALSE)

## Create the PBC object with Gumbel linking family
myPBC <- pbc(g, model="gumbel") # or:
myPBC <- pbcGumbel(g)

## Plot PBC graph
pbcPlot(myPBC)

## Generate n observations from the model
theta <- 1:4
n <- 100
data <- rPBC(n, theta, myPBC)

## Estimate the parameter vector
init <- rep(5, 4) # the 'par' argument of \code{\link{optim}}.
# it's better if you can provide an estimate based on pairwise likelihood to 
# increase the chances to get a good minimizer.

## Use \code{\link{pbcOptim}}
fitPBC <- pbcOptim(init, data, myPBC, method='BFGS')
fitPBC # estimate

## You may use \code{\link{optim}} instead
fn <- function(theta)-sum(log((dPBC(data, theta, myPBC)))) # -log likelihood
gr.temp <- function(u, theta)mpAlgo(myPBC, u, theta)@gradient # gradient of likelihood
gr <- function(theta){ # gradient of -log likelihood
  ap <- t(apply(data, 1, gr.temp, theta=theta))
  ap2 <- dPBC(data, theta, myPBC)
  apply(-ap*ap2, 2, sum)
}
fitPBC2 <- optim(par=init, fn=fn, gr=gr, method="BFGS")
fitPBC2

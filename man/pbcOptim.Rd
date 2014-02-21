\name{pbcOptim}
\alias{pbcOptim}
\title{
Maximum likelihood estimation for a PBC copula
}
\description{
This function performs maximum-likelihood inference in a PBC copula. 
}
\usage{
pbcOptim(par, data, pbcObj, method, lower = -Inf, upper = Inf, ...)
}
\arguments{
  \item{par}{the argument \code{par} in \code{optim}.}
  \item{data}{data matrix. }
  \item{pbcObj}{an object of class PBC. }
  \item{method}{method for the \code{optim} function. Two methods are available: Broyden-Fletcher
  -Goldfarb-Shanno (BFGS) and limited memory BFGS with bounds (L-BFGS-B). }
  \item{lower,upper}{bounds on the variables for the "L-BFGS-B" method. }
  \item{\dots}{ currently nothing. }
}

\value{
A numeric vector giving the estimate.
}

\details{
The R routine \code{\link{optim}} is used to minimize -log likelihood. Compared to brute use of \code{optim}, \code{pbcOptim} saves one call to \code{mpAlgo}. See \code{\link{optim}}.
}

\source{
The code for pbcOptim is based on that of \code{\link{optim}}. In particular, pbcOptim calls /src/lbfgsb.cpp which is a slight adaptation of lbfgsb.c and optim.c part of the R software (http://www.r-project.org/). See \code{\link{optim}} and ?optim.

The code for BFGS method is based on Pascal code in J.C. Nash, 'Compact Numerical Methods for Computers', 2nd edition, converted by p2c then re-crafted by B.D. Ripley. It is modified by V.T. Pham for pbcOptim to increase speed of this function.

The code for L-BFGS-B method is based on Fortran code by Zhu, Byrd, Lu-Chen and Nocedal obtained from Netlib, then modified by V.T. Pham.
}

\seealso{
\code{\link{optim}} 
}

\examples{
## Example with normal linking family
g <- graph.formula(X1-X4, X4-X2, X2-X3, X4-X5, simplify = FALSE)
pbcNormal <- pbcNormal(g)
theta <- runif(4)
pbcDataNormal <- rPBC(5, theta, pbcNormal)

## estimation 
# L-BFGS-B method with Lower and upper bound 
init <- rep(.5,4) # it's better if you can provide an estimate based 
# on pairwise likelihood to increase the chances to get a good minimizer.

opt = pbcOptim(init, pbcDataNormal, pbcNormal, 
               method="L-BFGS-B", lower=rep(0,4), upper=rep(0.99,4))
# BFGS method
opt2 = pbcOptim(init, pbcDataNormal, pbcNormal, method="BFGS")
}

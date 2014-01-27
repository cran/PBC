\name{PBC}
\alias{PBC}
\alias{rPBC}
\alias{pPBC}
\alias{dPBC}
\title{
Random number generation, distibution function and density for the PBC model
}
\description{
Random generation (rPBC), distribution function (pPBC) and density (dPBC) for the PBC model.
}
\usage{
rPBC(n, theta, pbcObj, \dots)
pPBC(u, theta, pbcObj, \dots)
dPBC(u, theta, pbcObj, \dots)
}
\arguments{
  \item{pbcObj}{an object of class PBC.}
  \item{n}{number of observations to be generated.}
  \item{theta}{parameter vector.}
  \item{u}{a vector or a matrix at which the function needs to be evaluated.}
  \item{\dots}{currently nothing.}
}
\details{
  The density is computed via the function \code{\link{mpAlgo}}. The parameter vector has length equal to the number of variables minus one. If \code{u} is a matrix, it has dimension (n,d) where n is the number of vectors the function is evaluated at, and d is the number of variables in the model.
}
\value{
  \code{rPBC()} generates random data, \code{pPBC()} computes the distribution function and
  \code{dPBC()} computes the density.
}
%\keyword{multivariate }
\seealso{
\code{\link{mpAlgo}},\code{\link{PBC-class}}.
}
\examples{
## set a parameter vector
theta <- runif(4)
## construct the graph 
g <- graph.formula(X1-X4, X4-X2, X2-X3, X4-X5, simplify = FALSE)
## create the PBC object with linking family "AMH"
myPBC.AMH <- pbc(g, model="amh")
## alternatively:
## myPBC.AMH <- pbcAMH(g)
## Generate 5 random data vectors
r1 <- rPBC(5, theta, myPBC.AMH)
## Compute the distribution function 
p1 <- pPBC(r1, theta, myPBC.AMH)
## Compute the density 
d1 <- dPBC(r1, theta, myPBC.AMH)
}


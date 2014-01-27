\name{mpAlgo}
\alias{mpAlgo}
\title{
Message-Passing algorithm
}
\description{
The function \code{mpAlgo} implements a message-passing algorithm to compute the likelihood and its gradient (with respect to parameter vector) for the PBC model.
}
\usage{
mpAlgo(pbcObj, u, theta, output="both", ...)
}
\arguments{
  \item{pbcObj}{an object of class PBC. }
  \item{u}{point at which is evaluated the density and the gradient.  }
  \item{theta}{parameter vector. }
  \item{output}{desired output: "density" or "gradient" or "both". The default is "both". }
  \item{\dots}{currently nothing. }
}

\value{
An object of class \code{\link{PBC}} The likelihood is found in the slot \code{@density} and its gradient in the slot \code{@gradient} (see example below).
}

\details{
\code{\link{mpAlgo}} implements the algorithm of Table 1 in J. C. Huang and N. Jojic, 2010, for the special case of a PBC model.
}
\references{
J. C. Huang and N. Jojic. Maximum-likelihood learning of cumulative distribution functions on graphs. Journal of Machine Learning Research  W&CP Series, 9:342--349, 2010.
\cr
G. Mazo G, S. Girard and F. Forbes. A class of high dimensional copulas based on products of bivariate copulas. http://hal.archives-ouvertes.fr/hal-00910775.
}
%\keyword{message-passing}
\seealso{
\code{\link{PBC-class}}
}
\examples{
## Example with Gumbel linking family 
g <- graph.formula(X1-X3,X2-X3,X3-X4,X4-X5,X4-X6,X6-X7,X6-X8,simplify = FALSE)
pbcGumbel <- pbc(g, model="gumbel")
u <- runif(8)
theta <- 1/runif(7)
pbcOut <- mpAlgo(pbcGumbel, u, theta)
pbcOut@density
pbcOut@gradient
## Example with a user defined dsitribution (must contain 'x' and 'y')
f <- expression(exp(-((-log(x))^(theta)+(-log(y))^(theta))^(1/theta))) # Gumbel
pbcUser <- pbc(g, model=f)
pbcOut2 <- mpAlgo(pbcUser, u, theta)
pbcOut2@density
pbcOut2@gradient
}

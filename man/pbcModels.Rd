\name{pbcModels}
\alias{pbcModels}
\alias{pbcGumbel}
\alias{pbcFGM}
\alias{pbcFrank}
\alias{pbcNormal}
\alias{pbcAMH}
\alias{pbcJoe}
\title{
Linking copula families for the PBC model
}
\description{
Linking copula families implemented in the \bold{PBC} package.
}
\usage{
pbcGumbel(graph)
pbcFGM(graph)
pbcFrank(graph)
pbcNormal(graph)
pbcAMH(graph)
pbcJoe(graph)
}
\arguments{
  \item{graph}{the graph (of class \code{\link{igraph}}) associated to the PBC copula.}
}
\details{
  A pair \eqn{(U_i,U_j)} of the PBC model has copula 
  \deqn{C_{ij}(u,v) = u^{ 1 - 1/n_i } * v^{ 1 - 1/n_j } * D_{ij}( u^{ 1/n_i }, v^{ 1/n_j } ),}
  where \eqn{n_i} and \eqn{n_j} are the number of neighbors in the graph for the variables \eqn{U_i} and \eqn{U_j} respectively (G. Mazo G, S. Girard and F. Forbes). The copula families implemented for \eqn{D_{ij}(u,v)} are given below.
  
  \describe{
  \item{\code{pbcGumbel}:}{family of Gumbel copulas:
         \deqn{ \exp(-((-\ln(u)))^\theta + (-\ln(v))^\theta)^{1/\theta})}{
        exp(-((-ln(u))^theta + (-ln(v))^theta)^(1/theta)) }
       with
     \eqn{\theta\in[1,\infty)}{theta in [1,Inf)}. }
  \item{\code{pbcFGM}:}{family of Farlie-Gumbel-Morgenstern (FGM) copulas:
           \deqn{u * v * (1 + \theta * (1 - u) * (1 - v))}{
        u * v * (1 + theta * (1 - u) * (1 - v)) }
       with
     \eqn{\theta \in [-1,1]}{theta in [-1,1]}. }
  \item{\code{pbcFrank}:}{family of Frank copulas:
           \deqn{- \ln(1 + (\exp(- \theta * u) - 1)*(\exp(- \theta * v) - 1)/(\exp(-\theta) - 1))/\theta}{
      - ln(1 + (exp(- theta * u) - 1)*(exp(- theta * v) - 1)/(exp(-theta) - 1))/theta }
       with
     \eqn{\theta\in(0,\infty)}{theta in (0,Inf)}. }
  \item{\code{pbcNormal}:}{family of normal copulas:
           \deqn{\exp(((\theta *  q(u))^2 + (\theta * q(v))^2 - 2 * \theta
           *  q(u) * q(v)) / (2 * (- 1 + \theta^2))) / ( 1 - \theta^2)^0.5}{
           exp(((theta *  q(u))^2 + (theta * q(v))^2 - 2 * theta
           *  q(u) * q(v)) / (2 * (- 1 + theta^2))) / ( 1 - theta^2)^0.5
            }
       with
     \eqn{\theta \in [-1,1])}{theta in [-1,1])}, q is the inverse of the standard normal distribution function. }
  \item{\code{pbcAMH}:}{Family of Ali-Mikhail-Haq (AMH) copulas:
         \deqn{u * v / (1 - \theta * (1 -u) * (1 -v))}{
        u * v / (1 - theta * (1 -u) * (1 -v)) }
       with
     \eqn{\theta \in [0,1)}{theta in [0,1)}. }
  \item{\code{pbcJoe}:}{Family of Joe copulas:
         \deqn{1 - ((1 - u)^\theta + (1 - v)^\theta - (1 - u)^\theta * (1 - v)^\theta)^{1/\theta}}{
        1 - ((1 - u)^theta + (1 - v)^theta - (1 - u)^theta * (1 - v)^theta)^(1/theta) }
       with
     \eqn{\theta\in[1,\infty)}{theta in [1,Inf)}. }
  }
}
\value{
  A \code{"\link{PBC}"} object.
}
\references{
G. Mazo G, S. Girard and F. Forbes. A class of high dimensional copulas based on products of bivariate copulas. http://hal.archives-ouvertes.fr/hal-00910775.
R. B. Nelsen. An Introduction to Copulas. Springer, 1999. \cr
}

\keyword{ models }
\seealso{
\code{\link{pbc}} 
}
\examples{
## Example for the FGM family
graph <- graph.formula(X1-X2,X2-X3, simpify = FALSE)
## Create a PBC object
fgmObject <- pbcFGM(graph)
## alternatively
g <- graph.formula(X1-X2,X2-X3, simpify = FALSE)
fgmObject <- pbc(g, model="fgm")
}

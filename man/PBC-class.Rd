\name{PBC-class}
\Rdversion{1.0}
\title{Class "PBC" for the PBC model}
\docType{class}
\alias{getRoot,PBC-method}
\alias{getGraph,PBC-method}
\alias{getF,PBC-method}
\alias{getNIteration,PBC-method}
\alias{getDxf,PBC-method}
\alias{getDxdyf,PBC-method}
\alias{getGraf,PBC-method}
\alias{getGradxf,PBC-method}
\alias{getGradxdyf,PBC-method}
\alias{getBINMAT,PBC-method}
\alias{getModel,PBC-method}
\alias{setDensity,PBC-method}
\alias{setGradient,PBC-method}
\alias{pbcPlot,PBC-method}
\alias{pbcPlot}
\alias{getDxdyf}
\alias{getDxf}
\alias{getF}
\alias{getGradxdyf}
\alias{getGradxf}
\alias{getGraf}
\alias{setGradient}
\alias{setDensity}
\alias{phi.student1}
\alias{phi.student}
\alias{phi.norm}
\alias{margin}
\alias{compute}
\alias{gradxf2.student}
\alias{gradxf.norm}
\alias{gradxdyphi2.student}
\alias{gradxdyphi.norm}
\alias{getModel}
\alias{getRoot}
\alias{getNIteration}
\alias{getGraph}
\alias{dxf.student}
\alias{dxf.norm}
\alias{dxdyphi.student}
\alias{dxdyphi.norm}
\alias{getBINMAT}
\alias{draw}
\alias{getBinMat}
\alias{pbc}
\alias{igraph-class}
\alias{PBC-class}
\description{
  The class \code{"PBC"} provides a function to create PBC objects. 
}
\usage{
pbc(g, model, ...)
}
\arguments{
  \item{g}{an \code{\link{igraph}} object.}
  \item{model}{a family of copulas among "gumbel", "fgm", "frank", "normal", "amh", "joe".}
  \item{\dots}{currently nothing.}
}
\section{Slots}{
  \describe{
    \item{\code{graph}:}{a graph of class \code{\link{igraph}} describing
      the PBC where the nodes represent the variables and 
      the edges represent the linking bivariate copulas.}
    \item{\code{root}:}{the root for the message-passing algorithm (center of the graph associated to the PBC model).}
    \item{\code{nIteration}:}{number of iterations for the message-passing algorithm.}
    \item{\code{BINMAT}:}{a matrix encoding the links between
    the variables (represented by nodes) and the bivariate copulas (represented by edges) in the
    graph. An element of the matrix is set to 1 if a variable and a bivariate copula are adjacent, 0 otherwise. }
    \item{\code{model}:}{a copula family to link the variable nodes. Choices include "gumbel", "fgm", "frank", "normal", "amh", "joe" (see \code{\link{pbcModels}} for details).}
    \item{\code{density}, \code{gradient}:}{the density and gradient (with respect to the parameter vector) obtained from the message-passing algorithm.}
  }
  For more details about the linking copula families, see \code{\link{pbcModels}}.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature("PBC")}: set the slot values, construct the graph and compute the encoding matrix. }
    \item{pbcPlot}{\code{signature("PBC")}: display the graph of the PBC model. }
  }
}

\value{
An object of class \code{PBC}.
}

\seealso{
  Linking copula families are detailed in \code{\link{pbcModels}}.
}
\examples{
## PBC class information
showClass("PBC")
## Create a PBC object with linking family "Gumbel"
g <- graph.formula(X1-X3,X2-X3,X3-X4,X4-X5,simplify = FALSE)
pbcObj <- pbc(g, model="gumbel")
}
%\keyword{multivariate}
%\keyword{classes}

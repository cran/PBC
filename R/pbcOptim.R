## Find parameter value that minimize function given
## Two methods available
## Broyden–Fletcher–Goldfarb–Shanno (BFGS) and 
## Limited-memory BFGS with bounds (L-BFGS-B)
pbcOptim <- function(par, data, pbcObj, method,
                     lower = -Inf, upper = Inf, ...) {
  theta <- par # for compatibility with Trung's code
  data <- t(data)
  margins <- expression(x^(1/n))
  dxg <- expression(x^(1/n-1)/n)
  graph <- getGraph(pbcObj)
  root <- getRoot(pbcObj)
  f <- pbcObj@f
  dxf <- pbcObj@dxf
  dxdyf <- pbcObj@dxdyf
  graf <- pbcObj@graf
  gradxf <- pbcObj@gradxf
  gradxdyf <- pbcObj@gradxdyf
  rootID <- as.numeric(substr(root,2,nchar(root)))
  BINMAT <- getBINMAT(pbcObj)
  model <- getModel(pbcObj)
  type = 0
  if (model == "gumbel")
    type = 2
  if (model == "fgm")
    type = 3
  if (model == "frank")
    type = 4
  if (model == "normal")
    type = 5
  if (model == "amh")
    type = 6
  if (model == "joe")
    type = 7
  nIteration <- getNIteration(pbcObj)
  out <- 2
  ret <- 0
  if (method == "BFGS")
    ret = .Call('bfgs', PACKAGE = 'PBC', theta, BINMAT, data, rootID, f, nIteration, dxf, 
                            dxdyf, graf, gradxf, gradxdyf, type, margins, dxg, out)
  if (method == "L-BFGS-B") {
    ret = .Call('lbfgsb', PACKAGE = 'PBC', theta, lower, upper, 
                BINMAT, data, rootID, f, nIteration, dxf, 
                dxdyf, graf, gradxf, gradxdyf, type, margins, dxg, out)
  }
  return (ret)
}
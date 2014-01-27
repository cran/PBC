## Create PBC object
pbc <- function(g, model, ...) {  
  ## Delete 'simplify=FALSE'
  g <- delete.vertices(g,V(g)[length(V(g))])
  ## Get root in center of tree
  root <- V(g)[get.diameter(g)[diameter(g)/2+1]]$name
  ## Get the longest path from root
  nItertaion <- as.integer((diameter(g)+1)/2)
  ## Check whether root given is in graph
  isValidRoot <- FALSE
  for (i in 1:length(V(g)))
    if (get.vertex.attribute(g, 'name')[i] == root)
      isValidRoot <- TRUE
  ## Check whether graph given is a tree
  if ((length(E(g)) != length(V(g)) - 1) || (!is.connected(g, mode=c("weak", "strong")))) 
    ret <- "Error! Graph is not a tree!"
  ## Check whether root given is in graph
  else if (!isValidRoot)
    ret <- "Error! Invalid root!"
  else {
    ## Initialize object
    #if ((model == "Gumbel") && (!is.null(model)))
    #  distribution <- expression(0)
    distribution <- model
    if (is.character(model) && (model == "gumbel") && (!is.null(model)))
      distribution <- expression(1)
    if (is.character(model) && (model == "fgm") && (!is.null(model)))
      distribution <- expression(2)
    if (is.character(model) && (model == "frank") && (!is.null(model)))
      distribution <- expression(3)
    if (is.character(model) && (model == "normal") && (!is.null(model)))
      distribution <- expression(4)
    if (is.character(model) && (model == "amh") && (!is.null(model)))
      distribution <- expression(5)
    if (is.character(model) && (model == "joe") && (!is.null(model)))
      distribution <- expression(6)
    #if ((model == "student") && (!is.null(model)))
    #  distribution <- expression(7)
    ret <- new("PBC", root, g, distribution, nItertaion)
  }    
  return (ret);
}

## Learn model with parameters given
mpAlgo <- function(pbcObj, u, theta, output="both", ...){
  
  
  ## Get parameters from object
  x <- u
  graph <- getGraph(pbcObj)
  root <- getRoot(pbcObj)
  f <- getF(pbcObj)
  dxf <- getDxf(pbcObj)
  dxdyf <- getDxdyf(pbcObj)
  graf <- getGraf(pbcObj)
  gradxf <- getGradxf(pbcObj)
  gradxdyf <- getGradxdyf(pbcObj)
  rootID <- as.numeric(substr(root,2,nchar(root)))
  BINMAT <- getBINMAT(pbcObj)
  model <- getModel(pbcObj)
  type = 0
  
  
  # ajout Gildas
  numberfunc = dim(BINMAT)[2]
  switch(model,
         "gumbel" = {copulafamily = gumbelCopula()},
         "fgm" = {copulafamily = fgmCopula()},
         "frank" = {copulafamily = frankCopula()},
         "normal" = {copulafamily = normalCopula()},
         "amh" = {copulafamily = amhCopula()},
         "joe" = {copulafamily = joeCopula()},
         "given" = {copulafamily = NA},
         {print("Model not found!")}
  )
  if(!model=="given"){
    for(j in 1:numberfunc){
      if(theta[j] < copulafamily@param.lowbnd | copulafamily@param.upbnd < theta[j] ){
        warning("theta is out of bounds")
      }
    }  
  }
  # fin ajout Gildas
  
  
  
  
  
  
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
  margins <- expression(x^(1/n))
  dxg <- expression(x^(1/n-1)/n)
  out <- 0
  if (output == "density")
    out <- 0
  if (output == "gradient")
    out <- 1
  if (output == "both")
    out <- 2
  ret = .Call('pbc', PACKAGE = 'PBC', BINMAT, theta, x, rootID, f, nIteration, dxf, dxdyf, graf, gradxf, gradxdyf, type, margins, dxg, out)
  if ((output == "density") || (output == "both"))
    pbcTrained <- setDensity(pbcObj, ret[1])
  if ((output == "gradient") || (output == "both")) {
    ## For student copula
    if (type == 8){
      ret1 <- numeric((length(ret)-1)/2)
      for (i in 1:length(ret1))
        ret1[i] <- ret[2 * i]
      if (output == "both")
        pbcTrained <- setGradient(pbcTrained, ret1)
      else 
        pbcTrained <- setGradient(pbcObj, ret1)
    }    
    else {
      if (output == "both")
        pbcTrained <- setGradient(pbcTrained, ret[2:length(ret)])
      else 
        pbcTrained <- setGradient(pbcObj, ret[2:length(ret)])
    }
  }
  return (pbcTrained)
}
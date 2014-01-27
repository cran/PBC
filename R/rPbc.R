## Random Number Generation for Cumulative Distribution Networks
rPBC <- function(n, theta, pbcObj, ...){
  
  
  BINMAT <- pbcObj@BINMAT
  d = dim(BINMAT)[1] # number of variables
  numberfunc = dim(BINMAT)[2] # number of functions
  u = matrix(ncol = d, nrow = n)
  model <- pbcObj@model
  switch(model,
         "gumbel" = {copulafamily = gumbelCopula()},
         "fgm" = {copulafamily = fgmCopula()},
         "frank" = {copulafamily = frankCopula()},
         "normal" = {copulafamily = normalCopula()},
         "amh" = {copulafamily = amhCopula()},
         "joe" = {copulafamily = joeCopula()},
         {print("Model not found!")}
    )
  for(N in 1:n){
    M = BINMAT  
    # fill the matrix M
    for(j in 1:numberfunc){
      
      
      if(theta[j] < copulafamily@param.lowbnd | copulafamily@param.upbnd < theta[j] ){
        stop("theta is out of bounds")
      }
      
      copulafamily@parameters = theta[j]
      M[M[,j] == 1,j] = rCopula(1, copulafamily)
    }
    # generate u via lemma 2.1 from Liebscher 2008
    for(i in 1:d){
      neighbors = which(BINMAT[i,] == 1) # function neighbors of variable i
      numberneighbors = length(neighbors) # number of neighbors
      u[N,i] = max(M[i,neighbors]^numberneighbors) # i-th variable generation
    }
  }
  return(u)
}

## Compute distribution
pPBC <- function(u, theta, pbcObj, ...){
  BINMAT <- pbcObj@BINMAT
  
  if(!is.matrix(u)){
    dim(u) <- c(1,length(u))
  }
  
  
  n = dim(u)[1] # number of observations
  d = length(theta) # dimension of parameters (or number of functions)
  ret <- numeric(n)
  for (i in 1:n)
    ret[i] = 1
  model <- pbcObj@model
  switch(model,
         "gumbel" = {copulafamily = gumbelCopula()},
         "fgm" = {copulafamily = fgmCopula()},
         "frank" = {copulafamily = frankCopula()},
         "normal" = {copulafamily = normalCopula()},
         "amh" = {copulafamily = amhCopula()},
         "joe" = {copulafamily = joeCopula()},
         {print("Model not found!")}
  )
  for (i in 1:n){
    for (j in 1:d){
      
      if(theta[j] < copulafamily@param.lowbnd | copulafamily@param.upbnd < theta[j] ){
        stop("theta is out of bounds")
      }
      
      copulafamily@parameters = theta[j]
      neighbors = which(BINMAT[,j] == 1) # variable neighbors of function j
      # number of neighbor functions of variable 1
      nbNeighbors1 = length(which(BINMAT[neighbors[1],] == 1))
      # number of neighbor functions of variable 2
      nbNeighbors2 = length(which(BINMAT[neighbors[2],] == 1))
      ret[i] = ret[i] * pCopula(c(u[i,neighbors[1]]^(1/nbNeighbors1),u[i,neighbors[2]]^(1/nbNeighbors2)),
                                copulafamily)
    }
  }
  return (ret)
}

## Compute density
dPBC <- function(u, theta, pbcObj, ...){
  
  if(!is.matrix(u)){
    dim(u) <- c(1,length(u))
  }
  
  n = dim(u)[1] # number of observations
  d = dim(u)[2] # number of variables
  ret <- numeric(n)
  for (i in 1:n)
    ret[i] <- mpAlgo(pbcObj, u[i,], theta, output = "density")@density
  return (ret)
}
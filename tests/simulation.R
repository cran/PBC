require(PBC)

isKept = function(x){
  if(sum(x<=u)==5){
    out = 1
  }else{
    out = 0
  }
  return(out)
}

###########################
## Testing Gumbel Copula ##
###########################
theta = runif(4)
u <- matrix(ncol = 5, nrow = 1)
u[1,] = runif(5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCGumbel <- pbcGumbel(g)
n=1000
pbcDataGumbel = rPBC(n, 1/theta, myPBCGumbel)
whoAreKept = apply(pbcDataGumbel,1,isKept)
FGumbel = pPBC(u, 1/theta, myPBCGumbel)

## Comparision
sum(whoAreKept)/n
FGumbel


##########################
## Testing Frank Copula ##
##########################
theta = runif(4)
u <- matrix(ncol = 5, nrow = 1)
u[1,] = runif(5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCFrank <- pbcFrank(g)
n=1000
pbcDataFrank = rPBC(n, theta, myPBCFrank)
whoAreKept = apply(pbcDataFrank,1,isKept)
FFrank = pPBC(u, theta, myPBCFrank)

## Comparision
sum(whoAreKept)/n
FFrank



#############################
## Testing Gaussian Copula ##
#############################
theta = runif(4)
u <- matrix(ncol = 5, nrow = 1)
u[1,] = runif(5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCNormal <- pbcNormal(g)
n=1000
pbcDataNormal = rPBC(n, theta, myPBCNormal)
whoAreKept = apply(pbcDataNormal,1,isKept)
FNormal = pPBC(u, theta, myPBCNormal)

## Comparision
sum(whoAreKept)/n
FNormal



##############################
##### Testing JOE Copula #####
##############################
theta = runif(4)
u <- matrix(ncol = 5, nrow = 1)
u[1,] = runif(5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCJoe <- pbcJoe(g)
n=1000
pbcDataJoe = rPBC(n, 1/theta, myPBCJoe)
whoAreKept = apply(pbcDataJoe,1,isKept)
FJoe = pPBC(u, 1/theta, myPBCJoe)

## Comparision
sum(whoAreKept)/n
FJoe



##############################
##### Testing AMH Copula #####
##############################
theta = runif(4)
u <- matrix(ncol = 5, nrow = 1)
u[1,] = runif(5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCAMH <- pbcAMH(g)
n=1000
pbcDataAMH = rPBC(n, theta, myPBCAMH)
whoAreKept = apply(pbcDataAMH,1,isKept)
FAMH = pPBC(u, theta, myPBCAMH)

## Comparision
sum(whoAreKept)/n
FAMH

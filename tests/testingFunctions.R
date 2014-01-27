## compare symbolic and 'myPBC' densities and gradient at data 'x' and parameter 'theta'
## for various precisions (i.e. the number of digits)

testingPBCPackage = function(N=5,x=NULL,theta=NULL,numberVariables=3,
                             f.exp,deltaf.exp,n1=1,n2=2,n3=1,n4,n5,
                             myPBC, model){
  
  numberEdges = numberVariables-1
  
if(is.null(x)==TRUE){
  x = matrix(nrow=N,ncol=numberVariables,runif(N*numberVariables))
}

if(is.null(theta)==TRUE){
  if ((model == "amh") || (model == "frank") || (model == "fgm"))
    theta = matrix(nrow=N,ncol=numberEdges,runif(N*numberEdges))
  if ((model == "joe") || (model == "gumbel"))
    theta = matrix(nrow=N,ncol=numberEdges,1/runif(N*numberEdges))
}

myDens = c()
packageDens = c()
myGrad = matrix(nrow = N, ncol = numberEdges)
packageGrad = matrix(nrow = N, ncol = numberEdges)
  
for(i in 1:N){
  thisTheta = theta[i,]
  #print(paste("theta = ",theta,sep=""))
  theta1 = thisTheta[1]
  theta2 = thisTheta[2]
  theta3 = thisTheta[3]
  theta4 = thisTheta[4]
  thisX = x[i,]
  #print(paste("x = ",x,sep=""))
  x1 = thisX[1]
  x2 = thisX[2]
  x3 = thisX[3]
  x4 = thisX[4]
  x5 = thisX[5]
  myDens[i] = eval(f.exp)
  myGrad[i,] = attr(eval(deltaf.exp),"gradient")
  packageOut = mpAlgo(myPBC, u=thisX, theta=thisTheta)
  packageDens[i] = packageOut@density
  packageGrad[i,] = packageOut@gradient
}


## Is gradient good? Up to which precision? 
isEqualGrad = c()
for(precision in 1:20){
  howManyEqual = sum(round(myGrad,precision)==round(packageGrad,precision))
  isEqualGrad[precision] = (howManyEqual==(2*N))
}


## Is density good? Up to which precision?
isEqualDensity = c()
for(precision in 1:20){
  howManyEqual = sum(round(myDens,precision)==round(packageDens,precision))
  isEqualDensity[precision] = (howManyEqual==N)
}


return(list(symbolicDensity=myDens,packageDensity=packageDens,
            symbolicGradient=myGrad,packageGradient=packageGrad,
            equalDensities=isEqualDensity,equalGradients=isEqualGrad))

}

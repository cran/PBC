require(PBC)
source("./testingFunctions.R")

densPrecision <- function(testing){
  ret <- numeric(20)
  for(preci in 1:20){
    dd = abs(testing$symbolicDensity-testing$packageDensity)
    dd1 = dd[!is.nan(dd)]
    dd2 = dd1[!is.infinite(dd1)]
    ret[preci] = round(sum(dd2),preci)
  }
  return (ret)
}

graPrecision <- function(testing){
  ret <- numeric(20)
  for(preci in 1:20){
    dd = abs(testing$symbolicGradient-testing$packageGradient)
    dd1 = dd[!is.nan(dd)]
    dd2 = dd1[!is.infinite(dd1)]
    #dd3 = apply(abs(dd2),1,sum) # norm '1'
    ret[preci] = round(sum(dd2),preci)
  }
  return (ret)
}

###################################
## TESTING GUMBEL COPULA FOR d=3 ##
###################################

## package Gumbel copula (d=3)

g <- graph.formula(X1-X2,X2-X3,simplify = FALSE)
myPBCGumbel <- pbc(g, model="gumbel")
pbcPlot(myPBCGumbel)

## symbolic Gumbel copula (d=3)

F.exp <- expression(
  
  exp(-
        (
          (-log(x1)/n1)^(theta1) +
            (-log(x2)/n2)^(theta1)
        )^(1/theta1)
  )
  *
    exp(-
          (
            (-log(x2)/n2)^(theta2) +
              (-log(x3)/n3)^(theta2)
          )^(1/theta2)
    ) 
)

f.exp <- D(D(D(F.exp,"x1"),"x2"),"x3")
withRespectTo <- c("theta1","theta2")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,f.exp=f.exp,deltaf.exp=deltaf.exp,
                            myPBC=myPBCGumbel, model="gumbel")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

###################################
## TESTING GUMBEL COPULA FOR d=5 ##
###################################

## package Gumbel copula (d=5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCGumbel <- pbc(g, model="gumbel")
pbcPlot(myPBCGumbel)

## symbolic Gumbel copula 

F.exp <- expression(
  
  exp(-
        (
          (-log(x1)/n1)^(theta1) +
            (-log(x4)/n4)^(theta1)
        )^(1/theta1)
  )
  *
    exp(-
          (
            (-log(x4)/n4)^(theta2) +
              (-log(x2)/n2)^(theta2)
          )^(1/theta2)
    ) 
  *
    exp(-
          (
            (-log(x2)/n2)^(theta3) +
              (-log(x3)/n3)^(theta3)
          )^(1/theta3)
    )
  *
    exp(-
          (
            (-log(x4)/n4)^(theta4) +
              (-log(x5)/n5)^(theta4)
          )^(1/theta4)
    ) 
)

f.exp <- D(D(D(D(D(F.exp,"x1"),"x2"),"x3"),"x4"),"x5")
withRespectTo <- c("theta1","theta2","theta3","theta4")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,numberVariables=5,f.exp=f.exp,
                            deltaf.exp=deltaf.exp,myPBC=myPBCGumbel,
                            n1=1,n2=2,n3=1,n4=3,n5=1, model="gumbel")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")



##########################################################################################################
####################################### END OF GUMBEL COPULA TESTS #######################################
##########################################################################################################



################################
## TESTING FGM COPULA FOR d=3 ##
################################

## package FGM copula (d=3)

g <- graph.formula(X1-X2,X2-X3,simplify = FALSE)
myPBCFGM <- pbc(g, model="fgm")
pbcPlot(myPBCFGM)

## symbolic FGM copula (d=3)

F.exp <- expression(
  
  (x1^(1/n1) * x2^(1/n2) + theta1 * x1^(1/n1) * x2^(1/n2) * (1 - x1^(1/n1)) * (1 - x2^(1/n2)))
  *
    (x3^(1/n3) * x2^(1/n2) + theta2 * x3^(1/n3) * x2^(1/n2) * (1 - x3^(1/n3)) * (1 - x2^(1/n2)))
  
)

f.exp <- D(D(D(F.exp,"x1"),"x2"),"x3")
withRespectTo <- c("theta1","theta2")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,f.exp=f.exp,deltaf.exp=deltaf.exp,
                            myPBC=myPBCFGM, model="fgm")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")


################################
## TESTING FGM COPULA FOR d=5 ##
################################

## package FGM copula (d=5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCFGM <- pbc(g, model="fgm")
pbcPlot(myPBCFGM)

## symbolic FGM copula 

F.exp <- expression(
  
  (x1^(1/n1) * x4^(1/n4) + theta1 * x1^(1/n1) * x4^(1/n4) * (1 - x1^(1/n1)) * (1 - x4^(1/n4)))
  *
    (x4^(1/n4) * x2^(1/n2) + theta2 * x4^(1/n4) * x2^(1/n2) * (1 - x4^(1/n4)) * (1 - x2^(1/n2)))
  *
    (x3^(1/n3) * x2^(1/n2) + theta3 * x3^(1/n3) * x2^(1/n2) * (1 - x3^(1/n3)) * (1 - x2^(1/n2)))
  *
    (x4^(1/n4) * x5^(1/n5) + theta4 * x4^(1/n4) * x5^(1/n5) * (1 - x4^(1/n4)) * (1 - x5^(1/n5)))
  
)

f.exp <- D(D(D(D(D(F.exp,"x1"),"x2"),"x3"),"x4"),"x5")
withRespectTo <- c("theta1","theta2","theta3","theta4")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,numberVariables=5,f.exp=f.exp,
                            deltaf.exp=deltaf.exp,myPBC=myPBCFGM,
                            n1=1,n2=2,n3=1,n4=3,n5=1, model="fgm")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")



##########################################################################################################
######################################## END OF FGM COPULA TESTS #########################################
##########################################################################################################



##################################
## TESTING FRANK COPULA FOR d=3 ##
##################################

## package Frank copula (d=3)

g <- graph.formula(X1-X2,X2-X3,simplify = FALSE)
myPBCFrank <- pbc(g, model="frank")
pbcPlot(myPBCFrank)
f <- expression(-log(1+((exp(-theta*x)-1)*(exp(-theta*y)-1))/(exp(-theta)-1))/theta)
## symbolic Frank copula (d=3)

F.exp <- expression(
  
  (-log(1+((exp(-theta1*x1^(1/n1))-1)*(exp(-theta1*x2^(1/n2))-1))/(exp(-theta1)-1))/theta1)
  *
    (-log(1+((exp(-theta2*x2^(1/n2))-1)*(exp(-theta2*x3^(1/n3))-1))/(exp(-theta2)-1))/theta2)
  
)

f.exp <- D(D(D(F.exp,"x1"),"x2"),"x3")
withRespectTo <- c("theta1","theta2")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,f.exp=f.exp,deltaf.exp=deltaf.exp,
                            myPBC=myPBCFrank, model="frank")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")


##################################
## TESTING FRANK COPULA FOR d=5 ##
##################################

## package Frank copula (d=5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCFrank <- pbc(g, model="frank")
pbcPlot(myPBCFrank)

## symbolic Frank copula 

F.exp <- expression(
  
  (-log(1+((exp(-theta1*x1^(1/n1))-1)*(exp(-theta1*x4^(1/n4))-1))/(exp(-theta1)-1))/theta1)
  *
    (-log(1+((exp(-theta2*x2^(1/n2))-1)*(exp(-theta2*x4^(1/n4))-1))/(exp(-theta2)-1))/theta2)
  *
    (-log(1+((exp(-theta3*x2^(1/n2))-1)*(exp(-theta3*x3^(1/n3))-1))/(exp(-theta3)-1))/theta3)
  *
    (-log(1+((exp(-theta4*x4^(1/n4))-1)*(exp(-theta4*x5^(1/n5))-1))/(exp(-theta4)-1))/theta4)
  
)

f.exp <- D(D(D(D(D(F.exp,"x1"),"x2"),"x3"),"x4"),"x5")
withRespectTo <- c("theta1","theta2","theta3","theta4")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,numberVariables=5,f.exp=f.exp,
                            deltaf.exp=deltaf.exp,myPBC=myPBCFrank,
                            n1=1,n2=2,n3=1,n4=3,n5=1, model="frank")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")




##########################################################################################################
####################################### END OF FRANK COPULA TESTS ########################################
##########################################################################################################



##################################
### TESTING JOE COPULA FOR d=3 ###
##################################

## package Joe copula (d=3)

g <- graph.formula(X1-X2,X2-X3,simplify = FALSE)
myPBCJoe <- pbc(g, model="joe")
pbcPlot(myPBCJoe)

## symbolic Joe copula (d=3)

F.exp <- expression(
  
  (1-((1-x1^(1/n1))^theta1+(1-x2^(1/n2))^theta1-(1-x1^(1/n1))^theta1*(1-x2^(1/n2))^theta1)^(1/theta1))
  *
    (1-((1-x3^(1/n3))^theta2+(1-x2^(1/n2))^theta2-(1-x3^(1/n3))^theta2*(1-x2^(1/n2))^theta2)^(1/theta2))
  
)

f.exp <- D(D(D(F.exp,"x1"),"x2"),"x3")
withRespectTo <- c("theta1","theta2")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,f.exp=f.exp,deltaf.exp=deltaf.exp,
                            myPBC=myPBCJoe, model="joe")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")


##################################
### TESTING JOE COPULA FOR d=5 ###
##################################

## package Joe copula (d=5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCJoe <- pbc(g, model="joe")
pbcPlot(myPBCJoe)

## symbolic Joe copula 

F.exp <- expression(
  
  (1-((1-x1^(1/n1))^theta1+(1-x4^(1/n4))^theta1-(1-x1^(1/n1))^theta1*(1-x4^(1/n4))^theta1)^(1/theta1))
  *
    (1-((1-x4^(1/n4))^theta2+(1-x2^(1/n2))^theta2-(1-x4^(1/n4))^theta2*(1-x2^(1/n2))^theta2)^(1/theta2))
  *
    (1-((1-x3^(1/n3))^theta3+(1-x2^(1/n2))^theta3-(1-x3^(1/n3))^theta3*(1-x2^(1/n2))^theta3)^(1/theta3))
  *
    (1-((1-x4^(1/n4))^theta4+(1-x5^(1/n5))^theta4-(1-x4^(1/n4))^theta4*(1-x5^(1/n5))^theta4)^(1/theta4))  
  
)

f.exp <- D(D(D(D(D(F.exp,"x1"),"x2"),"x3"),"x4"),"x5")
withRespectTo <- c("theta1","theta2","theta3","theta4")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,numberVariables=5,f.exp=f.exp,
                            deltaf.exp=deltaf.exp,myPBC=myPBCJoe,
                            n1=1,n2=2,n3=1,n4=3,n5=1, model="joe")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")




##########################################################################################################
####################################### END OF JOE COPULA TESTS #########################################
##########################################################################################################



##################################
### TESTING AMH COPULA FOR d=3 ###
##################################

## package AMH copula (d=3)

g <- graph.formula(X1-X2,X2-X3,simplify = FALSE)
myPBCAMH <- pbc(g, model="amh")
pbcPlot(myPBCAMH)

## symbolic AMH copula (d=3)

F.exp <- expression(
  
  x1^(1/n1) * x2^(1/n2) / (1 - theta1 * (1 - x1^(1/n1)) * (1 - x2^(1/n2)))
  *
    x3^(1/n3) * x2^(1/n2) / (1 - theta2 * (1 - x3^(1/n3)) * (1 - x2^(1/n2)))
  
)

f.exp <- D(D(D(F.exp,"x1"),"x2"),"x3")
withRespectTo <- c("theta1","theta2")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,f.exp=f.exp,deltaf.exp=deltaf.exp,
                            myPBC=myPBCAMH, model="amh")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")


##################################
### TESTING AMH COPULA FOR d=5 ###
##################################

## package AMH copula (d=5)

g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)
myPBCAMH <- pbc(g, model="amh")
pbcPlot(myPBCAMH)

## symbolic AMH copula 

F.exp <- expression(
  
  x1^(1/n1) * x4^(1/n4) / (1 - theta1 * (1 - x1^(1/n1)) * (1 - x4^(1/n4)))
  *
    x4^(1/n4) * x2^(1/n2) / (1 - theta2 * (1 - x4^(1/n4)) * (1 - x2^(1/n2)))
  *
    x3^(1/n3) * x2^(1/n2) / (1 - theta3 * (1 - x3^(1/n3)) * (1 - x2^(1/n2)))
  *
    x4^(1/n4) * x5^(1/n5) / (1 - theta4 * (1 - x4^(1/n4)) * (1 - x5^(1/n5)))
)

f.exp <- D(D(D(D(D(F.exp,"x1"),"x2"),"x3"),"x4"),"x5")
withRespectTo <- c("theta1","theta2","theta3","theta4")
deltaf.exp <- deriv(f.exp,withRespectTo)

## testing

testing = testingPBCPackage(N=10,numberVariables=5,f.exp=f.exp,
                            deltaf.exp=deltaf.exp,myPBC=myPBCAMH,
                            n1=1,n2=2,n3=1,n4=3,n5=1, model="amh")

# discrepancy between the density symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- densPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")

# discrepancy between the gradient symbolic and package calculations
# evaluated at N random points (see testingPBCPackage)

yy <- graPrecision(testing)
plot(yy,xlab="number of digits after the dot", ylab="sum of the norms")




##########################################################################################################
####################################### END OF AMH COPULA TESTS #########################################
##########################################################################################################

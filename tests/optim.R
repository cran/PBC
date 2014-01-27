require(PBC)

################
## Parameters ##
################
theta <- runif(4)
g <- graph.formula(X1-X4,X4-X2,X2-X3,X4-X5,simplify = FALSE)

##########################
## PBC objects Creation ##
##########################
myPBCGumbel <- pbcGumbel(g)
myPBCFGM <- pbcFGM(g)
myPBCFrank <- pbcFrank(g)
myPBCNormal <- pbcNormal(g)
myPBCAMH <- pbcAMH(g)
myPBCJoe <- pbcJoe(g)

#####################
## Data Generation ##
#####################
n=100
pbcDataGumbel = rPBC(n, 1/theta, myPBCGumbel)
pbcDataFGM = rPBC(n, theta, myPBCFGM)
pbcDataFrank = rPBC(n, theta, myPBCFrank)
pbcDataNormal = rPBC(10, theta, myPBCNormal)
pbcDataAMH = rPBC(n, theta, myPBCAMH)
pbcDataJoe = rPBC(n, 1/theta, myPBCJoe)

###################################
## Maximum Likelihood Estimation ##
###################################
thetaGumbel <- pbcOptim(rep(2, 4), pbcDataGumbel, myPBCGumbel, method = 'BFGS')
thetaFGM <- pbcOptim(rep(0.5, 4), pbcDataFGM, myPBCFGM, method = 'L-BFGS-B',
                     lower = rep(0, 4), upper = rep(1, 4))
thetaFrank <- pbcOptim(rep(1, 4), pbcDataFrank, myPBCFrank, method = 'BFGS')
thetaNormal <- pbcOptim(rep(0.5, 4), pbcDataNormal, myPBCNormal, method = 'L-BFGS-B',
                        lower = rep(0, 4), upper = rep(0.99, 4))
thetaAMH <- pbcOptim(rep(0.5, 4), pbcDataAMH, myPBCAMH, method = 'L-BFGS-B',
                     lower = rep(0, 4), upper = rep(0.99, 4))
thetaJoe <- pbcOptim(rep(2, 4), pbcDataJoe, myPBCJoe, method = 'BFGS')



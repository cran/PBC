setOldClass("igraph")

###################
### Class 'PBC' ###
###################

setClass(
  Class="PBC",
  representation=representation(
    root = "character",
    graph = "igraph",  
    f = "expression",
    dxf = "expression",
    dxdyf = "expression",
    graf = "expression",
    gradxf = "expression",
    gradxdyf = "expression",
    BINMAT = "matrix",
    model = "character",
    density = "numeric",
    gradient = "vector",
    nIteration = "numeric"
  )
)

###############
### Methods ###
###############

### Constructor ###

setMethod(
  f="initialize",
  signature="PBC",
  definition=function(.Object, root, g, f, nIteration){
    .Object@root <- root
    #.Object@graph <- graph
    .Object@f <- f
    .Object@nIteration <- nIteration
    .Object@dxf <- as.expression(D(f,"x"))
    .Object@dxdyf <- as.expression(D(D(f,"x"),"y"))
    .Object@graf <- as.expression(D(f,"theta"))
    .Object@gradxf <- as.expression(D(D(f,"x"),"theta"))
    .Object@gradxdyf <- as.expression(D(D(D(f,"x"),"y"),"theta"))
    #objGraph <- new("Graph", graph, root)
    ## Get root number
    #.Object@rootID <- getRootID(objGraph)
    #.Object@rootID <- as.numeric(substr(root,2,nchar(root)))
    ## Get binary matrix
    #.Object@BINMAT <- getBinMat(objGraph)
    #.Object@type <- 0
    #if (paste(f)=="0")
    #  .Object@type <- 1
    .Object@model <- "given"
    if (paste(f)=="1")
      .Object@model <- "gumbel"
    if (paste(f)=="2")
      .Object@model <- "fgm"
    if (paste(f)=="3")
      .Object@model <- "frank"
    if (paste(f)=="4")
      .Object@model <- "normal"
    if (paste(f)=="5")
      .Object@model <- "amh"
    if (paste(f)=="6")
      .Object@model <- "joe"
    #if (paste(f)=="7")
    #  .Object@model <- 8
    #compute <- function(x, y, theta, h){
    #  return (eval(h))
    #}
    
    .Object@graph <- graph.empty(n=0, directed=FALSE)
    edgeList <- get.edgelist(g)
    nVariable <- length(V(g))
    nFunction <- length(E(g))
    ## Draw diamond shape for functions base on a example from igraph tutorial
    diamond <- function(coords, v=NULL, params) {
      ## Set color
      vertex.color <- params("vertex", "color")
      if (length(vertex.color) != 1 && !is.null(v)) {
        vertex.color <- vertex.color[v]
      }
      ## Set size
      vertex.size <- 1/200 * params("vertex", "size")
      if (length(vertex.size) != 1 && !is.null(v)) {
        vertex.size <- vertex.size[v]
      }
      ## Set shape
      symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
              stars=cbind(vertex.size, vertex.size, vertex.size, vertex.size),
              add=TRUE, inches=FALSE)
    }    
    add.vertex.shape("diamond", clip=vertex.shapes("circle")$clip,
                     plot=diamond)
    ## Add variable vertex to full graph
    for (i in 1:nVariable){
      #xi <- paste("X", i, sep="")
      xi <- get.vertex.attribute(g, 'name')[i]
      if (xi != .Object@root) ## White color for not-root
        .Object@graph <- .Object@graph + vertex(xi, shape = "circle", color = "white", size = 10)
      else ## Red color for root
        .Object@graph <- .Object@graph + vertex(xi, shape = "circle", color = "green", size = 10)
    }
    ## Add function vertex and adjacent edges to full graph
    for (i in 1:nFunction){
      #Fi <- paste(expression(Phi), i, sep = "")
      if ((i == 1)&&(nVariable < 21))
        Fi <- expression(paste(Phi,1))
      if ((i == 2)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 2))
      if ((i == 3)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 3))
      if ((i == 4)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 4))
      if ((i == 5)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 5))
      if ((i == 6)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 6))
      if ((i == 7)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 7))
      if ((i == 8)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 8))
      if ((i == 9)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 9))
      if ((i == 10)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 10))
      if ((i == 11)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 11))
      if ((i == 12)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 12))
      if ((i == 13)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 13))
      if ((i == 14)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 14))
      if ((i == 15)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 15))
      if ((i == 16)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 16))
      if ((i == 17)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 17))
      if ((i == 18)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 18))
      if ((i == 19)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 19))
      if ((i == 20)&&(nVariable < 21))
        Fi <- expression(paste(Phi, 20))
      if (nVariable > 20)
        Fi <- paste("Phi", i)
      ## Add function vertex
      .Object@graph <- .Object@graph + vertex(Fi,  shape = "diamond", color = "blue", size = 10)
      ## Add adjacent edges
      #.Object@graph <- .Object@graph + edge(Fi, edgeList[i])
      #.Object@graph <- .Object@graph + edge(Fi, edgeList[i + nFunction]) 
      .Object@graph <- .Object@graph + edge(substr(Fi,0,99), edgeList[i])
      .Object@graph <- .Object@graph + edge(substr(Fi,0,99), edgeList[i + nFunction]) 
    }
    ## Get binary matrix from full graph
    nEdge <- length(get.edgelist(.Object@graph))/2
    #.Object@BINMAT <-  mat.or.vec(nVariable, nFunction)
    .Object@BINMAT <-  matrix(data=0,nrow=nVariable,ncol=nFunction)
    for (i in 1:nEdge){  
      ## Get matrix element position whose value is equal to 1
      vertexID <- get.edges(.Object@graph, E(.Object@graph))[i + nEdge]
      #pos_hor = as.numeric(str_extract(V(.Object@graph)[vertexID]$name,"\\d"))
      nameID <- V(.Object@graph)[vertexID]$name
      pos_hor = as.numeric(substr(nameID,2,nchar(nameID)))
      #pos_hor = get.edges(.Object@graph, E(.Object@graph))[i + nEdge]
      pos_ver = get.edges(.Object@graph, E(.Object@graph))[i] - nVariable
      .Object@BINMAT[(pos_ver - 1) * nVariable + pos_hor] = 1
    }
    return(.Object) 
  }
)


### Plot ###

setGeneric(name = "pbcPlot",
           def = function(object){standardGeneric ("pbcPlot")})
setMethod("pbcPlot","PBC",
          function(object){
            plot(object@graph, xlim=c(-1,1), ylim=c(-1,1), vertex.label.cex=0.5, vertex.label.color="red")
          }
)


###############
### Getters ###
###############

setGeneric(name = "getRoot",
           def = function(object){standardGeneric ("getRoot")})
setMethod("getRoot","PBC",
          function(object){
            return(object@root)
          }
)

setGeneric(name = "getGraph",
           def = function(object){standardGeneric ("getGraph")})
setMethod("getGraph","PBC",
          function(object){
            return(object@graph)
          }
)

setGeneric(name = "getF",
           def = function(object){standardGeneric ("getF")})
setMethod("getF","PBC",
          function(object){
            return(object@f)
          }
)

setGeneric(name = "getNIteration",
           def = function(object){standardGeneric ("getNIteration")})
setMethod("getNIteration","PBC",
          function(object){
            return(object@nIteration)
          }
)

setGeneric(name = "getDxf",
           def = function(object){standardGeneric ("getDxf")})
setMethod("getDxf","PBC",
          function(object){
            return(object@dxf)
          }
)

setGeneric(name = "getDxdyf",
           def = function(object){standardGeneric ("getDxdyf")})
setMethod("getDxdyf","PBC",
          function(object){
            return(object@dxdyf)
          }
)

setGeneric(name = "getGraf",
           def = function(object){standardGeneric ("getGraf")})
setMethod("getGraf","PBC",
          function(object){
            return(object@graf)
          }
)

setGeneric(name = "getGradxf",
           def = function(object){standardGeneric ("getGradxf")})
setMethod("getGradxf","PBC",
          function(object){
            return(object@gradxf)
          }
)

setGeneric(name = "getGradxdyf",
           def = function(object){standardGeneric ("getGradxdyf")})
setMethod("getGradxdyf","PBC",
          function(object){
            return(object@gradxdyf)
          }
)


setGeneric(name = "getBINMAT",
           def = function(object){standardGeneric ("getBINMAT")})
setMethod("getBINMAT","PBC",
          function(object){
            return(object@BINMAT)
          }
)

setGeneric(name = "getModel",
           def = function(object){standardGeneric ("getModel")})
setMethod("getModel","PBC",
          function(object){
            return(object@model)
          }
)

###############
### Setters ###
###############

setGeneric(name = "setDensity",
           def = function(object, density){standardGeneric ("setDensity")})
setMethod("setDensity","PBC",
          function(object, density){
            object@density <- density
            return(object)
          }
)

setGeneric(name = "setGradient",
           def = function(object, gradient){standardGeneric ("setGradient")})
setMethod("setGradient","PBC",
          function(object, gradient){
            object@gradient <- gradient
            return(object)
          }
)


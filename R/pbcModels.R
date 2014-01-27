## Pre-existing models for Cumulative Distribution Networks

## Gumbel copula
pbcGumbel <- function(graph)
{
  .Object <- pbc(graph, model="gumbel")
  return (.Object)
}

## Farlie-Gumbel-Morgenstern (FGM) copula
pbcFGM <- function(graph)
{
  .Object <- pbc(graph, model="fgm")
  return (.Object)
}

## Frank Copula
pbcFrank <- function(graph)
{
  .Object <- pbc(graph, model="frank")
  return (.Object)
}

## Normal copula
pbcNormal <- function(graph)
{
  .Object <- pbc(graph, model="normal")
  return (.Object)
}

## Ali-Mikhail-Haq (AMH) copula
pbcAMH <- function(graph)
{
  .Object <- pbc(graph, model="amh")
  return (.Object)
}

## Joe copula
pbcJoe <- function(graph)
{
  .Object <- pbc(graph, model="joe")
  return (.Object)
}
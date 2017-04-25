source("genalgCustom.R")
source("SA.R")

dataInput <- read.csv(file="1_dataSetSAPS2.csv",header = TRUE, sep=",")

getPSample <- function(){
  
  #x = dataInput[, c(2:ncol(dataInput))]
  #y = dataInput[, c(1)]
  
  #initialization of variables
  detP = 0
  nFeatures = ncol(dataInput)-1
  
  #creation of a symetric matrix with det>0
  while (detP <= 0){
    pMatrix = matrix(runif(nFeatures^2, min = 0, max = 1), nrow = nFeatures, ncol = nFeatures)
    pMatrix[lower.tri(pMatrix)] = t(pMatrix)[lower.tri(pMatrix)]
    #pMatrix = pMatrix %*% t(pMatrix)
    diag(pMatrix) = 1
    detP = det(pMatrix)
  }
  pMatrix[lower.tri(pMatrix, diag = TRUE)] = NA
  pSample = as.vector(t(pMatrix))
  pSample = pSample[!is.na(pSample)]
  return(pSample)
}

objectiveFunction <- function(solution){
  
  pMatrix = matrix(1L, nrow = 80, ncol = 80, byrow = TRUE)
  pMatrix[lower.tri(pMatrix, diag = FALSE)] <- solution
  pMatrix = t(pMatrix)
  pMatrix[lower.tri(pMatrix)] = t(pMatrix)[lower.tri(pMatrix)]
  x = c(2:ncol(dataInput))
  
  # parameter: length of neighborhoods
  k = 5
  
  compactNeighborhood = matrix(0L, nrow = nrow(dataInput))
  scatterNeighborhood = matrix(0L, nrow = nrow(dataInput))
  for (i in 1:21){
    neighbour =1
    xi = as.numeric(dataInput[i,x])
    generalMahalanobis = matrix(data = 0L, nrow = nrow(dataInput), ncol = 2)
    for (j in 1:21){
      if (i != j){
        xj = as.numeric(dataInput[j,x])
        generalMahalanobis[j,1] = t((xi - xj)) %*% pMatrix %*% (xi-xj)
        generalMahalanobis[j,2] = dataInput[j,1]
      }
    }
    #orderedMahalanobis = matrix(-1L, nrow = nrow(generalMahalanobis), ncol = ncol(generalMahalanobis))
    orderedMahalanobis = generalMahalanobis[order(generalMahalanobis[,1], decreasing = TRUE),]
    while (neighbour < k*2){
      if (orderedMahalanobis[neighbour, 2] == dataInput[i,1]){
        compactNeighborhood[i] = compactNeighborhood[i] + orderedMahalanobis[neighbour, 1]
        neighbour = neighbour + 1
      }else{
        scatterNeighborhood[i] = scatterNeighborhood[i] + orderedMahalanobis[neighbour, 1]
        neighbour = neighbour + 1
      }
    }
  }
  
  discriminabilityFunction = sum(compactNeighborhood)/sum(scatterNeighborhood)
  return(discriminabilityFunction)
}

geneticAlg = function(){
  #Genetic algorithm
  library(ggplot2)
  
  
  pMatrix <- getPMatrix()
  #minFeature = 0
  #maxFeature = 1
  minGen = rep(0,3160)
  maxGen = rep(1, 3160)
  
  
  
  
  results = rbgaCustom(pSample = getPSample, minGen, maxGen, evalFunc = objectiveFunction, popSize = 6000, iters = 5, mutationChance = 0.01, verbose = FALSE)
  #cat(genalg:::summary.rbga(results))
  genalg:::plot.rbga(results)
  plot(results, type = "hist")
  plot(results, type = "vars")
}


system.time(geneticAlg())
#system.time(simulated_annealing(getPSample, objectiveFunction))
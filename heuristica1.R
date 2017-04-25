source("genalgCustom.R")
source("SA.R")

#Reading dataset
dataInput <- read.csv(file="1_dataSetSAPS2.csv",header = TRUE, sep=",")

#return a vector with the upper triangle matrix elements of P matrix
getPSample <- function(currentSample = NULL){
  
  #initialization of variables
  detP = 0
  nFeatures = ncol(dataInput)-1
  
  #creation of a symetric matrix with det>0
  while (detP <= 0){
    #creating a square matrix with random values between 0, 1 
    pMatrix = matrix(runif(nFeatures^2, min = 0, max = 1), nrow = nFeatures, ncol = nFeatures)
    if (!is.null(currentSample)){
      solution = rnorm(length(currentSample), currentSample, 1)
      pMatrix[lower.tri(pMatrix, diag = FALSE)] <- solution
      pMatrix = t(pMatrix)
      pMatrix[lower.tri(pMatrix)] = t(pMatrix)[lower.tri(pMatrix)]
    }else{
      #setting symetric matrix
      pMatrix[lower.tri(pMatrix)] = t(pMatrix)[lower.tri(pMatrix)]
    }
    
    #setting 1 in main diagonal
    diag(pMatrix) = 1
    #validating determinant > 0
    detP = det(pMatrix)
  }
  #getting the upper triangle matrix from P matrix
  pMatrix[lower.tri(pMatrix, diag = TRUE)] = NA
  pSample = as.vector(t(pMatrix))
  pSample = pSample[!is.na(pSample)]
  return(pSample)
}

objectiveFunction <- function(solution){
  
  #transforming sample in P matrix
  pMatrix = matrix(1L, nrow = 80, ncol = 80, byrow = TRUE)
  pMatrix[lower.tri(pMatrix, diag = FALSE)] <- solution
  pMatrix = t(pMatrix)
  pMatrix[lower.tri(pMatrix)] = t(pMatrix)[lower.tri(pMatrix)]
  x = c(2:ncol(dataInput))
  
  # parameter: length of neighborhoods
  k = 5
  
  #Neighborhoods
  compactNeighborhood = matrix(0L, nrow = nrow(pMatrix))
  scatterNeighborhood = matrix(0L, nrow = nrow(pMatrix))
  for (i in 1:1000){
    neighbour =1
    xi = as.numeric(dataInput[i,x])
    generalMahalanobis = matrix(data = 0L, nrow = nrow(dataInput), ncol = 2)
    for (j in 1:1000){
      if (i != j){
        xj = as.numeric(dataInput[j,x])
        #Mahalanobis distance metric 
        generalMahalanobis[j,1] = t((xi - xj)) %*% pMatrix %*% (xi-xj)
        generalMahalanobis[j,2] = dataInput[j,1]
      }
    }
    orderedMahalanobis = generalMahalanobis[order(generalMahalanobis[,1], decreasing = TRUE),]
    #filling the nieghborhoods
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
  
  
  pMatrix <- getPSample()
  minGen = rep(0,3160)
  maxGen = rep(1, 3160)
  
  
  
  #Execution of genetic algorithms with
  #population size = 3500, iterations = 10, mutation probability = 0.01 and elitism = 20% by defect.
  results = rbgaCustom(pSample = getPSample, minGen, maxGen, evalFunc = objectiveFunction, 
                       popSize = 3500, iters = 10, mutationChance = 0.01, verbose = FALSE)
  #cat(genalg:::summary.rbga(results))
  genalg:::plot.rbga(results)
  plot(results, type = "hist")
  plot(results, type = "vars")
}


#system.time(geneticAlg())
system.time(results <- simulated_annealing(getPSample, objectiveFunction))

#dataInput <- read.csv(file="~/Documentos/Metaheuristics/LocallySupervisedMetricLearningHeuristics/1_dataSetSAPS2.csv",header = TRUE, sep=",")
dataInput <- read.csv("~/LSMLH/1_dataSetSAPS2.csv",header = TRUE, sep=",")

getPMatrix <- function(){
  
  #x = dataInput[, c(2:ncol(dataInput))]
  #y = dataInput[, c(1)]
  
  #initialization of variables
  detP = 0
  x = c(2:ncol(dataInput))
  nFeatures = ncol(dataInput)-1
  
  #creation of a symetric matrix with det>0
  while (detP <= 0){
    pMatrix = matrix(runif(nFeatures^2, min = 0, max = 1), nrow = nFeatures, ncol = nFeatures)
    pMatrix[lower.tri(pMatrix)] = t(pMatrix)[lower.tri(pMatrix)]
    #pMatrix = pMatrix %*% t(pMatrix)
    detP = det(pMatrix)
  }
  return(pMatrix)
}

objectiveFunction <- function(chromosome){
  
  pMatrix = matrix(chromosome, nrow = 80, ncol = 80, byrow = TRUE)
  
  # parameter: length of neighborhoods
  k = 5
  
  compactNeighborhood = matrix(0L, nrow = nrow(dataInput))
  scatterNeighborhood = matrix(0L, nrow = nrow(dataInput))
  for (i in 1:201){
    neighbour =1
    xi = as.numeric(dataInput[i,x])
    generalMahalanobis = matrix(data = 0L, nrow = nrow(dataInput), ncol = 2)
    for (j in 1:201){
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


#Genetic algorithm
library(genalg)
library(ggplot2)


pMatrix <- getPMatrix()
chromosome = as.vector(t(pMatrix))
minFeature = 0
maxFeature = 1
minGen = matrix(0L, nrow = 80, ncol = 80)
maxGen = matrix(1L, nrow = 80, ncol = 80)



results = rbga(minGen, maxGen, evalFunc = objectiveFunction, popSize = 120, iters = 20, mutationChance = 0.01, verbose = TRUE, showSettings = TRUE)
cat(genalg:::summary.rbga(results))
genalg:::plot.rbga(results)
plot(results, type = "hist")
plot(results, type = "vars")


simulated_annealing <- function(objectiveFunction, s0, niter = 10, step = 0.1) {
  
  # Initialize
  ## s stands for state
  ## f stands for function value
  ## b stands for best
  ## c stands for current
  ## n stands for neighbor
  s_b <- s_c <- s_n <- s0
  f_b <- f_c <- f_n <- func(s_n)
  message("It\tBest\tCurrent\tNeigh\tTemp")
  message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", 0L, f_b, f_c, f_n, 1))
  
  for (k in 1:niter) {     
    Temp <- (1 - step)^k
    # consider a random neighbor
    s_n <- rnorm(2, s_c, 1)
    f_n <- func(s_n)
    # update current state
    if (f_n < f_c || runif(1, 0, 1) < exp(-(f_n - f_c) / Temp)) {
      s_c <- s_n
      f_c <- f_n
    }
    # update best state
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n         
    }
    message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", k, f_b, f_c, f_n, Temp))
  }
  return(list(iterations = niter, best_value = f_b, best_state = s_b))
}







 


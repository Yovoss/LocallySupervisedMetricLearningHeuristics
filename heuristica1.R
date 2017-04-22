dataInput <- read.csv(file="~/Documentos/Metaheuristics/LocallySupervisedMetricLearningHeuristics/1_dataSetSAPS2.csv",header = TRUE, sep=",")
#dataInput <- read.csv("~/LSMLH/1_dataSetSAPS2.csv",header = TRUE, sep=",")
#x = dataInput[, c(2:ncol(dataInput))]
#y = dataInput[, c(1)]

# parameter: length of neighborhoods
k = 5

#initialization of variables
detP = 0
x = c(2:ncol(dataInput))
nFeatures = ncol(dataInput)-1

#creation of a symetric matrix with det>0
while (detP <= 0){
  pMatrix = matrix(rexp(nFeatures^2, rate=10), nrow = nFeatures, ncol = nFeatures)
  pMatrix = pMatrix %*% t(pMatrix)
  detP = det(pMatrix)
}

compactNeighborhood = matrix(0L, nrow = nrow(dataInput))
scatterNeighborhood = matrix(0L, nrow = nrow(dataInput))
for (i in 1:11){
  neighbour = 1
  xi = as.numeric(dataInput[i,x])
  generalMahalanobis = matrix(data = NA, nrow = nrow(dataInput), ncol = 2)
  for (j in 1:11){
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

objectiveFunction = sum(compactNeighborhood)/sum(scatterNeighborhood)

 


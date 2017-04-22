#dataInput <- read.csv(file="~/Documentos/Metaheuristics/LocallySupervisedMetricLearningHeuristics/1_dataSetSAPS2.csv",header = TRUE, sep=",")
dataInput <- read.csv("~/LSMLH/1_dataSetSAPS2.csv",header = TRUE, sep=",")
#x = dataInput[, c(2:ncol(dataInput))]
#y = dataInput[, c(1)]
x = c(2:ncol(dataInput))
k = 5
neighbour = 1
detP = 0
nFeatures = ncol(dataInput)-1
while (detP <= 0){
  pMatrix = matrix(as.numeric(rexp(nFeatures^2, rate=.1),nFeatures))
  pMatrix = pMatrix %*% t(pMatrix)
  detP = det(pMatrix)
}

compactNeighborhood = matrix(0L, nrow = nrow(dataInput), ncol = k)
scatterNeighborhood = matrix(0L, nrow = nrow(dataInput), ncol = k)
for (i in 1:2){
  xi = dataInput[i,x]
  generalMahalanobis = matrix(data = NA, nrow = nrow(dataInput), ncol = 1)
  for (j in 1:2){
    if (i != j){
      xj = dataInput[j,x]
      generalMahalanobis[j,1] = t((xi - xj)) %*% pMatrix %*% (xi-xj)
      generalMahalanobis[j,2] = dataInput[i,1]
    }
  }
  order(generalMahalanobis$V1,decreasing = FALSE)
  while (neighbour < k*2){
    if (generalMahalanobis[neighbour,2] == dataInput[i,1]){
      compactNeighborhood[i] = compactNeighborhood[i] + generalMahalanobis[neighbour, 1]
      neighbour = neighbour + 1
    }else{
      scatterNeighborhood[i] = scatterNeighborhood[i] + generalMahalanobis[neighbour, 1]
      neighbour = neighbour + 1
    }
  }
}

objectiveFunction = sum(compactNeighborhood)/sum(scatterNeighborhood)

 


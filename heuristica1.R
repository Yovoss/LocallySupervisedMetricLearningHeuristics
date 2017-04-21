dataInput <- read.csv(file="~/Documentos/Metaheuristics/LocallySupervisedMetricLearningHeuristics/1_dataSetSAPS2.csv",header = TRUE, sep=",")
x = dataInput[, c(2:ncol(dataInput))]
y = dataInput[, c(1)]
k = 5
detP = 0
nFeatures = ncol(x)
while (detP <= 0){
  pMatrix = matrix(rexp(nFeatures^2, rate=.1),nFeatures)
  pMatrix = pMatrix %*% t(pMatrix)
  detP = det(pMatrix)
}


for (i in 1:nrow(x)){
  xi = x[1,]
  for (j in 1:nrow(x)){
    if (i != j){
      xj = x[2,]
      generalMahalanobis[i,j] = t((xi - xj)) %*% pMatrix %*% (xi-xj) 
    }
  }
  #compactNeighborhood = 
  #scatterNeighborhood =
}

 
ObjectiveFunction= t(input)

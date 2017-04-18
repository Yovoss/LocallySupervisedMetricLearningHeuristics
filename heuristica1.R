input <- read_delim("~/Documentos/MIMICIII/SAPSII.csv",";", escape_double = FALSE, trim_ws = TRUE)
k = 5
nFeatures = ncol(input)
pMatrix = matrix(rexp(nFeatures^2, rate=.1),nFeatures)

for (i in 1:nrow(input)){
  xi = input[1,]
  for (j in 1:nrow(input)){
    if (i != j){
      xj = input[2,]
      generalMahalanobis[i,j] = t(xi - xj) * pMatrix * (xi-xj) 
    }
  }
  #compactNeighborhood = 
  #scatterNeighborhood =
}

 
ObjectiveFunction= t(input)

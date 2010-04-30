source("sa-fns.R")
testSet <- as.matrix(read.table('mathematica-model.txt'))
xmodel <- testSet[,1]
training <- testSet[,2]

thetaTest <- c(0.007, 0.00383, 0.929)

cm <- makeCovMatrix(as.matrix(xmodel), thetaTest, length(xmodel))

hm <- makeHMatrix(as.matrix(xmodel), 19, 2)

bv <- makeBetaVector(hm, solve(cm), training)

ev <- makeEVector(bv, hm, training, solve(cm))

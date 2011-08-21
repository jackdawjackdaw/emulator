#dyn.load("~/Projects/Emulator/build/src/libRbind/libRBIND.dylib")
dyn.load("~/Projects/Emulator/build/src/libRbind/libRBIND.so")



nr <- 3
nc <- 4

mat <- matrix(0, nrow=nr, ncol=nc)
mat[,1] <- seq(1,3)
mat[,2] <- seq(4,6)
mat[,3] <- seq(7,9)

mat

.C("testConvert", as.double(mat), as.integer(nc), as.integer(nr))

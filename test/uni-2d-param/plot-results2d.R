xy.grid <- read.table("sample_locations.dat")

emu.mean <- matrix(as.matrix(read.table("emu_mean.dat")), 32, 32)
emu.var <- matrix(as.matrix(read.table("emu_var.dat")), 32, 32)


true.function <- function(z) {
  x <- z[1]
  y <- z[2]
  (x-0.5)**2 + (y-0.5)**2
}

true.values <- matrix(apply(xy.grid, 1, true.function), 32, 32)



plot.results <- function(){
  par(mfrow=c(2,2))
  image(emu.mean,sub="mean")
  contour(emu.mean, xlab="x", ylab="y",  add=TRUE)
  image(emu.var,sub="var")
  contour(emu.var, xlab="x", ylab="y",  add=TRUE)
  image(true.values,sub="true")
  contour(true.values, xlab="x", ylab="y", add=TRUE)

  implaus <- (true.values - emu.mean)^2 / emu.var
  image(implaus, sub="implaus")
  contour(implaus, xlab="x", ylab="y", add=TRUE)
}

plot.results()

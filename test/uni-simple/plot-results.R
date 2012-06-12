# read the output from sample-emulator.sh
samples <- read.table("./sample_locations.dat")
emu.mean <- read.table("emu_mean.dat")
emu.var <- read.table("emu_var.dat")

model <- read.table("input_model_file.dat", skip=3)

plot(samples$V1, emu.mean$V1, type="l", col="blue", ylim=c(-0.1,1.5), xlab="x", ylab="y")
points(t(model), pch=19)
lines(samples$V1, emu.mean$V1+1.956*sqrt(emu.var$V1), col="red")
lines(samples$V1, emu.mean$V1-1.956*sqrt(emu.var$V1), col="red")
legend("topright", c("mean", "confidence", "training"), lty=c(1,1, -1), col=c("blue", "red", "black"), pch=c(-1,-1,19))
title(main="simple 1d emulator test")

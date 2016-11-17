setwd("~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/Oct8_30mbp/Outputs_Oct")


dat <- read.csv("SummStats_deletsOnly_Oct8-N10000.csv")
bendat <- read.csv("SummStats_withBens_Oct8-N10000.csv")

cols.to.aggregate <- c("generation", "theta", "theta.neut", "pi", "pi_n", "pi_s", "pi_n.pi_s", "mean.delet.muts.per.ind.poly", "var.delet.muts.per.ind.poly", "mean.neut.muts.per.ind.poly", "var.neut.muts.per.ind.poly", "mean.total.muts.per.ind.poly", "var.total.muts.per.ind.poly", "mean.delet.muts.per.ind.all", "var.delet.muts.per.ind.all", "mean.neut.muts.per.ind.all", "var.neut.muts.per.ind.all", "mean.total.muts.per.ind.all", "var.total.muts.per.ind.all", "num.delet.muts.fixed", "num.neut.muts.fixed", "mean.fitness.poly", "var.fitness.poly", "mean.fitness.total", "var.fitness.total")
dat$split.ID <- matrix(unlist(strsplit(as.character(dat$file), split="rep")), ncol=2, byrow=TRUE)[,1]
dat$to.average <- paste(dat$split.ID, dat$generation, sep="")
avg.dat <- aggregate(dat[cols.to.aggregate], by=list(dat$to.average), FUN=mean, na.rm=TRUE)
avg.dat <- avg.dat[order(avg.dat$generation) ,]

cols.to.aggregate <- c("generation", "theta", "theta.neut", "pi", "pi_n", "pi_s", "pi_n.pi_s", "mean.delet.muts.per.ind.poly", "var.delet.muts.per.ind.poly", "mean.ben.muts.per.ind.poly", "var.ben.muts.per.ind.poly", "mean.neut.muts.per.ind.poly", "var.neut.muts.per.ind.poly", "mean.total.muts.per.ind.poly", "var.total.muts.per.ind.poly", "mean.delet.muts.per.ind.all", "var.delet.muts.per.ind.all", "mean.ben.muts.per.ind.all", "var.ben.muts.per.ind.all", "mean.neut.muts.per.ind.all", "var.neut.muts.per.ind.all", "mean.total.muts.per.ind.all", "var.total.muts.per.ind.all", "num.delet.muts.fixed", "num.ben.muts.fixed", "num.neut.muts.fixed", "mean.fitness.poly", "var.fitness.poly", "mean.fitness.total", "var.fitness.total")
bendat$split.ID <- matrix(unlist(strsplit(as.character(bendat$file), split="rep")), ncol=2, byrow=TRUE)[,1]
bendat$to.average <- paste(bendat$split.ID, bendat$generation, sep="")
avg.bendat <- aggregate(bendat[cols.to.aggregate], by=list(bendat$to.average), FUN=mean, na.rm=TRUE)
avg.bendat <- avg.bendat[order(avg.bendat$generation) ,]

outc.dat <- avg.dat[grep("outc", avg.dat$Group.1) ,]
self.dat99 <- avg.dat[grep("self0.99", avg.dat$Group.1) ,]
self.dat90 <- avg.dat[grep("self0.90", avg.dat$Group.1) ,]
asex.dat99 <- avg.dat[grep("asex0.99", avg.dat$Group.1) ,]
asex.dat90 <- avg.dat[grep("asex0.90", avg.dat$Group.1) ,]
bens.outc.dat <- avg.bendat[grep("outc", avg.bendat$Group.1) ,]
bens.self.dat99 <- avg.bendat[grep("self0.99", avg.bendat$Group.1) ,]
bens.self.dat90 <- avg.bendat[grep("self0.90", avg.bendat$Group.1) ,]
bens.asex.dat99 <- avg.bendat[grep("asex0.99", avg.bendat$Group.1) ,]
bens.asex.dat90 <- avg.bendat[grep("asex0.90", avg.bendat$Group.1) ,]

cols <- c("green3", "green3", "darkorchid2", "darkorchid2", "blue", "blue", "orange", "orange", "red", "red")
jitter <- 25

	
pdf("SummaryStatsResults_Oct8-N10000.pdf", width=10, height=9)
par(mar=c(4,4,1,0.5))	# mfrow=c(1,2), 
layout(matrix(c(1,2,3,3), nrow=2, ncol=2))

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 0.00025), xlab="Generation", ylab="Theta_W")
points(outc.dat$generation, outc.dat$theta, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$theta, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$theta, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$theta, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$theta, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$theta, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$theta, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$theta, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$theta, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$theta, type="o", pch=16, col="darkorchid2")
legend("topright", pch=rep(c(2,16), 5), lty=rep(c(2,1), 5), col=c("green3", "green3", "darkorchid2", "darkorchid2", "blue", "blue", "orange", "orange", "red", "red"), c("Outc del", "Outc ben-del", "90% asex del", "90% asex ben-del", "99% asex del", "99% asex ben-del", "90% self del", "90% self ben-del", "99% self del", "99% self ben-del"), bg="white", ncol=2, cex=0.95)


plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 3500), xlab="Generation", ylab="Mean number polymorphic muts per ind", main="Polymorphic mutations")
points(outc.dat$generation, outc.dat$mean.delet.muts.per.ind.poly, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$mean.delet.muts.per.ind.poly, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$mean.delet.muts.per.ind.poly, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$mean.delet.muts.per.ind.poly, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$mean.delet.muts.per.ind.poly, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$mean.delet.muts.per.ind.poly, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$mean.delet.muts.per.ind.poly, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$mean.delet.muts.per.ind.poly, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$mean.delet.muts.per.ind.poly, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$mean.delet.muts.per.ind.poly, type="o", pch=16, col="darkorchid2")

points(bens.outc.dat$generation + jitter, bens.outc.dat$mean.ben.muts.per.ind.poly, type="o", pch=2, col="green3")
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$mean.ben.muts.per.ind.poly, type="o", pch=2, col="red")
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$mean.ben.muts.per.ind.poly, type="o", pch=2, col="orange")
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$mean.ben.muts.per.ind.poly, type="o", pch=2, col="blue")
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$mean.ben.muts.per.ind.poly, type="o", pch=2, col="darkorchid2")
legend("topright", pch=c(16,2), col=c("black", "black"), c("Deleterious", "Beneficial"), cex=1, bg="white", ncol=1)



#_____________________________________________________________________________________________________________________#
#	fitness plots

plot(0, type="n", xlim=c(10000,100000), ylim=c(0,3), xlab="Generation", ylab="Mean fitness")
points(outc.dat$generation, outc.dat$mean.fitness.poly, type="o", col="green3", lty=2, pch=18)
points(bens.outc.dat$generation + jitter, bens.outc.dat$mean.fitness.poly, type="o", pch=18, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$mean.fitness.poly, type="o", col="red", lty=2, pch=18)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$mean.fitness.poly, type="o", pch=18, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$mean.fitness.poly, type="o", col="orange", lty=2, pch=18)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$mean.fitness.poly, type="o", pch=18, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$mean.fitness.poly, type="o", col="blue", lty=2, pch=18)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$mean.fitness.poly, type="o", pch=18, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$mean.fitness.poly, type="o", col="darkorchid2", lty=2, pch=18)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$mean.fitness.poly, type="o", pch=18, col="darkorchid2")

points(outc.dat$generation, outc.dat$mean.fitness.total, type="o", pch=5, col="green3", lty=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$mean.fitness.total, type="o", pch=5, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$mean.fitness.total, type="o", pch=5, col="red", lty=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$mean.fitness.total, type="o", pch=5, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$mean.fitness.total, type="o", pch=5, col="orange", lty=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$mean.fitness.total, type="o", pch=5, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$mean.fitness.total, type="o", pch=5, col="blue", lty=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$mean.fitness.total, type="o", pch=5, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$mean.fitness.total, type="o", pch=5, col="darkorchid2", lty=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$mean.fitness.total, type="o", pch=5, col="darkorchid2")
legend("topleft", c("Polymorphic loci", "Total"), pch=c(18,5), col=c("black"), cex=1, bg="white", ncol=1)
legend("topright", pch="", lty=rep(c(2,1),5), lwd=2, col=cols, c("Outc del", "Outc ben-del", "90% asex del", "90% asex ben-del", "99% asex del", "99% asex ben-del", "90% self del", "90% self ben-del", "99% self del", "99% self ben-del"), bg="white", ncol=2, cex=0.95)

#_____________________________________________________________________________________________________________________#
#	pi plots

par(mfrow=c(2,2), mar=c(4,4,1,0.5))
plot(0, type="n", xlim=c(10000,100000), ylim=c(0.25, 0.9), xlab="Generation", ylab="pi_n / pi_s", main="pi_n / pi_s")
points(outc.dat$generation, outc.dat$pi_n.pi_s, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$pi_n.pi_s, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$pi_n.pi_s, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$pi_n.pi_s, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$pi_n.pi_s, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$pi_n.pi_s, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$pi_n.pi_s, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$pi_n.pi_s, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$pi_n.pi_s, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$pi_n.pi_s, type="o", pch=16, col="darkorchid2")

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 0.001), xlab="Generation", ylab="pi", main="pi")
points(outc.dat$generation, outc.dat$pi, type="o", pch=2, col="green3", lty=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$pi, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$pi, type="o", pch=2, col="red", lty=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$pi, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$pi, type="o", pch=2, col="orange", lty=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$pi, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$pi, type="o", pch=2, col="blue", lty=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$pi, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$pi, type="o", pch=2, col="darkorchid2", lty=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$pi, type="o", pch=16, col="darkorchid2")
legend("topright", pch=rep(c(2,16), 5), lty=rep(c(2,1), 5), col=c("green3", "green3", "darkorchid2", "darkorchid2", "blue", "blue", "orange", "orange", "red", "red"), c("Outc del", "Outc ben-del", "90% asex del", "90% asex ben-del", "99% asex del", "99% asex ben-del", "90% self del", "90% self ben-del", "99% self del", "99% self ben-del"), bg="white", ncol=2, cex=0.95)

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 0.001), xlab="Generation", ylab="pi_n", main="pi_n")
points(outc.dat$generation, outc.dat$pi_n, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$pi_n, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$pi_n, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$pi_n, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$pi_n, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$pi_n, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$pi_n, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$pi_n, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$pi_n, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$pi_n, type="o", pch=16, col="darkorchid2")

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 0.001), xlab="Generation", ylab="pi_s", main="pi_s")
points(outc.dat$generation, outc.dat$pi_s, type="o", pch=2, col="green3", lty=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$pi_s, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$pi_s, type="o", pch=2, col="red", lty=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$pi_s, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$pi_s, type="o", pch=2, col="orange", lty=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$pi_s, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$pi_s, type="o", pch=2, col="blue", lty=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$pi_s, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$pi_s, type="o", pch=2, col="darkorchid2", lty=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$pi_s, type="o", pch=16, col="darkorchid2")


#_____________________________________________________________________________________________________________________#
#	mutation plots

# Ka/Ks =  ratio of the number of nonsynonymous substitutions per non-synonymous site (Ka), in a given period of time, to the number of synonymous substitutions per synonymous site (Ks), in the same period. The latter are assumed to be neutral, so that the ratio indicates the net balance between deleterious and beneficial mutations. Values of Ka/Ks significantly above 1 are unlikely to occur without at least some of the mutations being advantageous. If beneficial mutations are assumed to make little contribution, then Ks estimates the degree of evolutionary constraint.
# The ratio is also known as ω or dN/dS

geno.size <- 30000000	# 30 Mbp

#dN/dS total
out.bendel.dnds <- ((bens.outc.dat$num.delet.muts.fixed + bens.outc.dat$num.ben.muts.fixed)/(geno.size*0.75)) / (bens.outc.dat$num.neut.muts.fixed/(geno.size*0.25))
out.del.dnds <- (outc.dat$num.delet.muts.fixed/(geno.size*0.75)) / (outc.dat$num.neut.muts.fixed/(geno.size*0.25))
self99.bendel.dnds <- ((bens.self.dat99$num.delet.muts.fixed + bens.self.dat99$num.ben.muts.fixed)/(geno.size*0.75)) / (bens.self.dat99$num.neut.muts.fixed/(geno.size*0.25))
self99.del.dnds <- (self.dat99$num.delet.muts.fixed/(geno.size*0.75)) / (self.dat99$num.neut.muts.fixed/(geno.size*0.25))
self90.bendel.dnds <- ((bens.self.dat90$num.delet.muts.fixed + bens.self.dat90$num.ben.muts.fixed)/(geno.size*0.75)) / (bens.self.dat90$num.neut.muts.fixed/(geno.size*0.25))
self90.del.dnds <- (self.dat90$num.delet.muts.fixed/(geno.size*0.75)) / (self.dat90$num.neut.muts.fixed/(geno.size*0.25))
asex99.bendel.dnds <- ((bens.asex.dat99$num.delet.muts.fixed + bens.asex.dat99$num.ben.muts.fixed)/(geno.size*0.75)) / (bens.asex.dat99$num.neut.muts.fixed/(geno.size*0.25))
asex99.del.dnds <- (asex.dat99$num.delet.muts.fixed/(geno.size*0.75)) / (asex.dat99$num.neut.muts.fixed/(geno.size*0.25))
asex90.bendel.dnds <- ((bens.asex.dat90$num.delet.muts.fixed + bens.asex.dat90$num.ben.muts.fixed)/(geno.size*0.75)) / (bens.asex.dat90$num.neut.muts.fixed/(geno.size*0.25))
asex90.del.dnds <- (asex.dat90$num.delet.muts.fixed/(geno.size*0.75)) / (asex.dat90$num.neut.muts.fixed/(geno.size*0.25))

#dN/dS deleterious only (can be calculated in either the presence or absence of beneficial muts)
out.bendel.dnds.delet <- (bens.outc.dat$num.delet.muts.fixed/(geno.size*0.75)) / (bens.outc.dat$num.neut.muts.fixed/(geno.size*0.25))
out.del.dnds.delet <- (outc.dat$num.delet.muts.fixed/(geno.size*0.75)) / (outc.dat$num.neut.muts.fixed/(geno.size*0.25))
self99.bendel.dnds.delet <- (bens.self.dat99$num.delet.muts.fixed/(geno.size*0.75)) / (bens.self.dat99$num.neut.muts.fixed/(geno.size*0.25))
self99.del.dnds.delet <- (self.dat99$num.delet.muts.fixed/(geno.size*0.75)) / (self.dat99$num.neut.muts.fixed/(geno.size*0.25))
self90.bendel.dnds.delet <- (bens.self.dat90$num.delet.muts.fixed/(geno.size*0.75)) / (bens.self.dat90$num.neut.muts.fixed/(geno.size*0.25))
self90.del.dnds.delet <- (self.dat90$num.delet.muts.fixed/(geno.size*0.75)) / (self.dat90$num.neut.muts.fixed/(geno.size*0.25))
asex99.bendel.dnds.delet <- (bens.asex.dat99$num.delet.muts.fixed/(geno.size*0.75)) / (bens.asex.dat99$num.neut.muts.fixed/(geno.size*0.25))
asex99.del.dnds.delet <- (asex.dat99$num.delet.muts.fixed/(geno.size*0.75)) / (asex.dat99$num.neut.muts.fixed/(geno.size*0.25))
asex90.bendel.dnds.delet <- (bens.asex.dat90$num.delet.muts.fixed/(geno.size*0.75)) / (bens.asex.dat90$num.neut.muts.fixed/(geno.size*0.25))
asex90.del.dnds.delet <- (asex.dat90$num.delet.muts.fixed/(geno.size*0.75)) / (asex.dat90$num.neut.muts.fixed/(geno.size*0.25))

#dN/dS beneficial only (can be calculated only in the presence of beneficial muts)
out.bendel.dnds.ben <- (bens.outc.dat$num.ben.muts.fixed/(geno.size*0.75)) / (bens.outc.dat$num.neut.muts.fixed/(geno.size*0.25))
self99.bendel.dnds.ben <- (bens.self.dat99$num.ben.muts.fixed/(geno.size*0.75)) / (bens.self.dat99$num.neut.muts.fixed/(geno.size*0.25))
self90.bendel.dnds.ben <- (bens.self.dat90$num.ben.muts.fixed/(geno.size*0.75)) / (bens.self.dat90$num.neut.muts.fixed/(geno.size*0.25))
asex99.bendel.dnds.ben <- (bens.asex.dat99$num.ben.muts.fixed/(geno.size*0.75)) / (bens.asex.dat99$num.neut.muts.fixed/(geno.size*0.25))
asex90.bendel.dnds.ben <- (bens.asex.dat90$num.ben.muts.fixed/(geno.size*0.75)) / (bens.asex.dat90$num.neut.muts.fixed/(geno.size*0.25))

par(mfrow=c(2,2), mar=c(4,4,1,0.5))

# dN/dS total
plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 1.2), xlab="Generation", ylab="dN / dS", main="dN / dS")
points(outc.dat$generation, out.del.dnds, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, out.bendel.dnds, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self99.del.dnds, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), self99.bendel.dnds, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self90.del.dnds, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), self90.bendel.dnds, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex99.del.dnds, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), asex99.bendel.dnds, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex90.del.dnds, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), asex90.bendel.dnds, type="o", pch=16, col="darkorchid2")
legend("topright", pch=rep(c(2,16), 5), lty=rep(c(2,1), 5), col=c("green3", "green3", "darkorchid2", "darkorchid2", "blue", "blue", "orange", "orange", "red", "red"), c("Outc del", "Outc ben-del", "90% asex del", "90% asex ben-del", "99% asex del", "99% asex ben-del", "90% self del", "90% self ben-del", "99% self del", "99% self ben-del"), bg="white", ncol=2, cex=0.95)

# ALPHA
plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 1), xlab="Generation", ylab="alpha", main="proportion adaptive substitutions")
points(outc.dat$generation, 0/out.del.dnds, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, out.bendel.dnds.ben / out.bendel.dnds, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), 0/self99.del.dnds, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), self99.bendel.dnds.ben / self99.bendel.dnds, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), 0/self90.del.dnds, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), self90.bendel.dnds.ben / self90.bendel.dnds, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), 0/asex99.del.dnds, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), asex99.bendel.dnds.ben / asex99.bendel.dnds, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), 0/asex90.del.dnds, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), asex90.bendel.dnds.ben / asex90.bendel.dnds, type="o", pch=16, col="darkorchid2")

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 1.2), xlab="Generation", ylab="dN / dS", main="dN/dS - Deleterious only")
points(outc.dat$generation, out.del.dnds.delet, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, out.bendel.dnds.delet, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self99.del.dnds.delet, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), self99.bendel.dnds.delet, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self90.del.dnds.delet, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), self90.bendel.dnds.delet, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex99.del.dnds.delet, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), asex99.bendel.dnds.delet, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex90.del.dnds.delet, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), asex90.bendel.dnds.delet, type="o", pch=16, col="darkorchid2")

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 1.2), xlab="Generation", ylab="dN / dS", main="dN / dS - Beneficial only")
points(bens.outc.dat$generation + jitter, out.bendel.dnds.ben, type="o", pch=16, col="green3")
points(bens.self.dat99$generation + (jitter*3), self99.bendel.dnds.ben, type="o", pch=16, col="red")
points(bens.self.dat90$generation + (jitter*5), self90.bendel.dnds.ben, type="o", pch=16, col="orange")
points(bens.asex.dat99$generation + (jitter*7), asex99.bendel.dnds.ben, type="o", pch=16, col="blue")
points(bens.asex.dat90$generation + (jitter*9), asex90.bendel.dnds.ben, type="o", pch=16, col="darkorchid2")


par(mfrow=c(2,2), mar=c(4,4,1,0.5))
plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 7500), xlab="Generation", ylab="Mean # (poly + fixed) mutations per individual", main="Deleterious mutations only")
points(outc.dat$generation, outc.dat$mean.delet.muts.per.ind.all, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$mean.delet.muts.per.ind.all, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$mean.delet.muts.per.ind.all, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$mean.delet.muts.per.ind.all, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$mean.delet.muts.per.ind.all, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$mean.delet.muts.per.ind.all, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$mean.delet.muts.per.ind.all, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$mean.delet.muts.per.ind.all, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$mean.delet.muts.per.ind.all, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$mean.delet.muts.per.ind.all, type="o", pch=16, col="darkorchid2")

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 1000), xlab="Generation", ylab="Mean # (poly + fixed) mutations per individual", main="Beneficial mutations only")
points(bens.outc.dat$generation + jitter, bens.outc.dat$mean.ben.muts.per.ind.all, type="o", pch=16, col="green3")
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$mean.ben.muts.per.ind.all, type="o", pch=16, col="red")
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$mean.ben.muts.per.ind.all, type="o", pch=16, col="orange")
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$mean.ben.muts.per.ind.all, type="o", pch=16, col="blue")
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$mean.ben.muts.per.ind.all, type="o", pch=16, col="darkorchid2")

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 5000), xlab="Generation", ylab="Total # Delet Fixed Mutations", main="Deleterious mutations only")
points(outc.dat$generation, outc.dat$num.delet.muts.fixed, type="o", col="green3", lty=2, pch=2)
points(bens.outc.dat$generation + jitter, bens.outc.dat$num.delet.muts.fixed, type="o", pch=16, col="green3")
points(self.dat99$generation + (jitter*2), self.dat99$num.delet.muts.fixed, type="o", col="red", lty=2, pch=2)
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$num.delet.muts.fixed, type="o", pch=16, col="red")
points(self.dat90$generation + (jitter*4), self.dat90$num.delet.muts.fixed, type="o", col="orange", lty=2, pch=2)
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$num.delet.muts.fixed, type="o", pch=16, col="orange")
points(asex.dat99$generation + (jitter*6), asex.dat99$num.delet.muts.fixed, type="o", col="blue", lty=2, pch=2)
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$num.delet.muts.fixed, type="o", pch=16, col="blue")
points(asex.dat90$generation + (jitter*8), asex.dat90$num.delet.muts.fixed, type="o", col="darkorchid2", lty=2, pch=2)
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$num.delet.muts.fixed, type="o", pch=16, col="darkorchid2")

plot(0, type="n", xlim=c(10000,100000), ylim=c(0, 1000), xlab="Generation", ylab="Total # Beneficial Fixed Mutations", main="Beneficial mutations only")
points(bens.outc.dat$generation + jitter, bens.outc.dat$num.ben.muts.fixed, type="o", pch=16, col="green3")
points(bens.self.dat99$generation + (jitter*3), bens.self.dat99$num.ben.muts.fixed, type="o", pch=16, col="red")
points(bens.self.dat90$generation + (jitter*5), bens.self.dat90$num.ben.muts.fixed, type="o", pch=16, col="orange")
points(bens.asex.dat99$generation + (jitter*7), bens.asex.dat99$num.ben.muts.fixed, type="o", pch=16, col="blue")
points(bens.asex.dat90$generation + (jitter*9), bens.asex.dat90$num.ben.muts.fixed, type="o", pch=16, col="darkorchid2")
legend("topleft", pch=rep(c(2,16), 5), lty=rep(c(2,1), 5), col=c("green3", "green3", "darkorchid2", "darkorchid2", "blue", "blue", "orange", "orange", "red", "red"), c("Outc del", "Outc ben-del", "90% asex del", "90% asex ben-del", "99% asex del", "99% asex ben-del", "90% self del", "90% self ben-del", "99% self del", "99% self ben-del"), bg="white", ncol=2, cex=0.95)


dev.off()


setwd("~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/Nov3_TransitionMatingSystems/Outputs_Nov/RelativeLoadCalcs_TransitionMatingSystem")


## LOAD PLOTS

# where load is the relative fixation of deleterious mutations in the transitioned mating systems to the outcrossers


# Nov 10, transition mating systems:




# can make these raw count plots from the summ stats already calculated:
dat <- read.csv("../SummStats_deletsOnly_transMateSys_Nov10-N10000.csv")
bendat <- read.csv("../SummStats_withBens_transMateSys_Nov10-N10000.csv")
reldat <- read.csv("RelativeLoadCalculated_Nov10_TransMateSys.csv", row.names=NULL)

# because these are in subfolders by rep
temp.filenames <- matrix(unlist(strsplit(as.character(reldat$filename), split="/")), ncol=2, byrow=TRUE)[,2]
reldat$split.ID <- matrix(unlist(strsplit(temp.filenames, split="rep")), ncol=2, byrow=TRUE)[,1]



dat2 <- read.csv("../SummStats_deletsOnly_transMateSys_Nov17-N10000.csv")
bendat2 <- read.csv("../SummStats_withBens_transMateSys_Nov17-N10000.csv")
reldat2 <- read.csv("RelativeLoadCalculated_Nov17_TransMateSys.csv", row.names=NULL)

reldat2$split.ID <- matrix(unlist(strsplit(as.character(reldat2$filename), split="rep")), ncol=2, byrow=TRUE)[,1]

## split again to get rid of the 2 dates (Nov 10 and Nov 17)
reldat$split.ID <- matrix(unlist(strsplit(reldat$split.ID, split="_N10000_")), ncol=2, byrow=TRUE)[,2]


dat <- rbind(dat, dat2)
bendat <- rbind(bendat, bendat2)
reldat <- rbind(reldat, reldat2)



library(Hmisc)

	# also do 95% CIs:
	avg.bendat.ci95 <- aggregate(bendat[cols.to.aggregate], by=list(bendat$to.average), FUN=sd, na.rm=TRUE)
	avg.bendat.ci95 <- (avg.bendat.ci95[cols.to.aggregate]/sqrt(10)) * 2.26
	avg.bendat.ci95 <- cbind(avg.bendat$Group.1, avg.bendat.ci95)		# put the file names back in there
	names(avg.bendat.ci95) <- c("Group.1", "generation", "theta", "theta.neut", "pi", "pi_n", "pi_s", "pi_n.pi_s", "mean.delet.muts.per.ind.poly", "var.delet.muts.per.ind.poly", "mean.ben.muts.per.ind.poly", "var.ben.muts.per.ind.poly", "mean.neut.muts.per.ind.poly", "var.neut.muts.per.ind.poly", "mean.total.muts.per.ind.poly", "var.total.muts.per.ind.poly", "mean.delet.muts.per.ind.all", "var.delet.muts.per.ind.all", "mean.ben.muts.per.ind.all", "var.ben.muts.per.ind.all", "mean.neut.muts.per.ind.all", "var.neut.muts.per.ind.all", "mean.total.muts.per.ind.all", "var.total.muts.per.ind.all", "num.delet.muts.fixed", "num.ben.muts.fixed", "num.neut.muts.fixed", "mean.fitness.poly", "var.fitness.poly", "mean.fitness.total", "var.fitness.total")

	errbar(add=TRUE, asex.dat99$generation + (jitter*6), asex.dat99$theta.neut, yminus=asex.dat99$theta.neut-ci.asex.dat99$theta.neut, yplus=asex.dat99$theta.neut+ci.asex.dat99$theta.neut, errbar.col="blue", col="white", cex=0.1)



plot.load <- function(dat, bendat, reldat, gens.to.sample.at, pop.size, figure.basename, xlimits){
	
	
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
	
	## also want to make plots that incorporate s values of the fixed mutations
	# plot 1 - (fitness for delets mating system x) / (fitness for delets outcrosser)
		
	
	# first average everything per group
	cols.to.aggregate <- c(names(reldat)[-c(1:2,length(names(reldat)))])	
	
	reldat$to.average <- paste(reldat$split.ID, reldat$generation, sep="")
	avg.reldat <- aggregate(reldat[cols.to.aggregate], by=list(reldat$to.average), FUN=mean, na.rm=TRUE)
	avg.reldat <- avg.reldat[order(avg.reldat$generation) ,]
	
	# sections to parse out:
	#	total.s.window.count
	#	fixed.s.window.count
	#	poly.s.window.count
	#	s.window.fitness.total
	#	s.window.fitness.poly

	# FITNESS ONLY FOR DELETS IN WINDOWS
	avg.total.fitness <- cbind(avg.reldat[,c(1:2)], avg.reldat[, grep("s.window.delet.fitness.total", names(avg.reldat))])
	# rename column headers because these are now all total fitness columns

	# combine into the 4 bins we want for now:		IS FITNESS ADDITIVE HERE?????????? no, for now I'm doing it multiplicatively
	avg.total.fitness$s.window.delet.fitness.total_1_0.01 <- avg.total.fitness$s.window.delet.fitness.total_.1_.0.1 * avg.total.fitness$s.window.delet.fitness.total_.0.1_.0.01
	avg.total.fitness$s.window.delet.fitness.total_.1e.04_0 <- avg.total.fitness$s.window.delet.fitness.total_.1e.04_.1e.05 * avg.total.fitness$s.window.delet.fitness.total_.1e.05_0
	avg.total.fitness$s.window.delet.fitness.total_.1_0 <- avg.total.fitness$s.window.delet.fitness.total_.1_.0.1 * avg.total.fitness$s.window.delet.fitness.total_.0.1_.0.01 * avg.total.fitness$s.window.delet.fitness.total_.0.01_.0.001 * avg.total.fitness$s.window.delet.fitness.total_.0.001_.1e.04 * avg.total.fitness$s.window.delet.fitness.total_.1e.04_.1e.05 * avg.total.fitness$s.window.delet.fitness.total_.1e.05_0
	
	outc.bendel.delet.fitdat <- avg.total.fitness[grep("_ben-del_TransTooutc", avg.total.fitness$Group.1) ,]
	outc.del.delet.fitdat <- avg.total.fitness[grep("_del_TransTooutc", avg.total.fitness$Group.1) ,]
	self.bendel.delet.fitdat99 <- avg.total.fitness[grep("_ben-del_TransToself0.99", avg.total.fitness$Group.1) ,]
	self.del.delet.fitdat99 <- avg.total.fitness[grep("_del_TransToself0.99", avg.total.fitness$Group.1) ,]
	self.bendel.delet.fitdat90 <- avg.total.fitness[grep("_ben-del_TransToself0.90", avg.total.fitness$Group.1) ,]
	self.del.delet.fitdat90 <- avg.total.fitness[grep("_del_TransToself0.90", avg.total.fitness$Group.1) ,]
	asex.bendel.delet.fitdat99 <- avg.total.fitness[grep("_ben-del_TransToasex0.99", avg.total.fitness$Group.1) ,]
	asex.del.delet.fitdat99 <- avg.total.fitness[grep("_del_TransToasex0.99", avg.total.fitness$Group.1) ,]
	asex.bendel.delet.fitdat90 <- avg.total.fitness[grep("_ben-del_TransToasex0.90", avg.total.fitness$Group.1) ,]
	asex.del.delet.fitdat90 <- avg.total.fitness[grep("_del_TransToasex0.90", avg.total.fitness$Group.1) ,]
	
	# FULL FITNESS FOR DELETS IN WINDOWS AND ALL BENEFICIAL EFFECTS INCLUDED
	avg.bendel.total.fitness <- cbind(avg.reldat[,c(1:2)], avg.reldat[, grep("s.window.bendel.fitness.total", names(avg.reldat))])
	# rename column headers because these are now all total fitness columns
	
	# combine into the 4 bins we want for now:		IS FITNESS ADDITIVE HERE?????????? no, for now I'm doing it multiplicatively
	avg.bendel.total.fitness$s.window.bendel.fitness.total_1_0.01 <- avg.bendel.total.fitness$s.window.bendel.fitness.total_.1_.0.1 * avg.bendel.total.fitness$s.window.bendel.fitness.total_.0.1_.0.01
	avg.bendel.total.fitness$s.window.bendel.fitness.total_.1e.04_0 <- avg.bendel.total.fitness$s.window.bendel.fitness.total_.1e.04_.1e.05 * avg.bendel.total.fitness$s.window.bendel.fitness.total_.1e.05_0
	avg.bendel.total.fitness$s.window.bendel.fitness.total_.1_0 <- avg.bendel.total.fitness$s.window.bendel.fitness.total_.1_.0.1 * avg.bendel.total.fitness$s.window.bendel.fitness.total_.0.1_.0.01 * avg.bendel.total.fitness$s.window.bendel.fitness.total_.0.01_.0.001 * avg.bendel.total.fitness$s.window.bendel.fitness.total_.0.001_.1e.04 * avg.bendel.total.fitness$s.window.bendel.fitness.total_.1e.04_.1e.05 * avg.bendel.total.fitness$s.window.bendel.fitness.total_.1e.05_0
	
	outc.bendel.bendel.fitdat <- avg.bendel.total.fitness[grep("_ben-del_TransTooutc", avg.bendel.total.fitness$Group.1) ,]
	outc.del.bendel.fitdat <- avg.bendel.total.fitness[grep("_del_TransTooutc", avg.bendel.total.fitness$Group.1) ,]
	self.bendel.bendel.fitdat99 <- avg.bendel.total.fitness[grep("_ben-del_TransToself0.99", avg.bendel.total.fitness$Group.1) ,]
	self.del.bendel.fitdat99 <- avg.bendel.total.fitness[grep("_del_TransToself0.99", avg.bendel.total.fitness$Group.1) ,]
	self.bendel.bendel.fitdat90 <- avg.bendel.total.fitness[grep("_ben-del_TransToself0.90", avg.bendel.total.fitness$Group.1) ,]
	self.del.bendel.fitdat90 <- avg.bendel.total.fitness[grep("_del_TransToself0.90", avg.bendel.total.fitness$Group.1) ,]
	asex.bendel.bendel.fitdat99 <- avg.bendel.total.fitness[grep("_ben-del_TransToasex0.99", avg.bendel.total.fitness$Group.1) ,]
	asex.del.bendel.fitdat99 <- avg.bendel.total.fitness[grep("_del_TransToasex0.99", avg.bendel.total.fitness$Group.1) ,]
	asex.bendel.bendel.fitdat90 <- avg.bendel.total.fitness[grep("_ben-del_TransToasex0.90", avg.bendel.total.fitness$Group.1) ,]
	asex.del.bendel.fitdat90 <- avg.bendel.total.fitness[grep("_del_TransToasex0.90", avg.bendel.total.fitness$Group.1) ,]
	
	avg.reldat[is.na(avg.reldat)] <- 0

	avg.total.muts <- cbind(avg.reldat[,c(1:2)], avg.reldat[, grep("s.window.total.muts", names(avg.reldat))])
	avg.poly.muts <- cbind(avg.reldat[,c(1:2)], avg.reldat[, grep("s.window.poly.muts", names(avg.reldat))])
	avg.fixed.muts <- cbind(avg.reldat[,c(1:2)], avg.reldat[, grep("s.window.fixed.muts", names(avg.reldat))])
	# rename column headers because these are now all total fitness columns
	
	# sum to get all, not just in bins of S
	avg.total.muts$s.window.total.muts_.1_0 <- avg.total.muts$s.window.total.muts_.1_.0.1 + avg.total.muts$s.window.total.muts_.0.1_.0.01 + avg.total.muts$s.window.total.muts_.0.01_.0.001 + avg.total.muts$s.window.total.muts_.0.001_.1e.04 + avg.total.muts$s.window.total.muts_.1e.04_.1e.05 + avg.total.muts$s.window.total.muts_.1e.05_0
	avg.poly.muts$s.window.poly.muts_.1_0 <- avg.poly.muts$s.window.poly.muts_.1_.0.1 + avg.poly.muts$s.window.poly.muts_.0.1_.0.01 + avg.poly.muts$s.window.poly.muts_.0.01_.0.001 + avg.poly.muts$s.window.poly.muts_.0.001_.1e.04 + avg.poly.muts$s.window.poly.muts_.1e.04_.1e.05 + avg.poly.muts$s.window.poly.muts_.1e.05_0
	avg.fixed.muts$s.window.fixed.muts_.1_0 <- avg.fixed.muts$s.window.fixed.muts_.1_.0.1 + avg.fixed.muts$s.window.fixed.muts_.0.1_.0.01 + avg.fixed.muts$s.window.fixed.muts_.0.01_.0.001 + avg.fixed.muts$s.window.fixed.muts_.0.001_.1e.04 + avg.fixed.muts$s.window.fixed.muts_.1e.04_.1e.05 + avg.fixed.muts$s.window.fixed.muts_.1e.05_0

	# total number DELET muts (fixed plus poly)
	outc.bendel.totmuts <- avg.total.muts[grep("_ben-del_TransTooutc", avg.total.muts$Group.1) ,]
	outc.del.totmuts <- avg.total.muts[grep("_del_TransTooutc", avg.total.muts$Group.1) ,]
	self.bendel.totmuts99 <- avg.total.muts[grep("_ben-del_TransToself0.99", avg.total.muts$Group.1) ,]
	self.del.totmuts99 <- avg.total.muts[grep("_del_TransToself0.99", avg.total.muts$Group.1) ,]
	self.bendel.totmuts90 <- avg.total.muts[grep("_ben-del_TransToself0.90", avg.total.muts$Group.1) ,]
	self.del.totmuts90 <- avg.total.muts[grep("_del_TransToself0.90", avg.total.muts$Group.1) ,]
	asex.bendel.totmuts99 <- avg.total.muts[grep("_ben-del_TransToasex0.99", avg.total.muts$Group.1) ,]
	asex.del.totmuts99 <- avg.total.muts[grep("_del_TransToasex0.99", avg.total.muts$Group.1) ,]
	asex.bendel.totmuts90 <- avg.total.muts[grep("_ben-del_TransToasex0.90", avg.total.muts$Group.1) ,]
	asex.del.totmuts90 <- avg.total.muts[grep("_del_TransToasex0.90", avg.total.muts$Group.1) ,]
	
	# number POLY DELET muts
	outc.bendel.polymuts <- avg.poly.muts[grep("_ben-del_TransTooutc", avg.poly.muts$Group.1) ,]
	outc.del.polymuts <- avg.poly.muts[grep("_del_TransTooutc", avg.poly.muts$Group.1) ,]
	self.bendel.polymuts99 <- avg.poly.muts[grep("_ben-del_TransToself0.99", avg.poly.muts$Group.1) ,]
	self.del.polymuts99 <- avg.poly.muts[grep("_del_TransToself0.99", avg.poly.muts$Group.1) ,]
	self.bendel.polymuts90 <- avg.poly.muts[grep("_ben-del_TransToself0.90", avg.poly.muts$Group.1) ,]
	self.del.polymuts90 <- avg.poly.muts[grep("_del_TransToself0.90", avg.poly.muts$Group.1) ,]
	asex.bendel.polymuts99 <- avg.poly.muts[grep("_ben-del_TransToasex0.99", avg.poly.muts$Group.1) ,]
	asex.del.polymuts99 <- avg.poly.muts[grep("_del_TransToasex0.99", avg.poly.muts$Group.1) ,]
	asex.bendel.polymuts90 <- avg.poly.muts[grep("_ben-del_TransToasex0.90", avg.poly.muts$Group.1) ,]
	asex.del.polymuts90 <- avg.poly.muts[grep("_del_TransToasex0.90", avg.poly.muts$Group.1) ,]
	
	# number FIXED DELET muts
	outc.bendel.fixedmuts <- avg.fixed.muts[grep("_ben-del_TransTooutc", avg.fixed.muts$Group.1) ,]
	outc.del.fixedmuts <- avg.fixed.muts[grep("_del_TransTooutc", avg.fixed.muts$Group.1) ,]
	self.bendel.fixedmuts99 <- avg.fixed.muts[grep("_ben-del_TransToself0.99", avg.fixed.muts$Group.1) ,]
	self.del.fixedmuts99 <- avg.fixed.muts[grep("_del_TransToself0.99", avg.fixed.muts$Group.1) ,]
	self.bendel.fixedmuts90 <- avg.fixed.muts[grep("_ben-del_TransToself0.90", avg.fixed.muts$Group.1) ,]
	self.del.fixedmuts90 <- avg.fixed.muts[grep("_del_TransToself0.90", avg.fixed.muts$Group.1) ,]
	asex.bendel.fixedmuts99 <- avg.fixed.muts[grep("_ben-del_TransToasex0.99", avg.fixed.muts$Group.1) ,]
	asex.del.fixedmuts99 <- avg.fixed.muts[grep("_del_TransToasex0.99", avg.fixed.muts$Group.1) ,]
	asex.bendel.fixedmuts90 <- avg.fixed.muts[grep("_ben-del_TransToasex0.90", avg.fixed.muts$Group.1) ,]
	asex.del.fixedmuts90 <- avg.fixed.muts[grep("_del_TransToasex0.90", avg.fixed.muts$Group.1) ,]
	
	#___________________________________________________________________________________________________________________________________________________________#
	#
	#		DATA ALL READ IN
	#___________________________________________________________________________________________________________________________________________________________#
		
	pdf(paste(c("RelativeandRawLoad-", figure.basename, ".pdf"), collapse=""), width=8, height=16)
	par(mar=c(4,4,1,1), mfrow=c(5,2))
	
	# FITNESS JUST FOR DELS
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 1), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative (to outcrossers) load for delet muts", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.delet.fitdat99$generation, 1-(self.bendel.delet.fitdat99$s.window.delet.fitness.total_.1_0 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1_0), col="red", type="o", pch=16)
	points(self.del.delet.fitdat99$generation, 1-(self.del.delet.fitdat99$s.window.delet.fitness.total_.1_0 / outc.del.delet.fitdat$s.window.delet.fitness.total_.1_0), col="red", type="o", pch=2, lty=2)
	points(self.bendel.delet.fitdat90$generation, 1-(self.bendel.delet.fitdat90$s.window.delet.fitness.total_.1_0 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1_0), col="orange", type="o", pch=16)
	points(self.del.delet.fitdat90$generation, 1-(self.del.delet.fitdat90$s.window.delet.fitness.total_.1_0 / outc.del.delet.fitdat$s.window.delet.fitness.total_.1_0), col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.1_0 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1_0), col="blue", type="o", pch=16)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_.1_0 / outc.del.delet.fitdat$s.window.delet.fitness.total_.1_0), col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.delet.fitdat90$generation, 1-(asex.bendel.delet.fitdat90$s.window.delet.fitness.total_.1_0 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1_0), col="darkorchid2", type="o", pch=16)
	points(asex.del.delet.fitdat90$generation, 1-(asex.del.delet.fitdat90$s.window.delet.fitness.total_.1_0 / outc.del.delet.fitdat$s.window.delet.fitness.total_.1_0), col="darkorchid2", type="o", pch=2, lty=2)
	legend("topleft", c("Outcrossing", "Asex 90%", "Asex 99%", "Self 90%", "Self 99%", "Beneficials present", "No beneficials"), pch=c(15,15,15,15,15,16,2), col=c("green3", "darkorchid2", "blue", "orange", "red", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts", cex.main=0.75)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1_0, col="green3", type="o", pch=16)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.1_0, col="green3", type="o", pch=2, lty=2)
	points(self.bendel.delet.fitdat99$generation, 1-self.bendel.delet.fitdat99$s.window.delet.fitness.total_.1_0, col="red", type="o", pch=16)
	points(self.del.delet.fitdat99$generation, 1-self.del.delet.fitdat99$s.window.delet.fitness.total_.1_0, col="red", type="o", pch=2, lty=2)
	points(self.bendel.delet.fitdat90$generation, 1-self.bendel.delet.fitdat90$s.window.delet.fitness.total_.1_0, col="orange", type="o", pch=16)
	points(self.del.delet.fitdat90$generation, 1-self.del.delet.fitdat90$s.window.delet.fitness.total_.1_0, col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.1_0, col="blue", type="o", pch=16)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_.1_0, col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.delet.fitdat90$generation, 1-asex.bendel.delet.fitdat90$s.window.delet.fitness.total_.1_0, col="darkorchid2", type="o", pch=16)
	points(asex.del.delet.fitdat90$generation, 1-asex.del.delet.fitdat90$s.window.delet.fitness.total_.1_0, col="darkorchid2", type="o", pch=2, lty=2)
	
	#____________________#
	# FULL FITNESS WITH DELS AND BENS
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.5, 1), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative (to outcrossers) load for ALL (ben & del) muts", cex.main=0.75)
	points(self.bendel.bendel.fitdat99$generation, 1-(self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1_0 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1_0), col="red", type="o", pch=16)
	points(self.del.bendel.fitdat99$generation, 1-(self.del.bendel.fitdat99$s.window.bendel.fitness.total_.1_0 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1_0), col="red", type="o", pch=2, lty=2)
	points(self.bendel.bendel.fitdat90$generation, 1-(self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.1_0 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1_0), col="orange", type="o", pch=16)
	points(self.del.bendel.fitdat90$generation, 1-(self.del.bendel.fitdat90$s.window.bendel.fitness.total_.1_0 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1_0), col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1_0 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1_0), col="blue", type="o", pch=16)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.1_0 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1_0), col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.bendel.fitdat90$generation, 1-(asex.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.1_0 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1_0), col="darkorchid2", type="o", pch=16)
	points(asex.del.bendel.fitdat90$generation, 1-(asex.del.bendel.fitdat90$s.window.bendel.fitness.total_.1_0 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1_0), col="darkorchid2", type="o", pch=2, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.5, 1), xlab="Generation", ylab="Load for ALL muts", main="Load for ALL (ben & del) muts", cex.main=0.75)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1_0, col="green3", type="o", pch=16)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1_0, col="green3", type="o", pch=2, lty=2)
	points(self.bendel.bendel.fitdat99$generation, 1-self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1_0, col="red", type="o", pch=16)
	points(self.del.bendel.fitdat99$generation, 1-self.del.bendel.fitdat99$s.window.bendel.fitness.total_.1_0, col="red", type="o", pch=2, lty=2)
	points(self.bendel.bendel.fitdat90$generation, 1-self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.1_0, col="orange", type="o", pch=16)
	points(self.del.bendel.fitdat90$generation, 1-self.del.bendel.fitdat90$s.window.bendel.fitness.total_.1_0, col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1_0, col="blue", type="o", pch=16)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.1_0, col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.bendel.fitdat90$generation, 1-asex.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.1_0, col="darkorchid2", type="o", pch=16)
	points(asex.del.bendel.fitdat90$generation, 1-asex.del.bendel.fitdat90$s.window.bendel.fitness.total_.1_0, col="darkorchid2", type="o", pch=2, lty=2)
	
	#____________________#
	sample.size <- 100
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 2), xlab="Generation", ylab="Relative (to outc) mean total number delet muts per ind", main="Relative (to outcrossers) mean TOTAL number delet muts per ind", cex.main=0.75)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_.1_0 / outc.bendel.totmuts$s.window.total.muts_.1_0, col="red", type="o", pch=16)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_.1_0 / outc.del.totmuts$s.window.total.muts_.1_0, col="red", type="o", pch=2, lty=2)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_.1_0 / outc.bendel.totmuts$s.window.total.muts_.1_0, col="orange", type="o", pch=16)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_.1_0 / outc.del.totmuts$s.window.total.muts_.1_0, col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_.1_0 / outc.bendel.totmuts$s.window.total.muts_.1_0, col="blue", type="o", pch=16)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_.1_0 / outc.del.totmuts$s.window.total.muts_.1_0, col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_.1_0 / outc.bendel.totmuts$s.window.total.muts_.1_0, col="darkorchid2", type="o", pch=16)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_.1_0 / outc.del.totmuts$s.window.total.muts_.1_0, col="darkorchid2", type="o", pch=2, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 50), xlab="Generation", ylab="Mean total number delet muts per ind", main="Mean TOTAL number delet muts per ind", cex.main=0.75)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.1_0/sample.size, col="green3", type="o", pch=16)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.1_0/sample.size, col="green3", type="o", pch=2, lty=2)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_.1_0/sample.size, col="red", type="o", pch=16)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_.1_0/sample.size, col="red", type="o", pch=2, lty=2)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_.1_0/sample.size, col="orange", type="o", pch=16)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_.1_0/sample.size, col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_.1_0/sample.size, col="blue", type="o", pch=16)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_.1_0/sample.size, col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_.1_0/sample.size, col="darkorchid2", type="o", pch=16)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_.1_0/sample.size, col="darkorchid2", type="o", pch=2, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(1, 11), xlab="Generation", ylab="Relative (to outc) mean number fixed delet muts per ind", main="Relative (to outcrossers) mean number FIXED delet muts per ind", cex.main=0.75)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_.1_0 / outc.bendel.fixedmuts$s.window.fixed.muts_.1_0, col="red", type="o", pch=16)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_.1_0 / outc.del.fixedmuts$s.window.fixed.muts_.1_0, col="red", type="o", pch=2, lty=2)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_.1_0 / outc.bendel.fixedmuts$s.window.fixed.muts_.1_0, col="orange", type="o", pch=16)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_.1_0 / outc.del.fixedmuts$s.window.fixed.muts_.1_0, col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_.1_0 / outc.bendel.fixedmuts$s.window.fixed.muts_.1_0, col="blue", type="o", pch=16)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_.1_0 / outc.del.fixedmuts$s.window.fixed.muts_.1_0, col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_.1_0 / outc.bendel.fixedmuts$s.window.fixed.muts_.1_0, col="darkorchid2", type="o", pch=16)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_.1_0 / outc.del.fixedmuts$s.window.fixed.muts_.1_0, col="darkorchid2", type="o", pch=2, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 50), xlab="Generation", ylab="Mean number fixed delet muts per ind", main="Mean number FIXED delet muts per ind", cex.main=0.75)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.1_0/sample.size, col="green3", type="o", pch=16)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.1_0/sample.size, col="green3", type="o", pch=2, lty=2)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_.1_0/sample.size, col="red", type="o", pch=16)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_.1_0/sample.size, col="red", type="o", pch=2, lty=2)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_.1_0/sample.size, col="orange", type="o", pch=16)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_.1_0/sample.size, col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_.1_0/sample.size, col="blue", type="o", pch=16)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_.1_0/sample.size, col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_.1_0/sample.size, col="darkorchid2", type="o", pch=16)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_.1_0/sample.size, col="darkorchid2", type="o", pch=2, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 2), xlab="Generation", ylab="Relative (to outc) mean number poly delet muts per ind", main="Relative (to outcrossers) mean number POLY delet muts per ind", cex.main=0.75)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_.1_0 / outc.bendel.polymuts$s.window.poly.muts_.1_0, col="red", type="o", pch=16)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_.1_0 / outc.del.polymuts$s.window.poly.muts_.1_0, col="red", type="o", pch=2, lty=2)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_.1_0 / outc.bendel.polymuts$s.window.poly.muts_.1_0, col="orange", type="o", pch=16)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_.1_0 / outc.del.polymuts$s.window.poly.muts_.1_0, col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_.1_0 / outc.bendel.polymuts$s.window.poly.muts_.1_0, col="blue", type="o", pch=16)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_.1_0 / outc.del.polymuts$s.window.poly.muts_.1_0, col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_.1_0 / outc.bendel.polymuts$s.window.poly.muts_.1_0, col="darkorchid2", type="o", pch=16)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_.1_0 / outc.del.polymuts$s.window.poly.muts_.1_0, col="darkorchid2", type="o", pch=2, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 50), xlab="Generation", ylab="Mean number poly delet muts per ind", main="Mean number POLY delet muts per ind", cex.main=0.75)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.1_0/sample.size, col="green3", type="o", pch=16)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.1_0/sample.size, col="green3", type="o", pch=2, lty=2)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_.1_0/sample.size, col="red", type="o", pch=16)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_.1_0/sample.size, col="red", type="o", pch=2, lty=2)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_.1_0/sample.size, col="orange", type="o", pch=16)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_.1_0/sample.size, col="orange", type="o", pch=2, lty=2)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_.1_0/sample.size, col="blue", type="o", pch=16)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_.1_0/sample.size, col="blue", type="o", pch=2, lty=2)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_.1_0/sample.size, col="darkorchid2", type="o", pch=16)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_.1_0/sample.size, col="darkorchid2", type="o", pch=2, lty=2)
	
	dev.off()
	
	
	
	
	pdf(paste(c("RelativeandRawLoad_DEL-", figure.basename, "_bySbins.pdf"), collapse=""), width=8, height=8)
	par(mar=c(4,4,1,1), mfrow=c(2,2))
	
	
	# *****
	# ADD TOGETHER EXISTING ONES TO GET THE 4 Nes bins wanted (but could do more bins since have 6 of them)
	
	
	# FITNESS JUST FOR DELS (i.e. FITNESS MEASURES DO NOT INCLudE BENEFICIAL MUTATION EFFECTS)
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 0.8), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative load for delet muts - Self 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.delet.fitdat99$generation, 1-(self.bendel.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1e.04_0), col="#fcae91", type="o", pch=16)
	points(self.bendel.delet.fitdat99$generation, 1-(self.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04), col="#fb6a4a", type="o", pch=16)
	points(self.bendel.delet.fitdat99$generation, 1-(self.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001), col="#de2d26", type="o", pch=16)
	points(self.bendel.delet.fitdat99$generation, 1-(self.bendel.delet.fitdat99$s.window.delet.fitness.total_1_0.01 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_1_0.01), col="#a50f15", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 0.8), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative load for delet muts - Self 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.delet.fitdat99$generation, 1-(self.del.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0 / outc.del.delet.fitdat$s.window.delet.fitness.total_.1e.04_0), col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat99$generation, 1-(self.del.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04 / outc.del.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04), col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat99$generation, 1-(self.del.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001 / outc.del.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001), col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat99$generation, 1-(self.del.delet.fitdat99$s.window.delet.fitness.total_1_0.01 / outc.del.delet.fitdat$s.window.delet.fitness.total_1_0.01), col="#a50f15", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 0.8), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts - self 99% - with bens", cex.main=0.75)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.delet.fitdat99$generation, 1-self.bendel.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0, col="#fcae91", type="o", pch=16)
	points(self.bendel.delet.fitdat99$generation, 1-self.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04, col="#fb6a4a", type="o", pch=16)
	points(self.bendel.delet.fitdat99$generation, 1-self.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001, col="#de2d26", type="o", pch=16)
	points(self.bendel.delet.fitdat99$generation, 1-self.bendel.delet.fitdat99$s.window.delet.fitness.total_1_0.01, col="#a50f15", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 0.8), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts - self 99% - dels only", cex.main=0.75)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat99$generation, 1-self.del.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0, col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat99$generation, 1-self.del.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04, col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat99$generation, 1-self.del.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001, col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat99$generation, 1-self.del.delet.fitdat99$s.window.delet.fitness.total_1_0.01, col="#a50f15", type="o", pch=17, lty=2)
	
	
	
	#_____ Selfing 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 0.8), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative load for delet muts - Self 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.delet.fitdat90$generation, 1-(self.bendel.delet.fitdat90$s.window.delet.fitness.total_.1e.04_0 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1e.04_0), col="#fdbe85", type="o", pch=16)
	points(self.bendel.delet.fitdat90$generation, 1-(self.bendel.delet.fitdat90$s.window.delet.fitness.total_.0.001_.1e.04 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04), col="#fd8d3c", type="o", pch=16)
	points(self.bendel.delet.fitdat90$generation, 1-(self.bendel.delet.fitdat90$s.window.delet.fitness.total_.0.01_.0.001 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001), col="#e6550d", type="o", pch=16)
	points(self.bendel.delet.fitdat90$generation, 1-(self.bendel.delet.fitdat90$s.window.delet.fitness.total_1_0.01 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_1_0.01), col="#a63603", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fdbe85", "#fd8d3c", "#e6550d", "#a63603", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 0.8), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative load for delet muts - Self 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.delet.fitdat90$generation, 1-(self.del.delet.fitdat90$s.window.delet.fitness.total_.1e.04_0 / outc.del.delet.fitdat$s.window.delet.fitness.total_.1e.04_0), col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat90$generation, 1-(self.del.delet.fitdat90$s.window.delet.fitness.total_.0.001_.1e.04 / outc.del.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04), col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat90$generation, 1-(self.del.delet.fitdat90$s.window.delet.fitness.total_.0.01_.0.001 / outc.del.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001), col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat90$generation, 1-(self.del.delet.fitdat90$s.window.delet.fitness.total_1_0.01 / outc.del.delet.fitdat$s.window.delet.fitness.total_1_0.01), col="#a63603", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 0.8), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts - self 90% - with bens", cex.main=0.75)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.delet.fitdat90$generation, 1-self.bendel.delet.fitdat90$s.window.delet.fitness.total_.1e.04_0, col="#fdbe85", type="o", pch=16)
	points(self.bendel.delet.fitdat90$generation, 1-self.bendel.delet.fitdat90$s.window.delet.fitness.total_.0.001_.1e.04, col="#fd8d3c", type="o", pch=16)
	points(self.bendel.delet.fitdat90$generation, 1-self.bendel.delet.fitdat90$s.window.delet.fitness.total_.0.01_.0.001, col="#e6550d", type="o", pch=16)
	points(self.bendel.delet.fitdat90$generation, 1-self.bendel.delet.fitdat90$s.window.delet.fitness.total_1_0.01, col="#a63603", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 0.8), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts - self 90% - dels only", cex.main=0.75)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat90$generation, 1-self.del.delet.fitdat90$s.window.delet.fitness.total_.1e.04_0, col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat90$generation, 1-self.del.delet.fitdat90$s.window.delet.fitness.total_.0.001_.1e.04, col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat90$generation, 1-self.del.delet.fitdat90$s.window.delet.fitness.total_.0.01_.0.001, col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.delet.fitdat90$generation, 1-self.del.delet.fitdat90$s.window.delet.fitness.total_1_0.01, col="#a63603", type="o", pch=17, lty=2)
	
	
	#_____ Asexual 99% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 1), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative load for delet muts - Asex 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1e.04_0), col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04), col="#6baed6", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001), col="#3182bd", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_1_0.01 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_1_0.01), col="#08519c", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 1), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative load for delet muts - Asex 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0 / outc.del.delet.fitdat$s.window.delet.fitness.total_.1e.04_0), col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04 / outc.del.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04), col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001 / outc.del.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001), col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_1_0.01 / outc.del.delet.fitdat$s.window.delet.fitness.total_1_0.01), col="#08519c", type="o", pch=17, lty=2)
	legend("bottomright", cex=0.8, c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#bdd7e7", "#6baed6", "#3182bd", "#08519c", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts - Asex 99% - with bens", cex.main=0.75)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0, col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04, col="#6baed6", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001, col="#3182bd", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_1_0.01, col="#08519c", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts - Asex 99% - dels only", cex.main=0.75)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0, col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04, col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001, col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_1_0.01, col="#08519c", type="o", pch=17, lty=2)
	
	
	#_____ Asexual 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 1), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative load for delet muts - Asex 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1e.04_0), col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04), col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001), col="#756bb1", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-(asex.bendel.delet.fitdat99$s.window.delet.fitness.total_1_0.01 / outc.bendel.delet.fitdat$s.window.delet.fitness.total_1_0.01), col="#54278f", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 1), xlab="Generation", ylab="Relative (to outcrossers) load for delet muts", main="Relative load for delet muts - Asex 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0 / outc.del.delet.fitdat$s.window.delet.fitness.total_.1e.04_0), col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04 / outc.del.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04), col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001 / outc.del.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001), col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-(asex.del.delet.fitdat99$s.window.delet.fitness.total_1_0.01 / outc.del.delet.fitdat$s.window.delet.fitness.total_1_0.01), col="#54278f", type="o", pch=17, lty=2)
	legend("bottomright", cex=0.9, c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts - Asex 90% - with bens", cex.main=0.75)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.delet.fitdat$generation, 1-outc.bendel.delet.fitdat$s.window.delet.fitness.total_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0, col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04, col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001, col="#756bb1", type="o", pch=16)
	points(asex.bendel.delet.fitdat99$generation, 1-asex.bendel.delet.fitdat99$s.window.delet.fitness.total_1_0.01, col="#54278f", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1), xlab="Generation", ylab="Load for delet muts", main="Load for ONLY delet muts - Asex 90% - dels only", cex.main=0.75)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.delet.fitdat$generation, 1-outc.del.delet.fitdat$s.window.delet.fitness.total_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_.1e.04_0, col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_.0.001_.1e.04, col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_.0.01_.0.001, col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.delet.fitdat99$generation, 1-asex.del.delet.fitdat99$s.window.delet.fitness.total_1_0.01, col="#54278f", type="o", pch=17, lty=2)
	
	dev.off()
	
	


	
	#____________________#
	# FULL FITNESS WITH DELS AND BENS
	
	
	pdf(paste(c("RelativeandRawLoad_BENDEL-", figure.basename, "_bySbins.pdf"), collapse=""), width=8, height=8)
	par(mar=c(4,4,1,1), mfrow=c(2,2))
	
	# *****
	# ADD TOGETHER EXISTING ONES TO GET THE 4 Nes bins wanted (but could do more bins since have 6 of them)
	
	
	# FITNESS FOR DELS AND BENEFICIALS
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 0.6), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative load for ALL muts - Self 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.bendel.fitdat99$generation, 1-(self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0), col="#fcae91", type="o", pch=16)
	points(self.bendel.bendel.fitdat99$generation, 1-(self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04), col="#fb6a4a", type="o", pch=16)
	points(self.bendel.bendel.fitdat99$generation, 1-(self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001), col="#de2d26", type="o", pch=16)
	points(self.bendel.bendel.fitdat99$generation, 1-(self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_1_0.01), col="#a50f15", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 0.6), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative load for ALL muts - Self 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.bendel.fitdat99$generation, 1-(self.del.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0), col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat99$generation, 1-(self.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04), col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat99$generation, 1-(self.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001), col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat99$generation, 1-(self.del.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_1_0.01), col="#a50f15", type="o", pch=17, lty=2)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.5, 0.15), xlab="Generation", ylab="load for ALL muts", main="Load for ALL muts - self 99% - with bens", cex.main=0.75)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.bendel.fitdat99$generation, 1-self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0, col="#fcae91", type="o", pch=16)
	points(self.bendel.bendel.fitdat99$generation, 1-self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04, col="#fb6a4a", type="o", pch=16)
	points(self.bendel.bendel.fitdat99$generation, 1-self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001, col="#de2d26", type="o", pch=16)
	points(self.bendel.bendel.fitdat99$generation, 1-self.bendel.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01, col="#a50f15", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.5, 0.15), xlab="Generation", ylab="load for ALL muts", main="Load for ALL muts - self 99% - dels only", cex.main=0.75)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat99$generation, 1-self.del.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0, col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat99$generation, 1-self.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04, col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat99$generation, 1-self.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001, col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat99$generation, 1-self.del.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01, col="#a50f15", type="o", pch=17, lty=2)
	
	
	
	#_____ Selfing 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 0.5), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative load for ALL muts - Self 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.bendel.fitdat90$generation, 1-(self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.1e.04_0 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0), col="#fdbe85", type="o", pch=16)
	points(self.bendel.bendel.fitdat90$generation, 1-(self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.0.001_.1e.04 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04), col="#fd8d3c", type="o", pch=16)
	points(self.bendel.bendel.fitdat90$generation, 1-(self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.0.01_.0.001 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001), col="#e6550d", type="o", pch=16)
	points(self.bendel.bendel.fitdat90$generation, 1-(self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_1_0.01 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_1_0.01), col="#a63603", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 0.5), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative load for ALL muts - Self 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.bendel.fitdat90$generation, 1-(self.del.bendel.fitdat90$s.window.bendel.fitness.total_.1e.04_0 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0), col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat90$generation, 1-(self.del.bendel.fitdat90$s.window.bendel.fitness.total_.0.001_.1e.04 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04), col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat90$generation, 1-(self.del.bendel.fitdat90$s.window.bendel.fitness.total_.0.01_.0.001 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001), col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat90$generation, 1-(self.del.bendel.fitdat90$s.window.bendel.fitness.total_1_0.01 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_1_0.01), col="#a63603", type="o", pch=17, lty=2)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fdbe85", "#fd8d3c", "#e6550d", "#a63603", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-1, 0.2), xlab="Generation", ylab="load for ALL muts", main="Load for ALL muts - self 90% - with bens", cex.main=0.75)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.bendel.fitdat90$generation, 1-self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.1e.04_0, col="#fdbe85", type="o", pch=16)
	points(self.bendel.bendel.fitdat90$generation, 1-self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.0.001_.1e.04, col="#fd8d3c", type="o", pch=16)
	points(self.bendel.bendel.fitdat90$generation, 1-self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_.0.01_.0.001, col="#e6550d", type="o", pch=16)
	points(self.bendel.bendel.fitdat90$generation, 1-self.bendel.bendel.fitdat90$s.window.bendel.fitness.total_1_0.01, col="#a63603", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-1, 0.2), xlab="Generation", ylab="load for ALL muts", main="Load for ALL muts - self 90% - dels only", cex.main=0.75)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat90$generation, 1-self.del.bendel.fitdat90$s.window.bendel.fitness.total_.1e.04_0, col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat90$generation, 1-self.del.bendel.fitdat90$s.window.bendel.fitness.total_.0.001_.1e.04, col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat90$generation, 1-self.del.bendel.fitdat90$s.window.bendel.fitness.total_.0.01_.0.001, col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.bendel.fitdat90$generation, 1-self.del.bendel.fitdat90$s.window.bendel.fitness.total_1_0.01, col="#a63603", type="o", pch=17, lty=2)
	
	
	#_____ Asexual 99% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0.1, 1), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative load for ALL muts - Asex 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0), col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04), col="#6baed6", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001), col="#3182bd", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_1_0.01), col="#08519c", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0.1, 1), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative load for ALL muts - Asex 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0), col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04), col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001), col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_1_0.01), col="#08519c", type="o", pch=17, lty=2)
	legend("bottomright", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#bdd7e7", "#6baed6", "#3182bd", "#08519c", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-1, 1), xlab="Generation", ylab="load for ALL muts", main="Load for ALL muts - Asex 99% - with bens", cex.main=0.75)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0, col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04, col="#6baed6", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001, col="#3182bd", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01, col="#08519c", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-1, 1), xlab="Generation", ylab="load for ALL muts", main="Load for ALL muts - Asex 99% - dels only", cex.main=0.75)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0, col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04, col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001, col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01, col="#08519c", type="o", pch=17, lty=2)
	
	
	#_____ Asexual 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0.2, 1), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative load for ALL muts - Asex 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0), col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04), col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001), col="#756bb1", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-(asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01 / outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_1_0.01), col="#54278f", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0.2, 1), xlab="Generation", ylab="Relative (to outcrossers) load for ALL muts", main="Relative load for ALL muts - Asex 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0), col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04), col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001), col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-(asex.del.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01 / outc.del.bendel.fitdat$s.window.bendel.fitness.total_1_0.01), col="#54278f", type="o", pch=17, lty=2)
	legend("bottomright", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-1, 1), xlab="Generation", ylab="load for ALL muts", main="Load for ALL muts - Asex 90% - with bens", cex.main=0.75)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.bendel.fitdat$generation, 1-outc.bendel.bendel.fitdat$s.window.bendel.fitness.total_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0, col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04, col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001, col="#756bb1", type="o", pch=16)
	points(asex.bendel.bendel.fitdat99$generation, 1-asex.bendel.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01, col="#54278f", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-1, 1), xlab="Generation", ylab="load for ALL muts", main="Load for ALL muts - Asex 90% - dels only", cex.main=0.75)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.bendel.fitdat$generation, 1-outc.del.bendel.fitdat$s.window.bendel.fitness.total_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.1e.04_0, col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.001_.1e.04, col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_.0.01_.0.001, col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.bendel.fitdat99$generation, 1-asex.del.bendel.fitdat99$s.window.bendel.fitness.total_1_0.01, col="#54278f", type="o", pch=17, lty=2)
	
	dev.off()



	
	pdf(paste(c("MutationCounts-POLY-load_bySbins_", figure.basename, ".pdf", width=8, height=8)
	par(mar=c(4,4,1,1), mfrow=c(2,2))
	
	# *****
	# ADD TOGETHER EXISTING ONES TO GET THE 4 Nes bins wanted (but could do more bins since have 6 of them)
	
	
	# FITNESS JUST FOR DELS (i.e. FITNESS MEASURES DO NOT INCLudE BENEFICIAL MUTATION EFFECTS)
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number POLY delet muts - Self 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_.1e.04_0 / outc.bendel.polymuts$s.window.poly.muts_.1e.04_0, col="#fcae91", type="o", pch=16)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_.0.001_.1e.04 / outc.bendel.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=16)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_.0.01_.0.001 / outc.bendel.polymuts$s.window.poly.muts_.0.01_.0.001, col="#de2d26", type="o", pch=16)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_1_0.01 / outc.bendel.polymuts$s.window.poly.muts_1_0.01, col="#a50f15", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number POLY delet muts - Self 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_.1e.04_0 / outc.del.polymuts$s.window.poly.muts_.1e.04_0, col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_.0.001_.1e.04 / outc.del.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_.0.01_.0.001 / outc.del.polymuts$s.window.poly.muts_.0.01_.0.001, col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_1_0.01 / outc.del.polymuts$s.window.poly.muts_1_0.01, col="#a50f15", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY POLY delet muts - self 99% - with bens", cex.main=0.75)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_.1e.04_0, col="#fcae91", type="o", pch=16)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=16)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_.0.01_.0.001, col="#de2d26", type="o", pch=16)
	points(self.bendel.polymuts99$generation, self.bendel.polymuts99$s.window.poly.muts_1_0.01, col="#a50f15", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY POLY delet muts - self 99% - dels only", cex.main=0.75)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_.1e.04_0, col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_.0.01_.0.001, col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.polymuts99$generation, self.del.polymuts99$s.window.poly.muts_1_0.01, col="#a50f15", type="o", pch=17, lty=2)
	
	
	
	#_____ Selfing 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number POLY delet muts - Self 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_.1e.04_0 / outc.bendel.polymuts$s.window.poly.muts_.1e.04_0, col="#fdbe85", type="o", pch=16)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_.0.001_.1e.04 / outc.bendel.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=16)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_.0.01_.0.001 / outc.bendel.polymuts$s.window.poly.muts_.0.01_.0.001, col="#e6550d", type="o", pch=16)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_1_0.01 / outc.bendel.polymuts$s.window.poly.muts_1_0.01, col="#a63603", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fdbe85", "#fd8d3c", "#e6550d", "#a63603", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number POLY delet muts - Self 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_.1e.04_0 / outc.del.polymuts$s.window.poly.muts_.1e.04_0, col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_.0.001_.1e.04 / outc.del.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_.0.01_.0.001 / outc.del.polymuts$s.window.poly.muts_.0.01_.0.001, col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_1_0.01 / outc.del.polymuts$s.window.poly.muts_1_0.01, col="#a63603", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY POLY delet muts - self 90% - with bens", cex.main=0.75)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_.1e.04_0, col="#fdbe85", type="o", pch=16)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=16)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_.0.01_.0.001, col="#e6550d", type="o", pch=16)
	points(self.bendel.polymuts90$generation, self.bendel.polymuts90$s.window.poly.muts_1_0.01, col="#a63603", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY POLY delet muts - self 90% - dels only", cex.main=0.75)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_.1e.04_0, col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_.0.01_.0.001, col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.polymuts90$generation, self.del.polymuts90$s.window.poly.muts_1_0.01, col="#a63603", type="o", pch=17, lty=2)
	
	
	
	
	#_____ Asexual 99% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 5), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number POLY delet muts - asex 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_.1e.04_0 / outc.bendel.polymuts$s.window.poly.muts_.1e.04_0, col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_.0.001_.1e.04 / outc.bendel.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=16)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_.0.01_.0.001 / outc.bendel.polymuts$s.window.poly.muts_.0.01_.0.001, col="#3182bd", type="o", pch=16)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_1_0.01 / outc.bendel.polymuts$s.window.poly.muts_1_0.01, col="#08519c", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#bdd7e7", "#6baed6", "#3182bd", "#08519c", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 5), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number POLY delet muts - asex 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_.1e.04_0 / outc.del.polymuts$s.window.poly.muts_.1e.04_0, col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_.0.001_.1e.04 / outc.del.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_.0.01_.0.001 / outc.del.polymuts$s.window.poly.muts_.0.01_.0.001, col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_1_0.01 / outc.del.polymuts$s.window.poly.muts_1_0.01, col="#08519c", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY POLY delet muts - asex 99% - with bens", cex.main=0.75)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_.1e.04_0, col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=16)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_.0.01_.0.001, col="#3182bd", type="o", pch=16)
	points(asex.bendel.polymuts99$generation, asex.bendel.polymuts99$s.window.poly.muts_1_0.01, col="#08519c", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY POLY delet muts - asex 99% - dels only", cex.main=0.75)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_.1e.04_0, col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_.0.01_.0.001, col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.polymuts99$generation, asex.del.polymuts99$s.window.poly.muts_1_0.01, col="#08519c", type="o", pch=17, lty=2)
	
	
	
	
	#_____ asexual 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number POLY delet muts - asex 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_.1e.04_0 / outc.bendel.polymuts$s.window.poly.muts_.1e.04_0, col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_.0.001_.1e.04 / outc.bendel.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_.0.01_.0.001 / outc.bendel.polymuts$s.window.poly.muts_.0.01_.0.001, col="#756bb1", type="o", pch=16)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_1_0.01 / outc.bendel.polymuts$s.window.poly.muts_1_0.01, col="#54278f", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number POLY delet muts - asex 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_.1e.04_0 / outc.del.polymuts$s.window.poly.muts_.1e.04_0, col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_.0.001_.1e.04 / outc.del.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_.0.01_.0.001 / outc.del.polymuts$s.window.poly.muts_.0.01_.0.001, col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_1_0.01 / outc.del.polymuts$s.window.poly.muts_1_0.01, col="#54278f", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY POLY delet muts - asex 90% - with bens", cex.main=0.75)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.polymuts$generation, outc.bendel.polymuts$s.window.poly.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_.1e.04_0, col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_.0.01_.0.001, col="#756bb1", type="o", pch=16)
	points(asex.bendel.polymuts90$generation, asex.bendel.polymuts90$s.window.poly.muts_1_0.01, col="#54278f", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY POLY delet muts - asex 90% - dels only", cex.main=0.75)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.polymuts$generation, outc.del.polymuts$s.window.poly.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_.1e.04_0, col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_.0.01_.0.001, col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.polymuts90$generation, asex.del.polymuts90$s.window.poly.muts_1_0.01, col="#54278f", type="o", pch=17, lty=2)
	
	
	dev.off()



	
	
	
	
	
	pdf(paste(c("MutationCounts-FIXED-load_bySbins_", figure.basename, ".pdf"), collapse=""), width=8, height=8)
	par(mar=c(4,4,1,1), mfrow=c(2,2))
	
	# *****
	# ADD TOGETHER EXISTING ONES TO GET THE 4 Nes bins wanted (but could do more bins since have 6 of them)
	
	
	# FITNESS JUST FOR DELS (i.e. FITNESS MEASURES DO NOT INCLudE BENEFICIAL MUTATION EFFECTS)
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 100), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number FIXED delet muts - Self 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_.1e.04_0 / outc.bendel.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#fcae91", type="o", pch=16)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_.0.001_.1e.04 / outc.bendel.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=16)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_.0.01_.0.001 / outc.bendel.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#de2d26", type="o", pch=16)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_1_0.01 / outc.bendel.fixedmuts$s.window.fixed.muts_1_0.01, col="#a50f15", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 100), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number FIXED delet muts - Self 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_.1e.04_0 / outc.del.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_.0.001_.1e.04 / outc.del.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_.0.01_.0.001 / outc.del.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_1_0.01 / outc.del.fixedmuts$s.window.fixed.muts_1_0.01, col="#a50f15", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY FIXED delet muts - self 99% - with bens", cex.main=0.75)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_.1e.04_0, col="#fcae91", type="o", pch=16)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=16)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_.0.01_.0.001, col="#de2d26", type="o", pch=16)
	points(self.bendel.fixedmuts99$generation, self.bendel.fixedmuts99$s.window.fixed.muts_1_0.01, col="#a50f15", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY FIXED delet muts - self 99% - dels only", cex.main=0.75)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_.1e.04_0, col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_.0.01_.0.001, col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.fixedmuts99$generation, self.del.fixedmuts99$s.window.fixed.muts_1_0.01, col="#a50f15", type="o", pch=17, lty=2)
	
	
	#_____ Selfing 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 25), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number FIXED delet muts - Self 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_.1e.04_0 / outc.bendel.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#fdbe85", type="o", pch=16)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_.0.001_.1e.04 / outc.bendel.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=16)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_.0.01_.0.001 / outc.bendel.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#e6550d", type="o", pch=16)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_1_0.01 / outc.bendel.fixedmuts$s.window.fixed.muts_1_0.01, col="#a63603", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fdbe85", "#fd8d3c", "#e6550d", "#a63603", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 25), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number FIXED delet muts - Self 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_.1e.04_0 / outc.del.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_.0.001_.1e.04 / outc.del.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_.0.01_.0.001 / outc.del.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_1_0.01 / outc.del.fixedmuts$s.window.fixed.muts_1_0.01, col="#a63603", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY FIXED delet muts - self 90% - with bens", cex.main=0.75)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_.1e.04_0, col="#fdbe85", type="o", pch=16)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=16)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_.0.01_.0.001, col="#e6550d", type="o", pch=16)
	points(self.bendel.fixedmuts90$generation, self.bendel.fixedmuts90$s.window.fixed.muts_1_0.01, col="#a63603", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY FIXED delet muts - self 90% - dels only", cex.main=0.75)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_.1e.04_0, col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_.0.01_.0.001, col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.fixedmuts90$generation, self.del.fixedmuts90$s.window.fixed.muts_1_0.01, col="#a63603", type="o", pch=17, lty=2)
	
	
	#_____ Asexual 99% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 350), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number FIXED delet muts - asex 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_.1e.04_0 / outc.bendel.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_.0.001_.1e.04 / outc.bendel.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=16)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_.0.01_.0.001 / outc.bendel.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#3182bd", type="o", pch=16)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_1_0.01 / outc.bendel.fixedmuts$s.window.fixed.muts_1_0.01, col="#08519c", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#bdd7e7", "#6baed6", "#3182bd", "#08519c", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 350), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number FIXED delet muts - asex 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_.1e.04_0 / outc.del.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_.0.001_.1e.04 / outc.del.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_.0.01_.0.001 / outc.del.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_1_0.01 / outc.del.fixedmuts$s.window.fixed.muts_1_0.01, col="#08519c", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY FIXED delet muts - asex 99% - with bens", cex.main=0.75)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_.1e.04_0, col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=16)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_.0.01_.0.001, col="#3182bd", type="o", pch=16)
	points(asex.bendel.fixedmuts99$generation, asex.bendel.fixedmuts99$s.window.fixed.muts_1_0.01, col="#08519c", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY FIXED delet muts - asex 99% - dels only", cex.main=0.75)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_.1e.04_0, col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_.0.01_.0.001, col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts99$generation, asex.del.fixedmuts99$s.window.fixed.muts_1_0.01, col="#08519c", type="o", pch=17, lty=2)
	
	
	#_____ asexual 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 20), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number FIXED delet muts - asex 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_.1e.04_0 / outc.bendel.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_.0.001_.1e.04 / outc.bendel.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_.0.01_.0.001 / outc.bendel.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#756bb1", type="o", pch=16)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_1_0.01 / outc.bendel.fixedmuts$s.window.fixed.muts_1_0.01, col="#54278f", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 20), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number FIXED delet muts - asex 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_.1e.04_0 / outc.del.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_.0.001_.1e.04 / outc.del.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_.0.01_.0.001 / outc.del.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_1_0.01 / outc.del.fixedmuts$s.window.fixed.muts_1_0.01, col="#54278f", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY FIXED delet muts - asex 90% - with bens", cex.main=0.75)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.fixedmuts$generation, outc.bendel.fixedmuts$s.window.fixed.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_.1e.04_0, col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_.0.01_.0.001, col="#756bb1", type="o", pch=16)
	points(asex.bendel.fixedmuts90$generation, asex.bendel.fixedmuts90$s.window.fixed.muts_1_0.01, col="#54278f", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY FIXED delet muts - asex 90% - dels only", cex.main=0.75)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.fixedmuts$generation, outc.del.fixedmuts$s.window.fixed.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_.1e.04_0, col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_.0.01_.0.001, col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.fixedmuts90$generation, asex.del.fixedmuts90$s.window.fixed.muts_1_0.01, col="#54278f", type="o", pch=17, lty=2)
	
	dev.off()
	





	
	
	pdf(paste(c("MutationCounts-TOTAL-load_bySbins_", figure.basename, ".pdf"), collapse=""), width=8, height=8)
	par(mar=c(4,4,1,1), mfrow=c(2,2))
	
	# *****
	# ADD TOGETHER EXISTING ONES TO GET THE 4 Nes bins wanted (but could do more bins since have 6 of them)
	
	
	# FITNESS JUST FOR DELS (i.e. FITNESS MEASURES DO NOT INCLudE BENEFICIAL MUTATION EFFECTS)
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number TOTAL delet muts - Self 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_.1e.04_0 / outc.bendel.totmuts$s.window.total.muts_.1e.04_0, col="#fcae91", type="o", pch=16)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_.0.001_.1e.04 / outc.bendel.totmuts$s.window.total.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=16)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_.0.01_.0.001 / outc.bendel.totmuts$s.window.total.muts_.0.01_.0.001, col="#de2d26", type="o", pch=16)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_1_0.01 / outc.bendel.totmuts$s.window.total.muts_1_0.01, col="#a50f15", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fcae91", "#fb6a4a", "#de2d26", "#a50f15", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number TOTAL delet muts - Self 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_.1e.04_0 / outc.del.totmuts$s.window.total.muts_.1e.04_0, col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_.0.001_.1e.04 / outc.del.totmuts$s.window.total.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_.0.01_.0.001 / outc.del.totmuts$s.window.total.muts_.0.01_.0.001, col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_1_0.01 / outc.del.totmuts$s.window.total.muts_1_0.01, col="#a50f15", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY TOTAL delet muts - self 99% - with bens", cex.main=0.75)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_.1e.04_0, col="#fcae91", type="o", pch=16)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=16)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_.0.01_.0.001, col="#de2d26", type="o", pch=16)
	points(self.bendel.totmuts99$generation, self.bendel.totmuts99$s.window.total.muts_1_0.01, col="#a50f15", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY TOTAL delet muts - self 99% - dels only", cex.main=0.75)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_.1e.04_0, col="#fcae91", type="o", pch=17, lty=2)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_.0.001_.1e.04, col="#fb6a4a", type="o", pch=17, lty=2)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_.0.01_.0.001, col="#de2d26", type="o", pch=17, lty=2)
	points(self.del.totmuts99$generation, self.del.totmuts99$s.window.total.muts_1_0.01, col="#a50f15", type="o", pch=17, lty=2)
	
	
	#_____ Selfing 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number TOTAL delet muts - Self 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_.1e.04_0 / outc.bendel.totmuts$s.window.total.muts_.1e.04_0, col="#fdbe85", type="o", pch=16)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_.0.001_.1e.04 / outc.bendel.totmuts$s.window.total.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=16)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_.0.01_.0.001 / outc.bendel.totmuts$s.window.total.muts_.0.01_.0.001, col="#e6550d", type="o", pch=16)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_1_0.01 / outc.bendel.totmuts$s.window.total.muts_1_0.01, col="#a63603", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#fdbe85", "#fd8d3c", "#e6550d", "#a63603", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number TOTAL delet muts - Self 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_.1e.04_0 / outc.del.totmuts$s.window.total.muts_.1e.04_0, col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_.0.001_.1e.04 / outc.del.totmuts$s.window.total.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_.0.01_.0.001 / outc.del.totmuts$s.window.total.muts_.0.01_.0.001, col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_1_0.01 / outc.del.totmuts$s.window.total.muts_1_0.01, col="#a63603", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY TOTAL delet muts - self 90% - with bens", cex.main=0.75)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_.1e.04_0, col="#fdbe85", type="o", pch=16)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=16)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_.0.01_.0.001, col="#e6550d", type="o", pch=16)
	points(self.bendel.totmuts90$generation, self.bendel.totmuts90$s.window.total.muts_1_0.01, col="#a63603", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY TOTAL delet muts - self 90% - dels only", cex.main=0.75)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_.1e.04_0, col="#fdbe85", type="o", pch=17, lty=2)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_.0.001_.1e.04, col="#fd8d3c", type="o", pch=17, lty=2)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_.0.01_.0.001, col="#e6550d", type="o", pch=17, lty=2)
	points(self.del.totmuts90$generation, self.del.totmuts90$s.window.total.muts_1_0.01, col="#a63603", type="o", pch=17, lty=2)
	
	
	#_____ Asexual 99% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 5), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number TOTAL delet muts - asex 99% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_.1e.04_0 / outc.bendel.totmuts$s.window.total.muts_.1e.04_0, col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_.0.001_.1e.04 / outc.bendel.totmuts$s.window.total.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=16)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_.0.01_.0.001 / outc.bendel.totmuts$s.window.total.muts_.0.01_.0.001, col="#3182bd", type="o", pch=16)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_1_0.01 / outc.bendel.totmuts$s.window.total.muts_1_0.01, col="#08519c", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#bdd7e7", "#6baed6", "#3182bd", "#08519c", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 5), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number TOTAL delet muts - asex 99% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_.1e.04_0 / outc.del.totmuts$s.window.total.muts_.1e.04_0, col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_.0.001_.1e.04 / outc.del.totmuts$s.window.total.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_.0.01_.0.001 / outc.del.totmuts$s.window.total.muts_.0.01_.0.001, col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_1_0.01 / outc.del.totmuts$s.window.total.muts_1_0.01, col="#08519c", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY TOTAL delet muts - asex 99% - with bens", cex.main=0.75)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_.1e.04_0, col="#bdd7e7", type="o", pch=16)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=16)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_.0.01_.0.001, col="#3182bd", type="o", pch=16)
	points(asex.bendel.totmuts99$generation, asex.bendel.totmuts99$s.window.total.muts_1_0.01, col="#08519c", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY TOTAL delet muts - asex 99% - dels only", cex.main=0.75)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_.1e.04_0, col="#bdd7e7", type="o", pch=17, lty=2)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_.0.001_.1e.04, col="#6baed6", type="o", pch=17, lty=2)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_.0.01_.0.001, col="#3182bd", type="o", pch=17, lty=2)
	points(asex.del.totmuts99$generation, asex.del.totmuts99$s.window.total.muts_1_0.01, col="#08519c", type="o", pch=17, lty=2)
	
	
	#_____ asexual 90% ______#
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number TOTAL delet muts - asex 90% to outc - with bens", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_.1e.04_0 / outc.bendel.totmuts$s.window.total.muts_.1e.04_0, col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_.0.001_.1e.04 / outc.bendel.totmuts$s.window.total.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_.0.01_.0.001 / outc.bendel.totmuts$s.window.total.muts_.0.01_.0.001, col="#756bb1", type="o", pch=16)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_1_0.01 / outc.bendel.totmuts$s.window.total.muts_1_0.01, col="#54278f", type="o", pch=16)
	legend("topleft", c("Ne*s 0 -- 10^-4", "Ne*s 10^-4 -- 0.001", "Ne*s 0.001 -- 0.01", "Ne*s 0.01 -- 1", "with bens", "delet only"), pch=c(15,15,15,15,16,17), col=c("#cbc9e2", "#9e9ac8", "#756bb1", "#54278f", "black", "black"))
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(-0.1, 2), xlab="Generation", ylab="Relative (to outcrossers) number delet muts", main="Relative number TOTAL delet muts - asex 90% to outc - dels only", cex.main=0.75)
	abline(h=0, col="gray", lty=3)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_.1e.04_0 / outc.del.totmuts$s.window.total.muts_.1e.04_0, col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_.0.001_.1e.04 / outc.del.totmuts$s.window.total.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_.0.01_.0.001 / outc.del.totmuts$s.window.total.muts_.0.01_.0.001, col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_1_0.01 / outc.del.totmuts$s.window.total.muts_1_0.01, col="#54278f", type="o", pch=17, lty=2)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number delet muts", main="Raw number for ONLY TOTAL delet muts - asex 90% - with bens", cex.main=0.75)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.1e.04_0, col="#bae4b3", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.0.001_.1e.04, col="#74c476", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_.0.01_.0.001, col="#31a354", type="o", pch=16)
	points(outc.bendel.totmuts$generation, outc.bendel.totmuts$s.window.total.muts_1_0.01, col="#006d2c", type="o", pch=16)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_.1e.04_0, col="#cbc9e2", type="o", pch=16)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=16)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_.0.01_.0.001, col="#756bb1", type="o", pch=16)
	points(asex.bendel.totmuts90$generation, asex.bendel.totmuts90$s.window.total.muts_1_0.01, col="#54278f", type="o", pch=16)
	
	plot(0, type="n", xlim=c(xlimits), ylim=c(0, 1000), xlab="Generation", ylab="Raw number for delet muts", main="Raw number for ONLY TOTAL delet muts - asex 90% - dels only", cex.main=0.75)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.1e.04_0, col="#bae4b3", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.0.001_.1e.04, col="#74c476", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_.0.01_.0.001, col="#31a354", type="o", pch=17, lty=2)
	points(outc.del.totmuts$generation, outc.del.totmuts$s.window.total.muts_1_0.01, col="#006d2c", type="o", pch=17, lty=2)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_.1e.04_0, col="#cbc9e2", type="o", pch=17, lty=2)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_.0.001_.1e.04, col="#9e9ac8", type="o", pch=17, lty=2)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_.0.01_.0.001, col="#756bb1", type="o", pch=17, lty=2)
	points(asex.del.totmuts90$generation, asex.del.totmuts90$s.window.total.muts_1_0.01, col="#54278f", type="o", pch=17, lty=2)
	
	dev.off()
}





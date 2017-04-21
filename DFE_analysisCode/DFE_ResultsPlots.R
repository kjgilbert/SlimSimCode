library(Hmisc) # for plotting the error bars
library(scales)


outfile.name <- "Results_Apr19_Messerlike_10Ngens"
ne.dat.file <- "NeEstimates_Apr19_MesserLikeParams_10Ngens-N10000.csv"
est.dat <- read.csv("est_dfe_results_Apr19_Messerlike_10Ngens-N10000.csv", header=TRUE)
dominance.val <- 0.3
gamma.shape <- 0.2


#______________________________________________________________________________________________________________#

## outputs from DFE

#______________________________________________________________________________________________________________#


# aggregate all the results into a single data frame (make into a file) then can subset as needed

setwd("~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/April2017_NewMesserLikeSims/OutcrossersTo10N/Outputs_10Ngens_Apr19")


data.rows <- as.numeric(system("ls DFE_Results/out*prop* | wc -l", intern=TRUE))
mut.props.dat <- data.frame(matrix(ncol=10, nrow=data.rows))	# 12 columns from the DFE outputs, but only 4 are the proportions in each category, rest are for naming/sorting purposes
names(mut.props.dat) <- c("group.ID", "pop.size", "genome.size", "mut.types", "mating.system", "replicate", "prop.0.1", "prop.1.10", "prop.10.100", "prop.100.inf")

files <- system("ls DFE_Results/out*prop*", intern=TRUE)
iterate <- 1
for(i in files){
	temp.dat <- read.table(i)
	row.info <- unlist(strsplit(i, split="_"))

	mut.props.dat[iterate,] <- c(paste(row.info[8:14], collapse="_"), row.info[11:15], temp.dat[c(3,6,9,12)])

	iterate <- iterate + 1
}

write.table(mut.props.dat, file=paste(c("2epoch_ProportionMutations_", outfile.name, ".csv"), collapse=""), sep=",", col.names=TRUE)
#______________________________________________________________________________________________________________#


data.rows <- as.numeric(system("ls DFE_Results/out*alpha_omega* | wc -l", intern=TRUE))
alpha.omega.dat <- data.frame(matrix(ncol=10, nrow=data.rows))	# 12 columns from the DFE outputs, but only 4 are the proportions in each category, rest are for naming/sorting purposes
names(alpha.omega.dat) <- c("group.ID", "pop.size", "genome.size", "mut.types", "mating.system", "replicate", "lambda", "selected.divergence", "alpha", "omega_a")

files <- system("ls DFE_Results/out*alpha_omega*", intern=TRUE)
iterate <- 1
for(i in files){
	temp.dat <- read.table(i)
	row.info <- unlist(strsplit(i, split="_"))

	alpha.omega.dat[iterate,] <- c(paste(row.info[6:12], collapse="_"), row.info[9:13], temp.dat[c(2,4,6,8)])
	
	iterate <- iterate + 1
}

write.table(alpha.omega.dat, file=paste(c("2epoch_AlphaOmega_", outfile.name, ".csv"), collapse=""), sep=",", col.names=TRUE)
#______________________________________________________________________________________________________________#



#______________________________________________________________________________________________________________#

##		ADD Ne DATA
Ne.dat <- read.csv(ne.dat.file, stringsAsFactors=FALSE)
Ne.dat$group.ID <- matrix(unlist(strsplit(Ne.dat$file, split="_rep")), ncol=2, byrow=TRUE)[,1]
mean.Ne.dat <- aggregate(Ne.dat$Ne.neutral, by=list(Ne.dat$group.ID), FUN=mean)
names(mean.Ne.dat) <- c("scenario", "Ne")

##		ADD and organize remaining DATA

mut.props.dat <- read.csv(paste(c("2epoch_ProportionMutations_", outfile.name, ".csv"), collapse=""), header=TRUE)
alpha.omega.dat <- read.csv(paste(c("2epoch_AlphaOmega_", outfile.name, ".csv"), collapse=""), header=TRUE)


# proportion mutations results:

cols.to.aggregate <- c("prop.0.1", "prop.1.10", "prop.10.100", "prop.100.inf")
mean.props.by.group <- aggregate(mut.props.dat[cols.to.aggregate], by=list(mut.props.dat$group.ID), FUN=mean)
var.props.by.group <- aggregate(mut.props.dat[cols.to.aggregate], by=list(mut.props.dat$group.ID), FUN=var)
sd_mean.props.by.group <- aggregate(mut.props.dat[cols.to.aggregate], by=list(mut.props.dat$group.ID), FUN=sd)
ci95_mean.props.by.group <- (sd_mean.props.by.group[cols.to.aggregate]/sqrt(10)) * 2.26
names(ci95_mean.props.by.group) <- c("ci95_prop.0.1", "ci95_prop.1.10", "ci95_prop.10.100", "ci95_prop.100.inf")
mean.proportion.muts <- cbind(mean.props.by.group, ci95_mean.props.by.group)


# alpha omega results:

cols.to.aggregate <- c("lambda", "selected.divergence", "alpha", "omega_a")
mean.alphaomega.by.group <- aggregate(alpha.omega.dat[cols.to.aggregate], by=list(alpha.omega.dat$group.ID), FUN=mean)
var.alphaomega.by.group <- aggregate(alpha.omega.dat[cols.to.aggregate], by=list(alpha.omega.dat$group.ID), FUN=var)
sd_mean.alphaomega.by.group <- aggregate(alpha.omega.dat[cols.to.aggregate], by=list(alpha.omega.dat$group.ID), FUN=sd)
ci95_mean.alphaomega.by.group <- (sd_mean.alphaomega.by.group[cols.to.aggregate]/sqrt(10)) * 2.26
names(ci95_mean.alphaomega.by.group) <- c("ci95_lambda", "ci95_selected.divergence", "ci95_alpha", "ci95_omega_a")
mean.alpha.omega <- cbind(mean.alphaomega.by.group, ci95_mean.alphaomega.by.group)

# est_dfe results:

cols.to.aggregate <- c("N1.0", "N2.0", "t2.0", "Nw.0", "f0.0", "L.0", "N1.1", "N2.1", "t2.1", "Nw.1", "b", "Es", "f0.1", "L.1")
mean.est.dat.by.group <- aggregate(est.dat[cols.to.aggregate], by=list(est.dat$group.ID), FUN=mean)
var.est.dat.by.group <- aggregate(est.dat[cols.to.aggregate], by=list(est.dat$group.ID), FUN=var)
sd_mean.est.dat.by.group <- aggregate(est.dat[cols.to.aggregate], by=list(est.dat$group.ID), FUN=sd)
ci95_mean.est.dat.by.group <- (sd_mean.est.dat.by.group[cols.to.aggregate]/sqrt(10)) * 2.26
names(ci95_mean.est.dat.by.group) <- c("ci95_N1.0","ci95_N2.0", "ci95_t2.0", "ci95_Nw.0","ci95_f0.0","ci95_L.0","ci95_N1.1", "ci95_N2.1", "ci95_t2.1", "ci95_Nw.1", "ci95_b", "ci95_Es", "ci95_f0.1", "ci95_L.1")
mean.est.dat <- cbind(mean.est.dat.by.group, ci95_mean.est.dat.by.group)

simulation.list <- c(paste(unique(mean.proportion.muts$Group.1)))
cols.list <- c("#a6cee3", "#ff7f00", "#1f78b4", "#b2df8a", "#e31a1c", "#33a02c", "#fb9a99", "#fdbf6f", "#cab2d6", "#006d2c", "#2ca25f", "#66c2a4", "#99d8c9", "#ccece6", "#edf8fb", "#810f7c", "#8856a7", "#8c96c6", "#9ebcda", "#bfd3e6", "#edf8fb", "#006837", "#31a354", "#78c679", "#addd8e", "#d9f0a3", "#ffffcc", "#0868ac", "#43a2ca", "#7bccc4", "#a8ddb5", "#ccebc5", "#f0f9e8", "#a50f15", "#de2d26", "#fb6a4a", "#fc9272", "#fcbba1", "#fee5d9")	



#______________________________________________________________________________________________________________#
##		PLOTTING FUNCTIONS
#___________________________________#
# from real Ne data:
ne.bar.heights <- function(Ne, h, gamma.mean){
	bar1 <- pgamma(((2*h)*1)/Ne, shape=gamma.shape, scale=gamma.mean/gamma.shape)
	temp <- pgamma(((2*h)*10)/Ne, shape=gamma.shape, scale=gamma.mean/gamma.shape)
	bar2 <- temp - bar1
	temp <- pgamma(((2*h)*100)/Ne, shape=gamma.shape, scale=gamma.mean/gamma.shape)
	bar3 <- temp - (bar2 + bar1)
	bar4 <- 1- (bar3 + bar2 + bar1)	
	return(c(bar1, bar2, bar3, bar4))
}
#___________________________________#
h.bar.heights <- function(h, gamma.mean){
	bar1 <- pgamma(((2*h)*1)/10000, shape=gamma.shape, scale=gamma.mean/gamma.shape)
	temp <- pgamma(((2*h)*10)/10000, shape=gamma.shape, scale=gamma.mean/gamma.shape)
	bar2 <- temp - bar1
	temp <- pgamma(((2*h)*100)/10000, shape=gamma.shape, scale=gamma.mean/gamma.shape)
	bar3 <- temp - (bar2 + bar1)
	bar4 <- 1- (bar3 + bar2 + bar1)	
	return(c(bar1, bar2, bar3, bar4))
}
#___________________________________#


#______________________________________________________________________________________________________________#

##		LINE GRAPH

#______________________________________________________________________________________________________________#


dat <- mean.proportion.muts

pdf(paste(c("DFE_LinePlot_", outfile.name,".pdf"), collapse=""), width=7, height=4.5); par(mar=c(4,4,1,1))

##	# make a plot for each generation sampled
##	generations <- unique(mut.props.dat$generation)
##	for(i in generations){
	# get this gen's dat
##		dat <- mean.proportion.muts[grep(i, mean.proportion.muts$Group.1) ,]

plot(0, type="n", xlim=c(0.75,4.25), ylim=c(0,0.9), xlab="Ne*s", ylab="Proportion mutations", xaxt="n", main=paste(outfile.name))
axis(side=1, at=1:4, labels=c("0-1", "1-10", "10-100", "100-inf"))

iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1),]
	if(length(grep("d0005", temp.dat$Group.1)) == 1){ s.d.mean <- 0.0005 ; line.type <- 2 ; point.type <- 6 }
	if(length(grep("d01", temp.dat$Group.1)) == 1){ s.d.mean <- 0.01 ; line.type <- 1 ; point.type <- 1 }
	
	predicted.gamma <- h.bar.heights(h=dominance.val, gamma.mean=s.d.mean)
	points(1:4, predicted.gamma, lty=line.type, col="gray60", type="o", pch=4, lwd=1.5)
	points(1:4, temp.dat[1, c("prop.0.1", "prop.1.10", "prop.10.100", "prop.100.inf")], col=cols.list[iterate], type="o", lwd=2.5, lty= line.type, pch= point.type)
	iterate <- iterate + 1
}
legend("topright", col=c(rep("gray60",3), cols.list[1:length(simulation.list)]), c("Modelled gamma", "s_d = 0.01", "s_d = 0.0005", simulation.list), pch=c(4, 1, 6, rep(15, length(simulation.list))), bg="white", ncol=1, cex=0.8, lty=c(0,1,2, rep(0,length(simulation.list))), pt.lwd=2)

dev.off()
#______________________________________________________________________________________________________________#







## 		BAR PLOT

#______________________________________________________________________________________________________________#


dat <- mean.proportion.muts


#___________________________________#
# make bar plot function
bars <- function(x, y, width, color){
	# y = height of bar
	# x = center location of bar
	polygon(x=c((x-width), (x-width), (x+width), (x+width)), y=c(0,y,y,0), col=color, border=TRUE)
}
wid=0.25
predict.wid <- wid/2

dist.between <- length(simulation.list)+1
bin1 <- 1
bin2 <- bin1 + dist.between
bin3 <- bin2 + dist.between
bin4 <- bin3 + dist.between
axis.at <- seq(bin1, bin4, length.out=4)

spots <- seq(bin1-length(simulation.list)/2,100, by=0.5)
opacity <- 0.8


pdf(paste(c("DFE_NesBins_", outfile.name,".pdf"), collapse=""), width=2*length(simulation.list), height=4.5)
par(mar=c(4,4,1,0.2))

plot(0, type="n", xlim=c(bin1-length(simulation.list)/2,bin4+length(simulation.list)/2), ylim=c(0,0.7), xlab="Ne*s", ylab="Proportion mutations", xaxt="n", main=paste(outfile.name))
axis(side=1, at= axis.at, labels=c("0-1", "1-10", "10-100", "100-inf"))
legend("topright", pch=c(45, rep(15,length(simulation.list))), col=c("darkgray", cols.list), c("Predicted", simulation.list), bg="white", ncol=1, cex=0.5)

iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1),]
	if(length(grep("d0005", temp.dat$Group.1)) == 1){ s.d.mean <- 0.0005 ; line.type <- 2 ; point.type <- 6 }
	if(length(grep("d01", temp.dat$Group.1)) == 1){ s.d.mean <- 0.01 ; line.type <- 1 ; point.type <- 1 }
	# true gamma:
	h.bars <- h.bar.heights(h=dominance.val, s.d.mean)
	# bins	
	bars(bin1+spots[iterate], temp.dat[, c("prop.0.1")], width=wid, col=cols.list[iterate])
	bars(bin2+spots[iterate], temp.dat[, c("prop.1.10")], width=wid, col=cols.list[iterate])
	bars(bin3+spots[iterate], temp.dat[, c("prop.10.100")], width=wid, col=cols.list[iterate])
	bars(bin4+spots[iterate], temp.dat[, c("prop.100.inf")], width=wid, col=cols.list[iterate])
	# confidence intervals
	segments(x0=bin1+spots[iterate], y0=temp.dat[,c("prop.0.1")]-temp.dat[,c("ci95_prop.0.1")], x1=bin1+spots[iterate], y1=temp.dat[,c("prop.0.1")]+temp.dat[,c("ci95_prop.0.1")], lwd=1, col="gray20")
	segments(x0=bin2+spots[iterate], y0=temp.dat[,c("prop.1.10")]-temp.dat[,c("ci95_prop.1.10")], x1=bin2+spots[iterate], y1=temp.dat[,c("prop.1.10")]+temp.dat[,c("ci95_prop.1.10")], lwd=1, col="gray20")
	segments(x0=bin3+spots[iterate], y0=temp.dat[,c("prop.10.100")]-temp.dat[,c("ci95_prop.10.100")], x1=bin3+spots[iterate], y1=temp.dat[,c("prop.10.100")]+temp.dat[,c("ci95_prop.10.100")], lwd=1, col="gray20")
	segments(x0=bin4+spots[iterate], y0=temp.dat[,c("prop.100.inf")]-temp.dat[,c("ci95_prop.100.inf")], x1=bin4+spots[iterate], y1=temp.dat[,c("prop.100.inf")]+temp.dat[,c("ci95_prop.100.inf")], lwd=1, col="gray20")
	# true gamma just plotted overtop instead:
	points(bin1+spots[iterate], h.bars[1], pch="-", cex=1.75, col=alpha("antiquewhite1", opacity))
	points(bin2+spots[iterate], h.bars[2], pch="-", cex=1.75, col=alpha("antiquewhite1", opacity))
	points(bin3+spots[iterate], h.bars[3], pch="-", cex=1.75, col=alpha("antiquewhite1", opacity))
	points(bin4+spots[iterate], h.bars[4], pch="-", cex=1.75, col=alpha("antiquewhite1", opacity))
	points(bin1+spots[iterate], h.bars[1], pch="-", cex=0.75, col="black")
	points(bin2+spots[iterate], h.bars[2], pch="-", cex=0.75, col="black")
	points(bin3+spots[iterate], h.bars[3], pch="-", cex=0.75, col="black")
	points(bin4+spots[iterate], h.bars[4], pch="-", cex=0.75, col="black")

	iterate <- iterate + 1
}

# with NE adjusted predictions:
plot(0, type="n", xlim=c(bin1-length(simulation.list)/2,bin4+length(simulation.list)/2), ylim=c(0,0.7), xlab="Ne*s", ylab="Proportion mutations", xaxt="n", main=paste(c("Ne-predicted", outfile.name), collapse=" "))
axis(side=1, at= axis.at, labels=c("0-1", "1-10", "10-100", "100-inf"))
legend("topright", pch=c(45, rep(15,length(simulation.list))), col=c("darkgray", cols.list), c("Ne Predicted", simulation.list), bg="white", ncol=1, cex=0.5)

iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1),]
	if(length(grep("d0005", temp.dat$Group.1)) == 1){ s.d.mean <- 0.0005 ; line.type <- 2 ; point.type <- 6 }
	if(length(grep("d01", temp.dat$Group.1)) == 1){ s.d.mean <- 0.01 ; line.type <- 1 ; point.type <- 1 }
	Ne <- mean.Ne.dat$Ne[grep(i, mean.Ne.dat$scenario)]
	# true gamma:
	Ne.bars <- ne.bar.heights(Ne, dominance.val, s.d.mean)
	# bins	
	bars(bin1+spots[iterate], temp.dat[, c("prop.0.1")], width=wid, col=cols.list[iterate])
	bars(bin2+spots[iterate], temp.dat[, c("prop.1.10")], width=wid, col=cols.list[iterate])
	bars(bin3+spots[iterate], temp.dat[, c("prop.10.100")], width=wid, col=cols.list[iterate])
	bars(bin4+spots[iterate], temp.dat[, c("prop.100.inf")], width=wid, col=cols.list[iterate])
	# confidence intervals
	segments(x0=bin1+spots[iterate], y0=temp.dat[,c("prop.0.1")]-temp.dat[,c("ci95_prop.0.1")], x1=bin1+spots[iterate], y1=temp.dat[,c("prop.0.1")]+temp.dat[,c("ci95_prop.0.1")], lwd=1, col="gray20")
	segments(x0=bin2+spots[iterate], y0=temp.dat[,c("prop.1.10")]-temp.dat[,c("ci95_prop.1.10")], x1=bin2+spots[iterate], y1=temp.dat[,c("prop.1.10")]+temp.dat[,c("ci95_prop.1.10")], lwd=1, col="gray20")
	segments(x0=bin3+spots[iterate], y0=temp.dat[,c("prop.10.100")]-temp.dat[,c("ci95_prop.10.100")], x1=bin3+spots[iterate], y1=temp.dat[,c("prop.10.100")]+temp.dat[,c("ci95_prop.10.100")], lwd=1, col="gray20")
	segments(x0=bin4+spots[iterate], y0=temp.dat[,c("prop.100.inf")]-temp.dat[,c("ci95_prop.100.inf")], x1=bin4+spots[iterate], y1=temp.dat[,c("prop.100.inf")]+temp.dat[,c("ci95_prop.100.inf")], lwd=1, col="gray20")
	# true gamma just plotted overtop instead:
	points(bin1+spots[iterate], Ne.bars[1], pch="-", cex=1.75, col=alpha("antiquewhite1", opacity))
	points(bin2+spots[iterate], Ne.bars[2], pch="-", cex=1.75, col=alpha("antiquewhite1", opacity))
	points(bin3+spots[iterate], Ne.bars[3], pch="-", cex=1.75, col=alpha("antiquewhite1", opacity))
	points(bin4+spots[iterate], Ne.bars[4], pch="-", cex=1.75, col=alpha("antiquewhite1", opacity))
	points(bin1+spots[iterate], Ne.bars[1], pch="-", cex=0.75, col="black")
	points(bin2+spots[iterate], Ne.bars[2], pch="-", cex=0.75, col="black")
	points(bin3+spots[iterate], Ne.bars[3], pch="-", cex=0.75, col="black")
	points(bin4+spots[iterate], Ne.bars[4], pch="-", cex=0.75, col="black")

	iterate <- iterate + 1
}
dev.off()
#______________________________________________________________________________________________________________#











## POP GEN PARAM PLOTS (Ne 2, t 2)

#______________________________________________________________________________________________________________#
dat <- mean.est.dat

	
pdf(paste(c("DFE_PopGenParamResults_", outfile.name,".pdf"), collapse=""), width=length(simulation.list)/1.25, height=4)
par(mfrow=c(1,2), mar=c(1,4,1,0.5))

plot(0, type="n", xlim=c(0.5,length(simulation.list)+0.5), ylim=c(0,1000), xlab="", ylab="Estimated Population Size", xaxt="n", main="Ne")
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1) ,]
	errbar(add=TRUE, iterate, temp.dat$N2.1, yplus=temp.dat$N2.1+temp.dat$ci95_N2.1, yminus=temp.dat$N2.1-temp.dat$ci95_N2.1, pch=20, errbar.col=cols.list[iterate], cex=0.5)
	points(iterate, temp.dat$N2.1, col=cols.list[iterate], pch=16, cex=1.5)
	iterate <- iterate + 1
}
plot(0, type="n", xlim=c(0.5,length(simulation.list)+0.5), ylim=c(0,2500), xlab="", ylab="Estimated time since pop size change", xaxt="n", main="time since change")
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1) ,]
	errbar(add=TRUE, iterate, temp.dat$t2.1, yplus=temp.dat$t2.1+temp.dat$ci95_t2.1, yminus=temp.dat$t2.1-temp.dat$ci95_t2.1, pch=20, errbar.col=cols.list[iterate], cex=0.5)
	points(iterate, temp.dat$t2.1, col=cols.list[iterate], pch=16, cex=1.5)
	iterate <- iterate + 1
}
dev.off()
#______________________________________________________________________________________________________________#











## ADAPTIVE RESULTS PLOTS (alpha, omega_a)

#______________________________________________________________________________________________________________#
dat <- mean.alpha.omega
dat2 <- mean.est.dat
		
pdf(paste(c("DFE_AdaptiveParamResults_", outfile.name,".pdf"), collapse=""), width=length(simulation.list)+4, height=4)
par(mfrow=c(1,5), mar=c(1,4,1,0.5))

plot(0, type="n", xlim=c(0.5,length(simulation.list)+0.5), ylim=c(-0.5,1), xlab="", ylab="Alpha", xaxt="n", main="alpha")
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1) ,]
	errbar(add=TRUE, iterate, temp.dat$alpha, yplus=temp.dat$alpha+temp.dat$ci95_alpha, yminus=temp.dat$alpha-temp.dat$ci95_alpha, pch=20, errbar.col=cols.list[iterate], cex=0.5)
	points(iterate, temp.dat$alpha, col=cols.list[iterate], pch=16, cex=1.5)
	iterate <- iterate + 1
}

plot(0, type="n", xlim=c(0.5,length(simulation.list)+0.5), ylim=c(-0.3,1), xlab="", ylab="omega_a", xaxt="n", main="omega")
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1) ,]
	errbar(add=TRUE, iterate, temp.dat$omega_a, yplus=temp.dat$omega_a+temp.dat$ci95_omega_a, yminus=temp.dat$omega_a-temp.dat$ci95_omega_a, pch=20, errbar.col=cols.list[iterate], cex=0.5)
	points(iterate, temp.dat$omega_a, col= cols.list[iterate], pch=16, cex=1.5)
	iterate <- iterate + 1
}
legend("topright", pch=15, col=cols.list, c(simulation.list), bg="white", ncol=1, cex=0.75)

plot(0, type="n", xlim=c(0.5,length(simulation.list)+0.5), ylim=c(-0.0025,0.005), xlab="", ylab="lambda", xaxt="n", main="lambda")
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1) ,]
	errbar(add=TRUE, iterate, temp.dat$lambda, yplus=temp.dat$lambda+temp.dat$ci95_lambda, yminus=temp.dat$lambda-temp.dat$ci95_lambda, pch=20, errbar.col=cols.list[iterate], cex=0.5)
	points(iterate, temp.dat$lambda, col= cols.list[iterate], pch=16, cex=1.5)
	iterate <- iterate + 1
}

plot(0, type="n", xlim=c(0.5,length(simulation.list)+0.5), ylim=c(-10,0.1), xlab="", ylab="Estimated s", xaxt="n", main="s")
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat2[grep(i, dat$Group.1) ,]
	errbar(add=TRUE, iterate, temp.dat$Es, yplus=temp.dat$Es+temp.dat$ci95_Es, yminus=temp.dat$Es-temp.dat$ci95_Es, pch=20, errbar.col=cols.list[iterate], cex=0.5)
	points(iterate, temp.dat$Es, col= cols.list[iterate], pch=16, cex=1.5)
	iterate <- iterate + 1
}

plot(0, type="n", xlim=c(0.5,length(simulation.list)+0.5), ylim=c(-0.1,10), xlab="", ylab="Estimated beta", xaxt="n", main="beta")
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat2[grep(i, dat$Group.1) ,]
	errbar(add=TRUE, iterate, temp.dat$b, yplus=temp.dat$b+temp.dat$ci95_b, yminus=temp.dat$b-temp.dat$ci95_b, pch=20, errbar.col=cols.list[iterate], cex=0.5)
	points(iterate, temp.dat$b, col= cols.list[iterate], pch=16, cex=1.5)
	iterate <- iterate + 1
}
dev.off()



library(Hmisc) # for plotting the error bars

setwd("~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/April2017_NewMesserLikeSims/OutcrossersTo10N/Outputs_10Ngens_Apr19")

outfile.name <- "Results_Apr19_Messerlike_10Ngens"
del.summ.stats <- "SummStats_deletsOnly_April19_MesserLikeParams_10Ngens-N10000.csv"
ben.summ.stats <- "SummStats_withBens_April19_MesserLikeParams_10Ngens-N10000.csv"
pop.size <- 10000
coding.geno.size <- 500000
last.gen.to.calculate.with <- 10*pop.size


#______________________________________________________________________________________________________#

dfe.ests <- read.csv(paste(c("2epoch_AlphaOmega_", outfile.name, ".csv"), collapse=""), header=TRUE)

cols.to.agg <- c("lambda", "selected.divergence", "alpha", "omega_a")
avg.dfe.ests <- aggregate(dfe.ests[cols.to.agg], by=list(dfe.ests$group.ID), FUN=mean, na.rm=TRUE)
sd.dfe.ests <- aggregate(dfe.ests[cols.to.agg], by=list(dfe.ests $group.ID), FUN=sd, na.rm=TRUE)
ci95_dfe.ests <- (sd.dfe.ests[cols.to.agg]/sqrt(10)) * 2.26
names(ci95_dfe.ests) <- c("ci95_lambda", "ci95_selected.divergence", "ci95_alpha", "ci95_omega_a")
avg.dfe.ests <- cbind(avg.dfe.ests, ci95_dfe.ests)

# simulation 'true' alpha values
del.cols.to.keep <- c("file", "generation", "num.delet.muts.fixed", "num.neut.muts.fixed")
ben.cols.to.keep <- c("file", "generation", "num.delet.muts.fixed", "num.ben.muts.fixed", "num.neut.muts.fixed")

# something went wrong in calculating summ stats - use gen 29000
true.del <- read.csv(del.summ.stats)
true.del <- true.del[which(true.del$generation == last.gen.to.calculate.with), del.cols.to.keep]	# get the numbers at the last generation in order to subtract

true.ben <- read.csv(ben.summ.stats)
true.ben <- true.ben[which(true.ben$generation == last.gen.to.calculate.with), ben.cols.to.keep]	# get the numbers at the last generation in order to subtract

# remnant from when I had multiple things to rbind here, keep consistency with names though
del.true <- true.del
ben.true <- true.ben

del.true$file <- matrix(unlist(strsplit(as.character(del.true$file), split="_rep")), ncol=2, byrow=TRUE)[,1]
ben.true$file <- matrix(unlist(strsplit(as.character(ben.true$file), split="_rep")), ncol=2, byrow=TRUE)[,1]

del.cols.to.agg <- c("num.delet.muts.fixed", "num.neut.muts.fixed")
ben.cols.to.agg <- c("num.delet.muts.fixed", "num.ben.muts.fixed", "num.neut.muts.fixed")

avg.del.true <- aggregate(del.true[del.cols.to.agg], by=list(del.true$file), FUN=mean, na.rm=TRUE)
sd.del.true <- aggregate(del.true[del.cols.to.agg], by=list(del.true$file), FUN=sd, na.rm=TRUE)
ci95_avg.del.true <- (sd.del.true[del.cols.to.agg]/sqrt(10)) * 2.26
names(ci95_avg.del.true) <- c("ci95_num.delet.muts.fixed", "ci95_num.neut.muts.fixed")
avg.del.true <- cbind(avg.del.true, ci95_avg.del.true)

avg.ben.true <- aggregate(ben.true[ben.cols.to.agg], by=list(ben.true$file), FUN=mean, na.rm=TRUE)
sd.ben.true <- aggregate(ben.true[ben.cols.to.agg], by=list(ben.true$file), FUN=sd, na.rm=TRUE)
ci95_avg.ben.true <- (sd.ben.true[ben.cols.to.agg]/sqrt(10)) * 2.26
names(ci95_avg.ben.true) <- c("ci95_num.delet.muts.fixed", "ci95_num.ben.muts.fixed", "ci95_num.neut.muts.fixed")
avg.ben.true <- cbind(avg.ben.true, ci95_avg.ben.true)

avg.del.true$geno.size <- coding.geno.size
avg.del.true$pop.size <- pop.size
avg.ben.true$geno.size <- coding.geno.size
avg.ben.true$pop.size <- pop.size


avg.ben.true$dNdS.total <- ((avg.ben.true$num.delet.muts.fixed + avg.ben.true$num.ben.muts.fixed) / (avg.ben.true$geno.size * 0.75)) / (avg.ben.true$num.neut.muts.fixed / (avg.ben.true$geno.size * 0.25))
avg.del.true$dNdS.total <- (avg.del.true$num.delet.muts.fixed / (avg.del.true$geno.size * 0.75)) / (avg.del.true$num.neut.muts.fixed / (avg.del.true$geno.size * 0.25))

avg.ben.true$dNdS.delet <- (avg.ben.true$num.delet.muts.fixed / (avg.ben.true$geno.size * 0.75)) / (avg.ben.true$num.neut.muts.fixed / (avg.ben.true$geno.size * 0.25))
avg.del.true$dNdS.delet <- (avg.del.true$num.delet.muts.fixed / (avg.del.true$geno.size * 0.75)) / (avg.del.true$num.neut.muts.fixed / (avg.del.true$geno.size * 0.25))

avg.ben.true$dNdS.ben <- (avg.ben.true$num.ben.muts.fixed / (avg.ben.true$geno.size * 0.75)) / (avg.ben.true$num.neut.muts.fixed / (avg.ben.true$geno.size * 0.25))
avg.del.true$dNdS.ben <- (0 / (avg.del.true$geno.size * 0.75)) / (avg.del.true$num.neut.muts.fixed / (avg.del.true$geno.size * 0.25))

avg.ben.true$alpha.true <- avg.ben.true$dNdS.ben / avg.ben.true$dNdS.total
avg.del.true$alpha.true <- avg.del.true$dNdS.ben / avg.del.true$dNdS.total

#_______________________________________________________________________________________________#
avg.ben.true$dN.ben.all <- avg.ben.true$num.ben.muts.fixed / (avg.ben.true$num.delet.muts.fixed + avg.ben.true$num.ben.muts.fixed)
avg.del.true$dN.ben.all <- 0
#_______________________________________________________________________________________________#

# add blank colmn so can merge everything:
avg.del.true$num.ben.muts.fixed <- 0
avg.del.true$ci95_num.ben.muts.fixed <- 0

avg.true.values <- rbind(avg.ben.true , avg.del.true)
# remove "SampleOutput_" from front so matches for merge
avg.true.values$Group.1 <- matrix(unlist(strsplit(avg.true.values$Group.1, split="put_")), ncol=2, byrow=TRUE)[,2]

dat <- merge(avg.dfe.ests, avg.true.values)

simulation.list <- c(paste(unique(dat$Group.1)))
cols.list <- c("#a6cee3", "#ff7f00", "#1f78b4", "#b2df8a", "#e31a1c", "#33a02c", "#fb9a99", "#fdbf6f", "#cab2d6", "#006d2c", "#2ca25f", "#66c2a4", "#99d8c9", "#ccece6", "#edf8fb", "#810f7c", "#8856a7", "#8c96c6", "#9ebcda", "#bfd3e6", "#edf8fb", "#006837", "#31a354", "#78c679", "#addd8e", "#d9f0a3", "#ffffcc", "#0868ac", "#43a2ca", "#7bccc4", "#a8ddb5", "#ccebc5", "#f0f9e8", "#a50f15", "#de2d26", "#fb6a4a", "#fc9272", "#fcbba1", "#fee5d9")	


##_____________________ FIGURE _____________________##

pdf(paste(c("AlphaValues_", outfile.name, ".pdf"), collapse=""), width=length(simulation.list)/2, height=5)
par(mar=c(2,4,2,1))

plot(0, type="n", xlim=c(0.75,length(simulation.list)+0.5), ylim=c(-0.1,0.85), xlab="", ylab="Alpha", xaxt="n", main=outfile.name)
legend("topright", pch=c(45, rep(16,30)), col=c("black", cols.list), c("observed alpha", simulation.list), bg="white", ncol=1, cex=0.4)
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1) ,]
	errbar(add=TRUE, iterate, temp.dat$alpha, yplus=temp.dat$alpha+temp.dat$ci95_alpha, yminus=temp.dat$alpha-temp.dat$ci95_alpha, pch=20, errbar.col=cols.list[iterate], cex=0.5)
	points(iterate, temp.dat$alpha, col="gray10", bg=cols.list[iterate], pch=21, cex=1.5)
	iterate <- iterate + 1
}
# true alphas
iterate <- 1
for(i in simulation.list){
	temp.dat <- dat[grep(i, dat$Group.1) ,]
	points(iterate, temp.dat$alpha.true, col="black", pch="-", cex=1.5)
	iterate <- iterate + 1
}
dev.off()







dir <- "~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/Sep28_NewRuns_LargerGenome_compareBens/Outputs_Sep29"
setwd(paste(c(dir, "/Outfiles"), collapse=""))



fixed.files <- c(
"FixedOutput_Sep29_N10000_24mbp_del_outc_rep1.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep2.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep3.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep4.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep5.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep6.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep7.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep8.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep9.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_outc_rep10.txt",
"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep1.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep2.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep3.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep4.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep5.txt",	
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep6.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep7.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep8.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep9.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.90_rep10.txt",
"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep1.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep2.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep3.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep4.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep5.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep6.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep7.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep8.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep9.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_asex0.99_rep10.txt",
"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep1.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep2.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep3.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep4.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep5.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep6.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep7.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep8.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep9.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.90_rep10.txt",
"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep1.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep2.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep3.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep4.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep5.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep6.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep7.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep8.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep9.txt",
	"FixedOutput_Sep29_N10000_24mbp_del_self0.99_rep10.txt",
	
"FixedOutput_Sep29_N10000_26mbp_del_outc_rep1.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep2.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep3.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep4.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep5.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep6.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep7.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep8.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep9.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_outc_rep10.txt",
"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep1.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep2.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep3.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep4.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep5.txt",	
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep6.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep7.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep8.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep9.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.90_rep10.txt",
"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep1.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep2.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep3.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep4.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep5.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep6.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep7.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep8.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep9.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_asex0.99_rep10.txt",
"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep1.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep2.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep3.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep4.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep5.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep6.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep7.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep8.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep9.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.90_rep10.txt",
"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep1.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep2.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep3.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep4.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep5.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep6.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep7.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep8.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep9.txt",
	"FixedOutput_Sep29_N10000_26mbp_del_self0.99_rep10.txt",

fixed.files <- c(

"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep1.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep2.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep3.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep4.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep5.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep6.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep7.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep8.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep9.txt",
	"FixedOutput_Sep29_N10000_24mbp_ben-del_outc_rep10.txt",

"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep1.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep2.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep3.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep4.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep5.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep6.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep7.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep8.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep9.txt",
	"FixedOutput_Sep29_N10000_26mbp_ben-del_outc_rep10.txt"
)

# find the effect sizes of each mutation type
library(scales)




# make a mixture gamma distribution

#_________ gamma distribution 1 ________
gamma.mean1 <- 0.01
gamma.shape1 <- 0.3
gamma.scale1 <- gamma.mean1/gamma.shape1

gamma.hist1 <- qgamma(0.0001*c(0:10000), shape=gamma.shape1, scale=gamma.scale1)
## hist1 <- hist(gamma.hist1, col="steelblue4", breaks=250, xlab="Selection coefficient")


#_________ gamma distribution 2 ________
gamma.mean2 <- 0.5
gamma.shape2 <- 10
gamma.scale2 <- gamma.mean2/gamma.shape2

gamma.hist2 <- qgamma(0.0001*c(0:10000), shape=gamma.shape2, scale=gamma.scale2)
## hist2 <- hist(gamma.hist2, col="steelblue4", breaks=250, xlab="Selection coefficient")


#_________ combined gamma distribution with weighting ________




pdf(paste(c(dir,"/EffectSizes_FixedMutations.pdf"), collapse=""), width=8, height=6)
bins <- seq(-0.5,0.5, by=0.0025)
par(mar=c(4,4,1,1))

combined.gamma.hist <- (0.95*gamma.hist1) + (0.05*gamma.hist2)
hist(-combined.gamma.hist, xlim=c(-0.2,0.2), ylim=c(0,1500), col="red", breaks=bins, main="Modelled distributions", xlab="Selection coefficient")
hist(gamma.hist1, col="green3", breaks=bins, add=TRUE)


hist(0.001, xlim=c(-0.1,0.5), ylim=c(0,200), xlab="Selection coefficient", main="Observed fixed mutations (10 replicates overlaid)")
for(i in fixed.files){
	dat <- read.table(i, skip=2)
	names(dat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
	temp <- hist(dat$seln_coeff[dat$seln_coeff != 0], breaks=bins, plot=FALSE)
	cuts <- cut(temp$breaks, c(-Inf, 0, Inf))
	plot(temp, xlim=c(-0.05,0.4), ylim=c(0,800), col=c(alpha("red", 0.1), alpha("green3", 0.1))[cuts], add=TRUE)
}

dev.off()





dir <- "/cap1/kgilbert/NewSims_MesserLikeParams"



# BENEFICIALS:

## source for ben output analyses

source('/cap1/kgilbert/NewSims_MesserLikeParams/CalculateSummStats_wBens.R', chdir = TRUE)


# 1,000 gen bneck

setwd("/cap1/kgilbert/DFE_alpha/Inputs")

sample.files <- system("ls SampleOutput_*_ben*", intern=TRUE)
fixed.files <- system("ls FixedOutput_*_ben*", intern=TRUE)
full.files <- system("ls FullOutput_*_ben*", intern=TRUE)


pop.size <- 10000
coding.genome.size <- 500000
gens.to.sample.at <- format(seq(pop.size*1, pop.size*50, by=pop.size), scientific=FALSE)  # may need to change this according to sims analyzed 
num.gens.sampled <- length(gens.to.sample.at)
samp.size <- 100

outfile.name <- "April18_MesserLikeParams_50Ngens-N10000.csv"

summ.stats(
	sample.output.files=sample.files, 
	fixed.output.files=fixed.files, 
	full.output.files=full.files,
	summ.stats.output.file=paste(c(dir,"/SummStats_withBens_", outfile.name), collapse=""), 
	num.gens.sampled=num.gens.sampled, 
	num.inds.sampled=samp.size, 
	sequence.length=coding.genome.size, 
	pop.size=pop.size, 
	sub.sample.final=TRUE)
	
	


# DELETERIOUS:

source('/cap1/kgilbert/NewSims_MesserLikeParams/CalculateSummStats_Delets.R', chdir = TRUE)

setwd("/cap1/kgilbert/DFE_alpha/Inputs")


sample.files <- system("ls SampleOutput_*_del*", intern=TRUE)
fixed.files <- system("ls FixedOutput_*_del*", intern=TRUE)
full.files <- system("ls FullOutput_*_del*", intern=TRUE)


summ.stats(
	sample.output.files=sample.files, 
	fixed.output.files=fixed.files, 
	full.output.files=full.files,
        summ.stats.output.file=paste(c(dir,"/SummStats_deletsOnly_", outfile.name), collapse=""),
        num.gens.sampled=num.gens.sampled,
        num.inds.sampled=samp.size,
        sequence.length=coding.genome.size,
        pop.size=pop.size,
        sub.sample.final=TRUE)	
	

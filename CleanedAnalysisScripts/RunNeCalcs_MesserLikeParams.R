

setwd("/cap1/kgilbert/NewSims_MesserLikeParams")

source('/cap1/kgilbert/NewSims_MesserLikeParams/Calculate_RealizedNe_Messer.R', chdir = TRUE)


sample.files <- system("ls /cap1/kgilbert/DFE*/Inputs/SampleOutput_*", intern=TRUE)
fixed.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FixedOutput_*", intern=TRUE)
full.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FullOutput_*", intern=TRUE)


sequence.length=(500000)
pop.size <- 10000
num.inds.sampled=100
ne.gen.to.estimate <- pop.size*50     ## MAY NEED TO CHANGE THIS DEPENDING ON SIMS RUN

output.file <- "_N10000_EstimatedNe.csv"


est.ne(sample.output.files=sample.files, 
	full.output.files=full.files, 
	fixed.output.files=fixed.files, 
	summ.stats.output.file=output.file, 
	sample.size=num.inds.sampled, 
	sequence.length=sequence.length, 
	pop.size=pop.size,
	last.gen=ne.gen.to.estimate)

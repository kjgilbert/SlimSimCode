source('/cap1/kgilbert/NewSims_MesserLikeParams/LoadCalculationScripts/Calculate_RealizedVsRecessiveLoad.R', chdir = TRUE)
setwd("/cap1/kgilbert/NewSims_MesserLikeParams/LoadCalculationScripts")

fixed.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FixedOutput_*", intern=TRUE)
full.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FullOutput_*", intern=TRUE)
sample.files <- system("ls /cap1/kgilbert/DFE*/Inputs/SampleOutput_*", intern=TRUE)



pop.size <- 10000
samp.size <- 100

from <- c(100000)
until <- c(250000)
gen.samps <- format(c((pop.size*10)+1, seq((pop.size*10.01),(pop.size*10.2), by=100), seq((pop.size*10.3), (pop.size*10.5), by=1000), seq((pop.size*11), (pop.size*25), by=10000)), scientific=FALSE)
gen.samps <- gen.samps[-39] # take out the last one here because I was set up to just do the subsamps and a separate line in the code handles the full files at the very last gen
filename <- "RealizedVsRecessiveLoad_Apr21_TransitionsAndBnecks_N10000.csv"


do.load.calcs(fixed.files=fixed.files, sample.files=sample.files, full.files=full.files, load.stats.output.file=filename, from.gen=from, current.gen=until, inds.sampled=samp.size, pop.size=pop.size, gens.sampled=gen.samps)

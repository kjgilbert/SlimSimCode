source('/cap1/kgilbert/NewSims_MesserLikeParams/LoadCalculationScripts/Calculate_RealizedVsRecessiveLoad.R', chdir = TRUE)
setwd("/cap1/kgilbert/NewSims_MesserLikeParams/LoadCalculationScripts")

fixed.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FixedOutput_*", intern=TRUE)
full.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FullOutput_*", intern=TRUE)
sample.files <- system("ls /cap1/kgilbert/DFE*/Inputs/SampleOutput_*", intern=TRUE)


from <- c(40000)
until <- c(100000)
filename <- "RealizedVsRecessiveLoad_Nov19_N10000.csv"

pop.size <- 10000
samp.size <- 100

do.load.calcs(fixed.files=fixed.files, sample.files=sample.files, full.files=full.files, load.stats.output.file=filename, from.gen=from, current.gen=until, inds.sampled=samp.size, pop.size=pop.size)

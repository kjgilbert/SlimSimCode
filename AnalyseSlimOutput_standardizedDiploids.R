

# m1 mutations neutral -> synonymous
# m2 and m3 mutations deleterious -> nonsynonymous

# output files are  
#	"Sample..." which are the samples over all generations except the last, 
#		as well as any fixed mutation outputs (eventally want that to only be at the last generation)
#	"Full..." which is all pop data at the last generation (not fixed muts)


## MUST REMEMBER THESE SCRIPTS ARE ONLY FOR m1 NEUTRAL and m2, m3 DELETERIOUS


# calculations for Watterson's theta:
source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/CalculateTheta.R', chdir = TRUE)
## where dat is the genomes output with sep="A"
#	calc.theta(dat, num.inds.sampled, sequence.length)
#	return(theta)


# calculate pi, pi_n, pi_s, pi_n/pi_s
source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/CalculatePi.R', chdir = TRUE)
## where mut.id.dat is the poly muts output
## where names(mut.id.dat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
## where muts.occurring.dat is the genomes output with sep="A"
#	calc.pi.stats(mut.id.dat, muts.occurring.dat, num.inds.sampled, sequence.length)
#	return(c(pi, pi_n, pi_s))	


# calculate number of mutations (mean and variance) per individual - will use both sample and fixed mut data
source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/DeletMutsPerInd.R', chdir = TRUE)
## where poly.mut.dat is the poly muts output
## where names(poly.mut.dat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
## where ind.dat is the genomes output with sep="A"
## where fixed.mut.dat is the fixed mutations output at the last generation (i.e. for all mutations that ever fixed)
## where names(fixed.mut.dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
#	mean.var.muts(poly.mut.dat, ind.dat, generation, fixed.mut.dat, num.inds.sampled)
#		return(c(
#		mean.delet.muts.per.ind.poly, var.delet.muts.per.ind.poly, 
#		mean.neut.muts.per.ind.poly, var.neut.muts.per.ind.poly,
#		mean.total.muts.per.ind.poly, var.total.muts.per.ind.poly,
#		mean.delet.muts.per.ind.all, var.delet.muts.per.ind.all, 
#		mean.neut.muts.per.ind.all, var.neut.muts.per.ind.all,
#		mean.total.muts.per.ind.all, var.total.muts.per.ind.all,
#		num.delet.muts.fixed, num.neut.muts.fixed))


# get mean and variance in fitness at the last time point ( the full output )
source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/MeanFitness.R', chdir = TRUE)
## where full.poly.muts.dat is the poly muts output from the final and full population sample
## where names(poly.mut.dat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
## where full.genomes.dat is the full, final genomes output with sep="A"
## where fixed.mut.dat is the fixed mutations output at the last generation (i.e. for all mutations that ever fixed)
## where names(fixed.mut.dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
#	calc.fitness(full.poly.muts.dat, full.genomes.dat, fixed.mut.dat)
# 	return(c(mean.fitness.poly, var.fitness.poly, mean.fitness.total, var.fitness.total))



# __________
num.gens.sampled <- 10
num.inds.sampled <- 100
sequence.length <- 100000000 * (2/10)	# because coding regions are only 200 out of every 1000
pop.size <- 10000
# __________


setwd("~/Documents/My_Documents/UofToronto/SLiM/Outputs")

sample.output.files <- c(
	"SampleOutput_Sep19_N10000_del_self99_rep5.txt"
	)
	
full.output.files <- c(
	"FullOutput_Sep19_N10000_del_self99_rep5.txt"
	)
	
fixed.output.files <- c(
	"FixedOutput_Sep19_N10000_del_self99_rep5.txt"
)


results <- data.frame(matrix(NA, ncol=25))
names(results) <- c("file", "generation", "theta", "pi", "pi_n", "pi_s", "pi_n.pi_s", "mean.delet.muts.per.ind.poly", "var.delet.muts.per.ind.poly", "mean.neut.muts.per.ind.poly", "var.neut.muts.per.ind.poly", "mean.total.muts.per.ind.poly", "var.total.muts.per.ind.poly", "mean.delet.muts.per.ind.all", "var.delet.muts.per.ind.all", "mean.neut.muts.per.ind.all", "var.neut.muts.per.ind.all", "mean.total.muts.per.ind.all", "var.total.muts.per.ind.all", "num.delet.muts.fixed", "num.neut.muts.fixed", "mean.fitness.poly", "var.fitness.poly", "mean.fitness.total", "var.fitness.total")


iterate <- 1
for(i in 1:length(sample.output.files)){
	# go through each file
	
	sample.file <- sample.output.files[i]
	full.file <- full.output.files[i]
	fixed.file <- fixed.output.files[i]


	## sample data output
	poly.mut.id.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*(num.gens.sampled-1)), by=2)])
	generations <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n ' GS ' ", sample.file), collapse=""), intern=TRUE), split=" ")), ncol=4, byrow=TRUE)[,2])
	genomes.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Genomes ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*(num.gens.sampled-1)), by=2)])


	## full data output
	full.samp.muts.start <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", full.file), collapse=""), intern=TRUE), split=":"))[1])
	full.samp.inds.start <- as.numeric(strsplit(system(paste(c("grep -n Individuals ", full.file), collapse=""), intern=TRUE), split=":")[[1]][1])
	full.samp.genomes.start <- as.numeric(strsplit(system(paste(c("grep -n Genomes ", full.file), collapse=""), intern=TRUE), split=":")[[1]][1])
	full.samp.file.end <- as.numeric(head(tail(unlist(strsplit(system(paste(c("wc -l ", full.file), collapse=""), intern=TRUE), split=" ")), n=2), n=1))
	

	## fixed data output
	fixed.mut.id.start <- 2
	fixeddat <- read.table(fixed.output.files, skip=fixed.mut.id.start)
	names(fixeddat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")


	for(j in 1:num.gens.sampled){
		# go through each time point of a given file
		
		# for sample time points:
		if(j < num.gens.sampled){
			polydat <- read.table(sample.file, skip=poly.mut.id.starts[j], nrow=((genomes.starts[j]-1) - poly.mut.id.starts[j]), sep=" ")
			names(polydat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
	
			genodat <- read.table(sample.file, skip=genomes.starts[j], nrow=(num.inds.sampled * 2), sep="A")
	
			gen <- generations[j]

			theta <- calc.theta(dat=genodat, num.inds.sampled, sequence.length)
			pi.stats <- calc.pi.stats(mut.id.dat=polydat, muts.occurring.dat=genodat, num.inds.sampled, sequence.length)
			mut.stats <- mean.var.muts(poly.mut.dat=polydat, ind.dat=genodat, generation=gen, fixed.mut.dat=fixeddat, num.inds.sampled)
			fitness.stats <- calc.fitness(diploid.poly.muts.dat=polydat, full.genomes.dat=genodat, fixed.mut.dat=fixeddat, pop.size=num.inds.sampled)

			results[iterate ,] <- c(sample.file, gen, theta, pi.stats, pi.stats[2]/pi.stats[3], mut.stats, fitness.stats)
			iterate <- iterate + 1
		}
print(c("spot 1", i,j))		
		# for last, full time point
		if(j == num.gens.sampled){
			polydat <- read.table(full.file, skip=full.samp.muts.start, nrow=((full.samp.inds.start-1) - full.samp.muts.start), sep=" ")
			names(polydat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
	
			genodat <- read.table(full.file, skip=full.samp.genomes.start, nrow=(pop.size*2), sep="A")

			gen <- j * pop.size

			theta <- calc.theta(dat=genodat, num.inds.sampled, sequence.length)
			pi.stats <- calc.pi.stats(mut.id.dat=polydat, muts.occurring.dat=genodat, num.inds.sampled, sequence.length)
			mut.stats <- mean.var.muts(poly.mut.dat=polydat, ind.dat=genodat, generation=gen, fixed.mut.dat=fixeddat, num.inds.sampled=(2*pop.size))
			fitness.stats <- calc.fitness(diploid.poly.muts.dat=polydat, full.genomes.dat=genodat, fixed.mut.dat=fixeddat, pop.size=pop.size)
			
			results[iterate ,] <- c(full.file, gen, theta, pi.stats, pi.stats[2]/pi.stats[3], mut.stats, fitness.stats)
			iterate <- iterate + 1
		}
print(c("spot 2", i,j))		
	}
}

write.csv(results, file="SlimSummaryStatsResults.csv")
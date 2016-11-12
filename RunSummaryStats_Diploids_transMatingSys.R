

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
## where poly.dat is the poly muts output
## where names(poly.dat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
## where muts.occurring.dat is the genomes output with sep="A"
#	calc.pi.stats(poly.dat, muts.occurring.dat, num.inds.sampled, sequence.length)
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







summ.stats <- function(sample.output.files, fixed.output.files, summ.stats.output.file, num.gens.sampled, num.inds.sampled, sequence.length, pop.size, sub.sample.final=TRUE){
	
	results <- data.frame(matrix(nrow=0, ncol=27))
	names(results) <- c("ignore", "file", "generation", "theta", "theta.neut", "pi", "pi_n", "pi_s", "pi_n.pi_s", "mean.delet.muts.per.ind.poly", "var.delet.muts.per.ind.poly", "mean.neut.muts.per.ind.poly", "var.neut.muts.per.ind.poly", "mean.total.muts.per.ind.poly", "var.total.muts.per.ind.poly", "mean.delet.muts.per.ind.all", "var.delet.muts.per.ind.all", "mean.neut.muts.per.ind.all", "var.neut.muts.per.ind.all", "mean.total.muts.per.ind.all", "var.total.muts.per.ind.all", "num.delet.muts.fixed", "num.neut.muts.fixed", "mean.fitness.poly", "var.fitness.poly", "mean.fitness.total", "var.fitness.total")
	write.table(results, append=FALSE, file=summ.stats.output.file, sep=",", col.names=TRUE)
	
	
	iterate <- 1
	for(i in 1:length(sample.output.files)){
		# go through each file
		
		sample.file <- sample.output.files[i]
		fixed.file <- fixed.output.files[i]
	
	
		## sample data output
		poly.mut.id.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*num.gens.sampled), by=2)])
		generations <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n ' GS ' ", sample.file), collapse=""), intern=TRUE), split=" ")), ncol=4, byrow=TRUE)[,2])
		genomes.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Genomes ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*num.gens.sampled), by=2)])
	
	
		## fixed data output
		if(length(readLines(fixed.file)) == 2){	# then no mutations fixed
			fixeddat <- NULL
		}else{	# otherwise read in fixed mutations as normal
			fixed.mut.id.start <- 2
			fixeddat <- read.table(fixed.file, skip=fixed.mut.id.start)
			names(fixeddat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
		}
	
	
		for(j in 1:num.gens.sampled){
			# go through each time point of a given file
			
			# for sample time points:
			polydat <- read.table(sample.file, skip=poly.mut.id.starts[j], nrow=((genomes.starts[j]-1) - poly.mut.id.starts[j]), sep=" ")
			names(polydat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
	
			genodat <- read.table(sample.file, skip=genomes.starts[j], nrow=(num.inds.sampled * 2), sep="A")
	
			gen <- generations[j]

			theta <- calc.theta(genome.dat=genodat, poly.dat=polydat, fixed.dat=fixeddat, generation=gen, num.inds.sampled, sequence.length)
			pi.stats <- calc.pi.stats(poly.dat=polydat, genome.dat=genodat, fixed.dat=fixeddat, generation=gen, num.inds.sampled, genome.size=sequence.length, use.manual.sample=FALSE)
			mut.stats <- mean.var.muts(poly.mut.dat=polydat, genome.dat=genodat, generation=gen, fixed.mut.dat=fixeddat, num.inds.sampled)
			fitness.stats <- calc.fitness(diploid.poly.muts.dat=polydat, full.genomes.dat=genodat, fixed.mut.dat=fixeddat, pop.size=num.inds.sampled, generation=gen)

			temp.results <- c(sample.file, gen, theta, pi.stats, pi.stats[2]/pi.stats[3], mut.stats, fitness.stats)
			write.table(t(temp.results), append=TRUE, file=summ.stats.output.file, sep=",", col.names=FALSE)
			iterate <- iterate + 1
			
			# clear previous data, see if this solves the weird plot results
			polydat <- NULL
			genodat <- NULL		
		}
		fixeddat <- NULL
	}
	print("Calculations complete :) ")
}

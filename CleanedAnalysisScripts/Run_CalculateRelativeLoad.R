source('/cap1/kgilbert/NewSims_MesserLikeParams/LoadCalculationScripts/Calculate_MeanFitness_deletOnlyAndFull_MeanNumMuts_sRange.R', chdir = TRUE)

# April 21 transition mating systems with and without 1000 gen bottneck to N=1000

setwd("/cap1/kgilbert/NewSims_MesserLikeParams/LoadCalculationScripts")


output.file <- "RelativeLoadCalculated_Dec14_transitionAndBneck.csv"


fixed.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FixedOutput_*", intern=TRUE)
sample.files <- system("ls /cap1/kgilbert/DFE*/Inputs/SampleOutput_*", intern=TRUE)

sequence.length=(500000)
pop.size <- 10000
num.inds.sampled=100

gens.to.sample.at <- format(c((pop.size*10)+1, seq((pop.size*10.01),(pop.size*10.2), by=100), seq((pop.size*10.3), (pop.size*10.5), by=1000), seq((pop.size*11), (pop.size*25), by=10000)), scientific=FALSE)
num.gens.sampled <- length(gens.to.sample.at)


effect.size.min <- c(-0.00001, -0.0001, -0.001, -0.01, -0.1, -1, -1, -1, -1, -1, -1, -0.1, -0.01, -0.001, -0.0001)                # the MOST negative in the range wanted
effect.size.max <- c(0, 0, 0, 0, 0, 0, -0.1, -0.01, -0.001, -0.0001, -0.00001, -0.01, -0.001, -0.0001, -0.00001)



# get: 
#	total number fixed delet muts in that range
#	avg total number delet muts per ind in that range



iterate <- 1
for(i in 1:length(sample.files)){
	# go through each file
	sample.file <- sample.files[i]
	fixed.file <- fixed.files[i]

	## sample data output
	poly.mut.id.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*num.gens.sampled), by=2)])
	generations <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n ' GS ' ", sample.file), collapse=""), intern=TRUE), split=" ")), ncol=4, byrow=TRUE)[,2])
	genomes.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Genomes ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*num.gens.sampled), by=2)])
		

	## fixed data output
	if(length(readLines(fixed.file)) == 2){		# if no mutations fixed
		fixeddat <- NULL
	}else{										# otherwise read in fixed mutations as normal
		fixed.mut.id.start <- 2
		fixeddat <- read.table(fixed.file, skip=fixed.mut.id.start)
		names(fixeddat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
	}
	

	# go through each time point of a given file and do the calculations at that point in time
	for(j in 1:num.gens.sampled){
		# lists mutation identities that are segregating in the ENTIRE pop, not necessarily just in this sample
		polydat <- read.table(sample.file, skip=poly.mut.id.starts[j], nrow=((genomes.starts[j]-1) - poly.mut.id.starts[j]), sep=" ")
		names(polydat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")

		# lists all mutations per individual in a sample of individuals
		genodat <- read.table(sample.file, skip=genomes.starts[j], nrow=(num.inds.sampled * 2), sep="A")

		gen <- generations[j]

		fitness.and.mut.window.stats <- calc.fitness.window(poly.mut.dat=polydat, full.genomes.dat=genodat, fixed.mut.dat=fixeddat, pop.size=num.inds.sampled, generation=gen, min.s=effect.size.min, max.s=effect.size.max)
		
		temp.results <- c(sample.file, gen, fitness.and.mut.window.stats)
		names(temp.results)[1] <- "filename"
		names(temp.results)[2] <- "generation"
		
		if(i==1 & j==1){
			write.table(t(temp.results), append=FALSE, file=paste(c(dir, "/", output.file), collapse=""), sep=",", col.names=TRUE)
		}else{
			write.table(t(temp.results), append=TRUE, file=paste(c(dir, "/", output.file), collapse=""), sep=",", col.names=FALSE)
		}
		iterate <- iterate + 1
		
		# clear previous data, see if this solves the weird plot results
		polydat <- NULL
		genodat <- NULL		
	}
	fixeddat <- NULL
}





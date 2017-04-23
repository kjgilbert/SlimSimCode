setwd("/cap1/kgilbert/NewSims_MesserLikeParams/LoadCalculationScripts")
source('/cap1/kgilbert/NewSims_MesserLikeParams/LoadCalculationScripts/CalculateHomHetLoad_andNumNeuts.R', chdir = TRUE)


fixed.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FixedOutput_*", intern=TRUE)
##	full.files <- system("ls /cap1/kgilbert/DFE*/Inputs/FullOutput_*", intern=TRUE)
sample.files <- system("ls /cap1/kgilbert/DFE*/Inputs/SampleOutput_*", intern=TRUE)




sequence.length=(500000)
pop.size <- 10000
num.inds.sampled=100


gens.to.sample.at <- format(c((pop.size*10)+1, seq((pop.size*10.01),(pop.size*10.2), by=100), seq((pop.size*10.3), (pop.size*10.5), by=1000), seq((pop.size*11), (pop.size*25), by=10000)), scientific=FALSE)
num.gens.sampled <- length(gens.to.sample.at)


effect.size.min <- c(-0.00001, -0.0001, -0.001, -0.01, -0.1, -1, -1, -1, -1, -1, -1, -0.1, -0.01, -0.001, -0.0001)
effect.size.max <- c(0, 0, 0, 0, 0, 0, -0.1, -0.01, -0.001, -0.0001, -0.00001, -0.01, -0.001, -0.0001, -0.00001)



output.file <- "HomHet_Load_Apr21_TransitionsAndBnecks.csv"



iterate <- 1
for(i in 1:length(sample.files)){						# go through simulation replicates
	# go through each file
	sample.file <- sample.files[i]
	fixed.file <- fixed.files[i]

	## sample data output
	poly.mut.id.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*num.gens.sampled), by=2)])
	generations <- as.numeric(matrix(unlist(strsplit(system(paste(c("grep -n ' GS ' ", sample.file), collapse=""), intern=TRUE), split=" ")), ncol=4, byrow=TRUE)[,2])
	genomes.starts <- as.numeric(unlist(strsplit(system(paste(c("grep -n Genomes ", sample.file), collapse=""), intern=TRUE), split=":"))[seq(1, (2*num.gens.sampled), by=2)])
		

	## fixed data output
	if(length(readLines(fixed.file)) == 2){		# if no mutations fixed
		fixed.mut.dat <- NULL
	}else{										# otherwise read in fixed mutations as normal
		fixed.mut.id.start <- 2
		fixed.mut.dat <- read.table(fixed.file, skip=fixed.mut.id.start)
		names(fixed.mut.dat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
	}

	# go through each time point of a given file and do the calculations at that point in time
	for(j in 1:num.gens.sampled){						# go through generations within simulation replicates
		# lists mutation identities that are segregating in the ENTIRE pop, not necessarily just in this sample
		poly.mut.dat <- read.table(sample.file, skip=poly.mut.id.starts[j], nrow=((genomes.starts[j]-1) - poly.mut.id.starts[j]), sep=" ")
		names(poly.mut.dat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")

		# lists all mutations per individual in a sample of individuals
		full.genomes.dat <- read.table(sample.file, skip=genomes.starts[j], nrow=(num.inds.sampled * 2), sep="A")

		generation <- generations[j]

		# get load due to only heterozygous muts (in polydat)
		# get load due to only homozygous muts (from polymorphic data only)
		# get load due to only homozygous muts (from polymorphic AND fixed data)


		hom.het.results <- calc.hom.het.load.windows(poly.mut.dat=poly.mut.dat, full.genomes.dat=full.genomes.dat, fixed.mut.dat=fixed.mut.dat, pop.size=num.inds.sampled, generation, min.s=effect.size.min, max.s=effect.size.max)


		
		temp.results <- c(sample.file, generation, hom.het.results)
		names(temp.results)[1] <- "filename"
		names(temp.results)[2] <- "generation"
		
		if(i==1 & j==1){
			write.table(t(temp.results), append=FALSE, file=output.file, sep=",", col.names=TRUE)
		}else{
			write.table(t(temp.results), append=TRUE, file=output.file, sep=",", col.names=FALSE)
		}
		iterate <- iterate + 1
		
		# clear previous data, see if this solves the weird plot results
		poly.mut.dat <- NULL
		full.genomes.dat <- NULL		
	}
	fixed.mut.dat <- NULL
}


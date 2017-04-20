source('/cap1/kgilbert/RunningR_SlimResults/Jan4_N1000_SlimOutputsFiles/SourceFilesR/MeanFitness_deletOnlyAndFull_MeanNumMuts_sRange.R', chdir = TRUE)


# Dec 14 transition mating systems with 1000 gen bottneck to N=1000

dir <- "/cap1/kgilbert/RunningR_SlimResults/Dec14_comboMateSysBneck"
setwd(paste(c(dir, "/SlimOutputs"), collapse=""))

sample.files <- system("ls SampleOutput_Dec14_N10000_25mbp_*", intern=TRUE)
fixed.files <- system("ls FixedOutput_Dec14_N10000_25mbp_*", intern=TRUE)



sequence.length=(100000000 * (2.5/10))
pop.size <- 10000
num.inds.sampled=100

gens.to.sample.at <- format(c((pop.size*10)+1, seq((pop.size*10)+2, (pop.size*10.09)+2, by=150), 
	(pop.size*10.1)+1,
	seq((pop.size*10.1)+2, (pop.size*10.15)+2, by=150), 
	seq((pop.size*10.2)+2, (pop.size*10.5)+2, by=1000),
	seq((pop.size*11)+2, (pop.size*15)+2, by=10000)), scientific=FALSE)
num.gens.sampled <- length(gens.to.sample.at)


# could either do fixed bins, or quantiles (do both) - quantiles more useful since most muts will be nearer to neutrality and very few extra large effect

# make quantile bins based on the main (95%) distribution of the deleterious mutations
gamma.mean <- 0.01
gamma.shape <- 0.3
gamma.scale <- gamma.mean/gamma.shape
quantiles <- qgamma(0.025*c(0:25),shape=gamma.shape,scale=gamma.scale)
# remove the 0 and the inf
quantiles <- quantiles[quantiles > 0 & quantiles < Inf]
quantiles <- -quantiles

effect.size.min <- c(-0.00001, -0.0001, -0.001, -0.01, -0.1, -1, -1, -1, -1, -1, -1, -0.1, -0.01, -0.001, -0.0001)                # the MOST negative in the range wanted
effect.size.max <- c(0, 0, 0, 0, 0, 0, -0.1, -0.01, -0.001, -0.0001, -0.00001, -0.01, -0.001, -0.0001, -0.00001)
	# the LEAST negative in the range wanted - will do everything from 0 (neutrals are automatically excluded in the calculation code)

output.file <- "RelativeLoadCalculated_Dec14_transitionAndBneck.csv"

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





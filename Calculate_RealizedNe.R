setwd("/cap1/kgilbert/WestgridOutputs")

# calculate pi:
source('/cap1/kgilbert/WestgridOutputs/CalculatePi.R', chdir = TRUE)
	# gives back: pi overall, pi_n, pi_s


##	sample.output.files <- system("ls SampleOutput_Nov19_N10000_25mbp_*", intern=TRUE)
##	full.output.files <- system("ls FullOutput_Nov19_N10000_25mbp_*", intern=TRUE)
##	fixed.output.files <- system("ls FixedOutput_Nov19_N10000_25mbp_*", intern=TRUE)
##	sequence.length <- 25000000
##	pop.size <- 10000
##	sample.size <- 100
##	summ.stats.output.file <- "test"



est.ne <- function(sample.output.files, full.output.files, fixed.output.files, summ.stats.output.file, sample.size, sequence.length, pop.size){
	
	results <- data.frame(matrix(nrow=0, ncol=4))
	names(results) <- c("ignore", "file", "generation", "Ne.neutral")
	write.table(results, append=FALSE, file=summ.stats.output.file, sep=",", col.names=TRUE)
	
	
	iterate <- 1
	for(i in 1:length(sample.output.files)){
		# go through each file
		
		sample.file <- sample.output.files[i]
		full.file <- full.output.files[i]
		fixed.file <- fixed.output.files[i]
	
		## full data output
		full.samp.muts.start <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", full.file), collapse=""), intern=TRUE), split=":"))[1])
		full.samp.inds.start <- as.numeric(strsplit(system(paste(c("grep -n Individuals ", full.file), collapse=""), intern=TRUE), split=":")[[1]][1])
		full.samp.genomes.start <- as.numeric(strsplit(system(paste(c("grep -n Genomes ", full.file), collapse=""), intern=TRUE), split=":")[[1]][1])
		full.samp.file.end <- as.numeric(head(tail(unlist(strsplit(system(paste(c("wc -l ", full.file), collapse=""), intern=TRUE), split=" ")), n=2), n=1))
		
	
		## fixed data output
		if(length(readLines(fixed.file)) == 2){	# then no mutations fixed
			fixeddat <- NULL
		}else{	# otherwise read in fixed mutations as normal
			fixed.mut.id.start <- 2
			fixeddat <- read.table(fixed.file, skip=fixed.mut.id.start)
			names(fixeddat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
		}
	

		# for last, full time point
		polydat <- read.table(full.file, skip=full.samp.muts.start, nrow=((full.samp.inds.start-1) - full.samp.muts.start), sep=" ")
		names(polydat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")

		genodat <- read.table(full.file, skip=full.samp.genomes.start, nrow=(pop.size*2), sep="A")

		
		# sample from a vector of odd numbers since all inds have 2 paired genomes (diploid) and they start on an odd line and end on an even line
		odd.nums <- seq(1, (pop.size * 2), by=2)
		sub.samp <- sample(odd.nums, size=sample.size, replace=FALSE)
		diploid.sub.samp <- sort(c(sub.samp, (sub.samp + 1)))
		
		genodat <- genodat[diploid.sub.samp ,]
		
		gen <- 10*pop.size	# doing last gen, 10N
		
		pi.stats <- calc.pi.stats(poly.dat=polydat, genome.dat=genodat, fixed.dat=fixeddat, generation=gen, num.inds.sampled=sample.size, genome.size=sequence.length, use.manual.sample=TRUE)
		names(pi.stats) <- c("pi", "pi_n", "pi_s")
		
		total.mut.rate <- 7*10^-9
		neut.mut.rate <- 0.25*total.mut.rate
		seln.mut.rate <- 0.75*total.mut.rate
		
		
		# pi_s = 4 Ne mu_s		# only want Ne from neutral data, so use synonymous
		Ne_s <- pi.stats["pi_s"]/(4*neut.mut.rate)
##		Ne_n <- pi.stats["pi_n"]/(4*seln.mut.rate)
##		Ne_total <- pi.stats["pi"]/(4*total.mut.rate)
	
		
		temp.results <- c(sample.file, gen, Ne_s)
		write.table(t(temp.results), append=TRUE, file=summ.stats.output.file, sep=",", col.names=FALSE)
		iterate <- iterate + 1
			
			# clear previous data, see if this solves the weird plot results
		polydat <- NULL
		genodat <- NULL		
		fixeddat <- NULL
	}
}


##				# Ne*s bins:
##				#	0-1
##				#	1-10
##				#	10-100
##				#	100-inf
##				
##				# therefore s bins:
##				
##				zero.s.boundary <- 0
##				one.s.boundary <- 1/Ne_s
##				ten.s.boundary <- 10/Ne_s
##				hundred.s.boundary <- 100/Ne_s
##				inf.s.boundary <- 1
##				
##				# how many muts in each bin:
##				
##					# FIXED
##					
##					# POLY
##				
##				
##		
##		
## GLEMIN:
##		geno.wide.U <- neut.mut.rate*sequence.length
##		sd <- -0.01
##		hd <- 0.3
##		recomb <- 5*10^-8
##		
##		alpha <- e^-()
##		
##		
##		pop.size <- 10000
##		
##		self.rate <- 0.99
##		Fval <- self.rate/(2-self.rate)
##		
##		Ne <- (alpha*pop.size)/(1+Fval)
		
		

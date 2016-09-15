## CALCULATE MEAN AND VARIANCE FOR FITNESS

# for one population, within that pop, may later extend to do multiple pops


# NEED DIPLOID DATA
# so use full output at end



# this is a full output file, that was exported to the file name of my choosing, so it has the sections that can be subset, but does not need trimming of other outputs or run info

file <- "~/Documents/My_Documents/UofToronto/SLiM/Outputs/Old_Outputs/Test_SmallerN/out_Aug16_N5000_deletOnly_partRec_outcr_gen8N-full.txt"

# find mutation IDs and their effects
## mut.id.start <- "Mutations:"
mut.id.start <- as.numeric(strsplit(system(paste(c("grep -n Mutations ", file), collapse=""), intern=TRUE), split=":")[[1]][1])

## inds.start <- "Individuals:"
inds.start <- as.numeric(strsplit(system(paste(c("grep -n Individuals ", file), collapse=""), intern=TRUE), split=":")[[1]][1])

## genomes.start <- "Genomes:"
genomes.start <- as.numeric(strsplit(system(paste(c("grep -n Genomes ", file), collapse=""), intern=TRUE), split=":")[[1]][1])

file.end <- as.numeric(strsplit(system(paste(c("wc -l ", file), collapse=""), intern=TRUE), split=" ")[[1]][1])


#read all lines of the data file in
## dat <- readLines(file)



#take only the lines starting 1 line after "Mutations:" and up to 1 line before "Individuals"
## mutdat <- dat[(which(dat == mut.id.start)+1):(which(dat == inds.start)-1)]
mutdat <- read.table(file, skip=mut.id.start, nrow=((inds.start-1) - mut.id.start), sep=" ")
names(mutdat) <- c("mut.ID", "mut.type", "base.pos", "seln.coef", "dom.coef", "pop.ID", "gen.ID", "mut.prev")
# columns =
#	(1) a unique identifying ID, 
#	(2) the mutation type, 
#	(3) the base position, 
#	(4) the selection coefficient (here always 0 since this is a neutral model), 
#	(5) the dominance coefficient (here always 0.5), 
#	(6) the identifier of the subpopulation in which the mutation first arose, 
#	(7) the generation in which it arose, and 
#	(8) the prevalence of the mutation (the number of genomes that contain the mutation, where – in the way that SLiM uses the term “genome” – there are two genomes per individual)


#take only the lines starting 1 line after "Individuals:" and up to 1 line before "Genomes"
## inddat <- dat[(which(dat == inds.start)+1):(which(dat == genomes.start)-1)]
#		inddat <- read.table(file, skip=inds.start, nrow=((genomes.start-1) - inds.start), sep=" ")
#		names(inddat) <- c("pop.ind", "sex", "genome.1", "genome.2")
# don't really need this info, it's redundant


#take only the lines starting 1 line after "Genomes:" and up to end of file
## gendat <- dat[(which(dat == mut.id.start)+1):length(dat)]
gendat <- read.table(file, skip=genomes.start, sep="A")
names(gendat) <- c("pop.genome", "muts")


# loop through every individual, i.e. every pair
# get its fitness based on the diploid combination of mutations it has
pop.size <- 5000

fitness.results <- data.frame(matrix(NA, ncol=2))
names(fitness.results) <- c("individual", "fitness")

iterate.inds <- 1
for(i in seq(1, (pop.size*2), by=2)){
	genome.1.muts <- unlist(strsplit(as.character(gendat$muts[i]), split=" "))[-1]
	genome.2.muts <- unlist(strsplit(as.character(gendat$muts[(i+1)]), split=" "))[-1]

	# mutations both chromosomes have make them homozygotes for that mutation
	hom.muts <- which(genome.1.muts %in% genome.2.muts)
	# mutations only one chromosome has make them heterozygotes for that mutation
	het.muts <- c(which(!(genome.1.muts %in% genome.2.muts)), which(!(genome.2.muts %in% genome.1.muts)))
	# no need/no way to detect wild type homozygots, so don't need to worry about for fitness
	
	ind.fitness <- 1
	# match up homozygous mutations to fitness effects - homozygous = s
	for(j in hom.muts){
		ind.fitness <- ind.fitness + mutdat$seln.coef[mutdat$mut.ID == j]
	}
	# match up heterozygous mutations to fitness effects - heterozygous = hs
	for(j in het.muts){
		ind.fitness <- ind.fitness + (mutdat$seln.coef[mutdat$mut.ID == j] * mutdat$dom.coef[mutdat$mut.ID == j])
	}
	
	fitness.results[iterate.inds,] <- c(iterate.inds, ind.fitness)
	iterate.inds <- iterate.inds + 1
}



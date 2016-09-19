## CALCULATE THE MEAN AND VARIANCE IN NUMBER OF DELETERIOUS MUTATIONS PER INDIVIDUAL



# do from each a file of sampled inds and their muts, then tack on data from fixed muts because that gives time mutation took to fix


mean.var.muts <- function(poly.mut.dat, ind.dat, generation, fixed.mut.dat, num.inds.sampled){
	
	# from sample data at a single generation and fixed data from one file for all generations
	
	names(pdat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
	
	neut.muts <- c(which(pdat$mut.type == "m1"))
	delet.mut.IDs <- pdat$mut.ID[-neut.muts]


	
	return(poly.muts.mean, poly.muts.var, num.fixed.muts, all.muts.mean, all.muts.var)
}


# first go through poly.muts file at that time point and find the ID numbers for all the deleterious mutations

pdat <- read.table("~/Documents/My_Documents/UofToronto/SLiM/Outputs/N_1thousand_Outputs/output_Aug29_N1000_outc_samp100/poly.muts.out.90000out_Aug29_N1000_outc_samp100.txt")
# columns =
#	(1) a unique identifying ID, 
#	(2) the mutation type, 
#	(3) the base position, 
#	(4) the selection coefficient (here always 0 since this is a neutral model), 
#	(5) the dominance coefficient (here always 0.5), 
#	(6) the identifier of the subpopulation in which the mutation first arose, 
#	(7) the generation in which it arose, and 
#	(8) the prevalence of the mutation (the number of genomes that contain the mutation, where – in the way that SLiM uses the term “genome” – there are two genomes per individual)
names(pdat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")

neut.muts <- c(which(pdat$mut.type == "m1"))
delet.mut.IDs <- pdat$mut.ID[-neut.muts]



# then go through the genomes file and count how many of those mutations at that time point are in a given ind, do for all inds, then get mean and var

gdat <- read.table("~/Documents/My_Documents/UofToronto/SLiM/Outputs/N_1thousand_Outputs/output_Aug29_N1000_outc_samp100/genomes.out.90000out_Aug29_N1000_outc_samp100.txt", sep="A")
# haploid inds in rows- pop ID, ind ID, A for autosome, mutation IDs possessed by that ind

genomes <- gdat[,2]
delet.muts.per.ind <- data.frame(matrix(NA, ncol=2))
names(delet.muts.per.ind) <- c("ind.ID", "num.delet.muts")
num.inds.sampled <- 100

for(i in 1: num.inds.sampled){
	num.delet.muts.ind <- table(as.integer(unlist(strsplit(as.character(genomes[i]), split=" "))) %in% delet.mut.IDs)[2]	# take only TRUEs
	delet.muts.per.ind[i,] <- c(i, num.delet.muts.ind)
}

mean.num.per.ind <- mean(as.numeric(delet.muts.per.ind$num.delet.muts))
var.num.per.ind <- var(as.numeric(delet.muts.per.ind$num.delet.muts))







##		dat <- read.table("poly.muts.out.5000TestOutputs_Sample100.txt", sep=" ")
# columns =
#	(1) a unique identifying ID, 
#	(2) the mutation type, 
#	(3) the base position, 
#	(4) the selection coefficient (here always 0 since this is a neutral model), 
#	(5) the dominance coefficient (here always 0.5), 
#	(6) the identifier of the subpopulation in which the mutation first arose, 
#	(7) the generation in which it arose, and 
#	(8) how many generations it took for the mutation to fix
fdat <- read.table("~/Documents/My_Documents/UofToronto/SLiM/Outputs/N_1thousand_Outputs/output_Aug29_N1000_outc_samp100/fixed.muts.out.90000out_Aug29_N1000_outc_samp100.txt")



# do anything with fixed mutations?



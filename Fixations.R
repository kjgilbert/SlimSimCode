## HOW MANY FIXATIONS ARE SYNONYMOUS VS NONSYNONYMOUS



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

names(fdat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "time.to.fix")


count.fixations <- function(dat, neut.IDs, seln.IDs){
	# give a "dat" dataframe that MUST have at least 2 cols, one with "mut.ID" and one with "mut.type"
	# neut.IDs are the labels for neutral (synonymous mutations) in quotations - usu. "m1" for me ("m2" is introns)
	summ <- summary(dat$mut.type)
	num.neutral.fixed.muts <- summ[c(neut.IDs)]
	num.neutral.fixed.muts <- sum(num.neutral.fixed.muts)
	
	num.seln.fixed.muts <- summ[c(seln.IDs)]
	num.seln.fixed.muts <- sum(num.seln.fixed.muts)
	
	answer <- data.frame(matrix(c(num.seln.fixed.muts, num.neutral.fixed.muts), ncol=2))
	names(answer) <- c("num.selected.muts_fixed", "num.neutral.muts_fixed")
	return(answer)
}

test <- count.fixations(fdat, neut.IDs="m1", seln.IDs=c("m3","m4","m5","m6","m7","m8","m9","m10","m11","m12","m13"))
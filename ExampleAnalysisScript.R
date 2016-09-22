source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/RunSummaryStats_Diploids.R', chdir = TRUE)




setwd("~/Documents/My_Documents/UofToronto/SLiM/Outputs")

sample.files <- c(
	"SampleOutput_Sep19_N10000_del_self0.99_rep5.txt"
	)
	
full.files <- c(
	"FullOutput_Sep19_N10000_del_self0.99_rep5.txt"
	)
	
fixed.files <- c(
	"FixedOutput_Sep19_N10000_del_self0.99_rep5.txt"
)

# run the calculations:

summ.stats(
	sample.output.files=sample.files, 
	full.output.files=full.files, 
	fixed.output.files=fixed.files, 
	summ.stats.output.file="SummStats_filename.csv", 
	num.gens.sampled=10, 
	num.inds.sampled=100, 
	sequence.length=(100000000 * (2/10)), 
	pop.size=10000, 
	sub.sample.final=TRUE)
setwd("~/Documents/My_Documents/UofToronto/SLiM/Outputs/Asex_N10000outputs")


# synonymous muts -> all the neutral ones in exons
# nonsynonymous muts -> all the selected ones

# m2 mutations are in introns, so ignore
# use m1 and all m3+ mutation types


# harmonic sum function for theta
sum.a_k <- function(k){
	harmsum <- 0 
	for(i in 1:(k-1)){
		harmsum <- harmsum + 1/i
		if (harmsum >= (k-1)) break
	}
	return(harmsum)
}




# *****
num.gens.sampled <- 3
num.inds.sampled <- 100
sequence.length <- 100000000 * (2/10)	# because coding regions are only 200 out of every 1000
# *****
temp.results <- data.frame(matrix(NA, ncol=7))
names(temp.results) <- c("file", "generation", "pi", "pi_n", "pi_s", "pi_n.pi_s", "theta")


files <- c(
	"out_Aug16_N10000_asex_sample100_10Ngens.txt",
	"output_Aug29_N10000_asex_samp100_10Ngens.txt",
	"N10000_ByGens_Asexual/output_Sep13_N10000_asex_samp100_rep1_2Ngens.txt",
	"N10000_ByGens_Asexual/output_Sep13_N10000_asex_samp100_rep1_3Ngens.txt"
)

# if files are all in the directory, then instead:
##	files <- do a system command to list all in the directory? ls out*

# go through each file and get the summary stats per generation

for(i in files){
	# when coming from a stdout output file, i.e. one that is created from slim outside of its execution and script via ">", then run these lines to clean that file up:
	# cut start of file off and take only things after the line "1 " which comes just after it says things about starting
	system(paste(c("sed '1,/^1 $/d' ", i, " > clean_", i), collapse=""))

	# check in case have to remove the last line that says job completed:
	##	system(paste(c("sed '$d' clean_", i, " > ", i), collapse=""))
	
	system(paste(c("grep -n '#OUT\\|Mutations:\\|Individuals:\\|Genomes' clean_", i, " > temp.txt"), collapse=""))
	system("sed -i -- 's/#//g' temp.txt")					# remove '#' so R doesn't read it in as a comment
	system(paste(c("wc -l clean_", i, " > num.lines.txt"), collapse=""))	# how long is the file (need the last line number for the end)
	num.lines <- unlist(read.table("num.lines.txt")[1])
	
	temp <- read.table("temp.txt", sep=":")
	# before trimming the file, record the generations we sampled
	gens.sampled <- matrix(unlist(strsplit(as.character(temp$V3[temp$V2 == "OUT"]), split=" ")), ncol=num.gens.sampled)[2,]
	
	# then take just the first colum which is the line numbers we want:
	temp <- temp[,1]

	# loop through time points sampled and get summary stats:
	j <- 2
	for(gen in 1:num.gens.sampled){

		#<><><><><><><><><><>#                             #<><><><><><><><><><>#
		#<><><><><><><><><><>#    use POLY MUTS output     #<><><><><><><><><><>#
		#<><><><><><><><><><>#                             #<><><><><><><><><><>#

		# mutation file = mutation line +1 to genome line - 1
		mutfile.lines <- c((temp[j]+1),(temp[j+1]-1))
		system(paste(c("awk 'NR >= ", mutfile.lines[1], " && NR <= ", mutfile.lines[2], "' clean_", i, " > poly.muts.out.", gens.sampled[gen], i, sep=""), collapse=""))
		
		# make data frames of all possible neutral and deleterious mutations at all time points recorded
		neut.muts <- NULL
		seln.muts <- NULL
		
		dat <- read.table(paste(c("poly.muts.out.", gens.sampled[gen], i), collapse=""))
		names(dat) <- c("mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
		# remove intron neutral muts
		exclude.dat <- dat[dat$mut.type == "m2" ,]
		dat <- dat[dat$mut.type != "m2" ,]
		neut.muts <- dat[dat$mut.type == "m1" ,]
		seln.muts <- dat[dat$mut.type != "m1" ,]
		exclude.neuts <- exclude.dat$mut.ID
				
		# all possible mut IDs of any type (need for later calcs):
		
		selected.mut.IDs <- c(seln.muts$mut.ID)
		neutral.mut.IDs <- c(neut.muts$mut.ID)
		exclude.neutral.introns <- c(exclude.neuts)
		
		
		#<><><><><><><><><><>#                             #<><><><><><><><><><>#
		#<><><><><><><><><><>#    use GENOMES output       #<><><><><><><><><><>#
		#<><><><><><><><><><>#                             #<><><><><><><><><><>#
	
		# genome file = genome line +1 to OUT line - 1
		genomefile.lines <- c((temp[j+1]+1),(temp[j+2]-1))
		system(paste(c("awk 'NR >= ", genomefile.lines[1], " && NR <= ", genomefile.lines[2], "' clean_", i, " > genomes.out.", gens.sampled[gen], i, sep=""), collapse=""))
	
		generation <- gens.sampled[gen]
		dat <- read.table(paste(c("genomes.out.", gens.sampled[gen], i), collapse=""), sep="A")
	
		#	column 1 = pop ID, then colon, then the ID of the individual, then A for Autosome and then the number of mutations with the identifiers 0 through ...
		
		# so take random pairs, find the number of non-overlapping mutations, i.e. unique, so polymorphic between them
		# sum across all pairs
		# divide by num genomes * length genome
		
		# get all possible line number pairs
		pairs <- combn(1:length(dat[,1]), 2)
		total.poly <- 0
		nonsyn.total.poly <- 0
		syn.total.poly <- 0
		for(k in 1:dim(pairs)[2]){
			# take corresponding row for the pairs, and take second column of data file which has the mutation info
			pair1 <- dat[pairs[1,k], 2]
			pair2 <- dat[pairs[2,k], 2]
			
			total.muts1 <- unlist(strsplit(as.character(pair1), split=" "))[-1]	# must remove 1 b/c splitting the data makes a leading space that comes up as a shared mutation
			total.muts2 <- unlist(strsplit(as.character(pair2), split=" "))[-1]	# must remove 1 b/c splitting the data makes a leading space that comes up as a shared mutation
			total.non.intron.muts1 <- as.numeric(total.muts1[!(total.muts1 %in% exclude.neutral.introns)])
			total.non.intron.muts2 <- as.numeric(total.muts2[!(total.muts2 %in% exclude.neutral.introns)])
					
			neut.muts1 <- as.numeric(total.muts1[as.numeric(total.muts1) %in% neutral.mut.IDs])
			neut.muts2 <- as.numeric(total.muts2[as.numeric(total.muts2) %in% neutral.mut.IDs])
			seln.muts1 <- as.numeric(total.muts1[as.numeric(total.muts1) %in% selected.mut.IDs])
			seln.muts2 <- as.numeric(total.muts2[as.numeric(total.muts2) %in% selected.mut.IDs])
			
			# overall for pi
			total.matches <- length(intersect(total.non.intron.muts1, total.non.intron.muts2))	# how many mutations are the same
			poly1 <- length(total.non.intron.muts1) - total.matches		# how many muts in ind 1 are not in ind 2
			poly2 <- length(total.non.intron.muts2) - total.matches		# how many muts in ind 2 are not in ind 1
			temp.total.poly <- poly1 + poly2			# how many muts total are nonoverlapping between the two, i.e. the polymorphic/segregating sites
			total.poly <- total.poly + temp.total.poly
	
			# nonsynonymous (selected) for pi_n
			nonsyn.total.matches <- length(intersect(seln.muts1, seln.muts2))
			nonsyn.poly1 <- length(seln.muts1) - nonsyn.total.matches
			nonsyn.poly2 <- length(seln.muts2) - nonsyn.total.matches
			nonsyn.temp.total.poly <- nonsyn.poly1 + nonsyn.poly2
			nonsyn.total.poly <- nonsyn.total.poly + nonsyn.temp.total.poly
	
			# synonymous (neutral) for pi_s
			syn.total.matches <- length(intersect(neut.muts1, neut.muts2))
			syn.poly1 <- length(neut.muts1) - syn.total.matches
			syn.poly2 <- length(neut.muts2) - syn.total.matches
			syn.temp.total.poly <- syn.poly1 + syn.poly2
			syn.total.poly <- syn.total.poly + syn.temp.total.poly
		}
		pi <- total.poly / (num.inds.sampled * sequence.length)
		pi_n <- nonsyn.total.poly / (num.inds.sampled * sequence.length)
		pi_s <- syn.total.poly / (num.inds.sampled * sequence.length)
		
		
		## CALCULATE THETA
		
		# get total num poly sites across all samples = S_k; i.e., how many muts are not fixed in the sample?
		# first pair shared mutations:
		matched.sites <- intersect(unlist(strsplit(as.character(dat[1,2]), split=" ")), unlist(strsplit(as.character(dat[2,2]), split=" ")))
		# those shared mutations with 3rd individual and onward:
		for(k in 3:num.inds.sampled){
			temp.matched.sites <- intersect(unlist(strsplit(as.character(dat[k,2]), split=" ")), matched.sites)
			matched.sites <- temp.matched.sites
		}
		num.matched.sites <- length(matched.sites) - 1	# must subtract 1 because when I split the data into the second part, there is a leading space, so when it's split again, the space comes up as a shared mutation
		# that many ARE fixed in the SAMPLE, so total number muts in sample minus this should be number polymorphic sites in the sample
		# to get S_k, the total number of polymorphic sites, count all mutations and subtract the number matched across all
		all.muts <- 0
		for(k in 1:num.inds.sampled){
			temp.all.muts <- unlist(strsplit(as.character(dat[k,2]), split=" "))
			all.muts <- c(all.muts, temp.all.muts)
		}
		num.all.muts <- length(unique(all.muts)) - 1
			# must remove 1 b/c splitting the data makes a leading space that comes up as a shared mutation
		num.poly.muts <- num.all.muts - num.matched.sites
			
		theta <- num.poly.muts / (sum.a_k(num.inds.sampled) * sequence.length)
				
		# output:
		temp.results[gen,] <- c(i, generation, pi, pi_n, pi_s, pi_n/pi_s, theta)
	
		j <- j+5
		system(paste("rm poly.muts.out.*"))
		system(paste("rm genomes.out.*"))
	}	

	system(paste("rm num.lines.*"))
	system(paste("rm temp.*"))
	system(paste("rm clean_*"))
}

write.csv(temp.results, file=paste("popgen_paramResults_asex_Sep14.csv")













# mean and variance in number of delet muts per ind

# mean and variance in fitness per ind


## ASK BEN

# ARE THE FIXED MUTS THAT SLIM OUTPUTS SEPARATELY ENTIRELY DIFFERENT ID NUMBERS? i.e. the 'Genomes' output does not include any ID numbers that correspond to the fixed mut output?
# and does FIXED include LOST, or only gone to fixation? (should be the latter!!)
# does stopping slim from removing fixed delet muts change that?
# why are they removed? does it mean ONLY once removed another mut might be able to happen at that site?






##		results <- data.frame(matrix(NA, ncol=16))
##		names(results) <- c("file_generation", "pi", "pi_n", "pi_s", "pi_n.pi_s", "theta", "mean_num_delet_muts", "var_num_delet_muts", "mean_fitness", "var_fitness", "mean_fixed_delet_muts", "var_fixed_delet_muts", "mean_fixed_neut_muts", "var_fixed_neut_muts", "mean_fixed_ben_muts", "var_fixed_ben_muts")






## NEED DIPLOID INFO TO RECREATE FITNESS VALUES - b/c of dominance

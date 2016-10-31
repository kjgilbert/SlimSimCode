


# RESET WD EACH TIME TO SUIT!!!!


setwd("/cap1/kgilbert/DFE_alpha/Oct11_2epoch_DFE_Results")

# setwd("~/Documents/My_Documents/UofToronto/SLiM/Running_SLiM/Sep28_NewRuns_LargerGenome_compareBens/Outputs_Sep29/test")


# list only the directories then go through their contents
class0folders <- system("ls -d OutputClass0*", intern=TRUE)
class1folders <- system("ls -d OutputClass1*", intern=TRUE)

results <- data.frame(matrix(nrow=length(class0folders), ncol=20))
names(results) <- c("group.ID", "pop.size", "genome.size", "mut.types", "mating.system", "replicate", "N1.0", "N2.0", "t2.0", "Nw.0", "f0.0", "L.0", "N1.1", "N2.1", "t2.1", "Nw.1", "b", "Es", "f0.1", "L.1")

# there should be a matching class 1 to every class 0 folder
for(i in 1:length(class0folders)){
	contents0 <- system(paste(c("cat ", class0folders[i], "/est_dfe.out"), collapse=""), intern=TRUE)
	contents1 <- system(paste(c("cat ", class1folders[i], "/est_dfe.out"), collapse=""), intern=TRUE)
	rep.ID <- unlist(strsplit(class0folders[i], split="SFS_"))[2]
	group.ID <- unlist(strsplit(rep.ID, split="_rep"))[1]
	spec.IDs <- c(unlist(strsplit(rep.ID, split="_"))[3:7])
	row <- c(group.ID, spec.IDs, unlist(strsplit(contents0, split=" "))[c(2,4,6,8,10,12)], unlist(strsplit(contents1, split=" "))[c(2,4,6,8,10,12,14,16)])
	results[i,] <- row
}

write.table(results, file="est_dfe_results_2epoch_Oct11.csv", sep=",", col.names=TRUE)



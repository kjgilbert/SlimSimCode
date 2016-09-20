source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/MakePBSfiles_WestgridSlim.R', chdir = TRUE)
source('~/Documents/My_Documents/UofToronto/SLiM/SlimSimCode/MakeInputFiles_Slim.R', chdir = TRUE)


dir <- "~/Documents/My_Documents/UofToronto/SLiM/InputScripts/Sep19_Inputs"

setwd(dir)

rand.seeds <- c(
"1234567890",
"2345678901",
"3456789012",
"4567890123",
"5678901234",
"6789012345",
"7890123456",
"8901234567",
"9012345678",
"0123456789"
)

reps <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

pop.size <- 10000
mut.rate <- "3.5e-8"
recomb.rate <- "5e-8"
sample.dips <- "diploid"
N.number.generations <- 10
inds.to.sample <- 100


iterate <- 1
for(i in rand.seeds){

	# outcrossing:
	make.slim.input(filename.start="Sep19", 
		rand.seed=i, pop.size=pop.size, mut.rate=mut.rate, recomb.rate=recomb.rate, 
		mate.sys="outc", 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate]
	)
	
	# asexual 90%:
	make.slim.input(filename.start="Sep19", 
		rand.seed=i, pop.size=pop.size, mut.rate=mut.rate, recomb.rate=recomb.rate, 
		mate.sys="asex", 
		prop.mate.type="0.90", 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate]
	)

	# asexual 99%:
	make.slim.input(filename.start="Sep19", 
		rand.seed=i, pop.size=pop.size, mut.rate=mut.rate, recomb.rate=recomb.rate, 
		mate.sys="asex", 
		prop.mate.type="0.99", 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate]
	)

	# selfing 90%:
	make.slim.input(filename.start="Sep19", 
		rand.seed=i, pop.size=pop.size, mut.rate=mut.rate, recomb.rate=recomb.rate, 
		mate.sys="self", 
		prop.mate.type="0.90", 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate]
	)

	# selfing 99%:
	make.slim.input(filename.start="Sep19", 
		rand.seed=i, pop.size=pop.size, mut.rate=mut.rate, recomb.rate=recomb.rate, 
		mate.sys="self", 
		prop.mate.type="0.99", 
		total.N.gens=N.number.generations, samp.size=inds.to.sample, samp.type=sample.dips, rep=reps[iterate]
	)

	iterate <- iterate + 1
}





# see if these will work on grex:

RAM <- "18gb"
walltime <- "168:00:00"
server <- "grex"

multi.pbs(input.file.directory=dir, RAM=RAM, walltime=walltime, server=server)



# do the selfers on jasper:

RAM <- "5gb"
walltime <- "50:00:00"
server <- "jasper"

dir <- "~/Documents/My_Documents/UofToronto/SLiM/InputScripts/Sep19_Inputs/selfers"

multi.pbs(input.file.directory= dir, RAM=RAM, walltime=walltime, server=server)






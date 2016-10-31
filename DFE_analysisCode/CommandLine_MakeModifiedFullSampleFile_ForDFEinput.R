
args <- commandArgs(trailingOnly=TRUE)

# args 1 is basename of file
# args 2 is working directory (inputs dir)

setwd(as.character(args[2]))


taking.in.file <- paste(c("SampleOutput_", as.character(args[1])), collapse="")

spitting.out.file.name <- paste(c("ModifiedSampleOutput_", as.character(args[1])), collapse="")



# remove everything up to just before the last sample point "Mutations" line
# should start on line: "#OUT: 1000000 GS 2000"

start.line <- system(paste(c("grep -n '#OUT: 1000000 GS 2000' ", taking.in.file), collapse=""), intern=TRUE)
start.line <- as.numeric(unlist(strsplit(start.line, split=":"))[1])


# make the file:
system(paste(c("tail -n+", start.line, " ", taking.in.file, " > ", spitting.out.file.name), collapse=""))




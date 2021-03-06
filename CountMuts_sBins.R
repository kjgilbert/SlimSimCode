
count.muts.windows <- function(poly.mut.dat, fixed.mut.dat, generation, min.s=-1, max.s=0){	# default window is all the negative s values
	## ALSO GET TOTAL NUMBER OF DELET UTS IN EACH WINDOW
	# for fixed and poly muts
	if(is.null(fixed.mut.dat)){
		fixed.mut.count <- "none"	
	}else{
		fixed.mut.count <- "go ahead"	
		#	(all fixed muts are always present at the last full generation time point sampled)
		fixed.mut.dat <- fixed.mut.dat[fixed.mut.dat$gen.fixed <= as.numeric(generation) ,]
	}
		
	window.counts.poly <- NULL
	window.counts.fixed <- NULL
	for(i in 1:length(min.s)){
		poly.mut.count <- length(poly.mut.dat$seln_coeff[poly.mut.dat$seln_coeff > min.s[i] & poly.mut.dat$seln_coeff < max.s[i]])
		if(fixed.mut.count == "none"){
			fixed.mut.count <- c(0,0,0)
		}else{
			fixed.mut.count <- length(fixed.mut.dat$seln_coeff[fixed.mut.dat$seln_coeff > min.s[i] & fixed.mut.dat$seln_coeff < max.s[i]])
		}
		
		window.counts.poly <- c(window.counts.poly, poly.mut.count)
		window.counts.fixed <- c(window.counts.fixed, fixed.mut.count)
	}
	
	window.counts <- c(window.counts.poly, window.counts.fixed, (window.counts.poly + window.counts.fixed))
	names(window.counts) <- c(paste("poly.s.window.count", min.s, max.s, sep="_"), paste("fixed.s.window.count", min.s, max.s, sep="_"), paste("total.s.window.count", min.s, max.s, sep="_"))

	return(window.counts)
}

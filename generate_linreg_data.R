#Ran once to generate 100 ind. sub samples from pop. of 10k with maf of 0.05
	#Will use the same MACS output
library(SNPRelate)
library(GWASTools)

#'@param n.pop Total population size
#'@param n.sample Sample size
#'@param maf.min Minimum minor allele frequency
#'@param file.start Character string for file prefix
#'@param n.snp How many SNPs to keep?
#'@param r 
#'@param length Length of segment in base-pairs
#'@param seed Optional seed
#'@param save.start 
make.one.site <- function(n.pop, n.sample, maf.min, 
                          file.start, n.snp, r, macs.loc, tempdir="./",
                          length=50000, seed=NULL){
  all_genos <- rcc:::makemacsdata(N=n.pop, seed=seed, t=0.001, length=length, r=r, macs.loc=macs.loc)
	S <- 1:n.sample
	s_genos <- all_genos[S,]
	full.file <- paste(file.start, "_fullgen.gds", sep="")
	
	full.file <- paste(file.start, "_fullgen.gds", sep="")
	snpgdsCreateGeno(full.file, t(all_genos), sample.id=1:dim(all_genos)[1], snp.id=1:dim(all_genos)[2])

	gdsobj <- snpgdsOpen(full.file)
	snpset <- snpgdsLDpruning(gdsobj, maf=maf.min, missing.rate=0.05, 
	                          method="corr", slide.max.bp=10e6, ld.threshold=0.8, sample.id=S)
	closefn.gds(gdsobj)

	if(length(snpset$chr1) < n.snp) stop("Try increasing the length. Not enough SNPs generated.\n")

	gds <- GdsGenotypeReader(full.file)
	all_genos <- t(getGenotype(gds))
	cat(dim(all_genos), "\n")
	close(gds)

	s_genos <- all_genos[S, snpset$chr1[1:n.snp]]
	all_genos <- all_genos[, snpset$chr1[1:n.snp]]

	all.file <- paste(file.start, "_all.RData", sep="")
	sample.file <- paste(file.start, "_sample.RData", sep="")
	
	save(all_genos, file=all.file)
	save(s_genos, file=sample.file)
}

make.genos <- function(macs.loc="~/Dropbox/rcc-temp/for_jcgs/rcc/macs/", tempdir="./"){
  set.seed(953625)
  r_vals <- c(0.001, 0.01, 0.002, 0.004, 0.001)
  lengths <- rep(150000, 5)

  for(i in 1:5){
	  file.start <- paste("subset_", i, "_r", r_vals[i], "_n100", sep="")
	  make.one.site(n.pop=10000, n.sample=100, maf.min=0.05, tempdir=tempdir,
	              file.start=file.start, n.snp=100, r=r_vals[i], seed=seed, 
	              length=lengths[i], macs.loc=macs.loc)
  }

	i <- 1
	file.start <- paste("subset_", i, "_r", r_vals[i], "_n100", sep="")
	all.file <- paste(file.start, "_all.RData", sep="")
	sample.file <- paste(file.start, "_sample.RData", sep="")
	load(all.file)
	all_total_genos <- all_genos
	load(sample.file)
	sample_total_genos <- s_genos

  for(i in 2:5){
    cat(i, "\n")
	  file.start <- paste("subset_", i, "_r", r_vals[i], "_n100", sep="")
	  all.file <- paste(file.start, "_all.RData", sep="")
	  sample.file <- paste(file.start, "_sample.RData", sep="")
	  load(all.file)
	  all_total_genos <- cbind(all_total_genos, all_genos)
	  load(sample.file)
	  sample_total_genos <- cbind(sample_total_genos, s_genos)
  }
  save(all_total_genos, file="population_genos.RData")
  save(sample_total_genos, file="sample_genos.RData")
  
  
  for(i in 1:5){
    file.start <- paste("subset_", i, "_r", r_vals[i], "_n100", sep="")
    unlink(paste0(file.start, "_fullgen.gds"))
    unlink(paste0(file.start, "_all.RData"))
    unlink(paste0(file.start, "_sample.RData"))
  }
}
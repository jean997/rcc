#N: Number of individuals to generate
#length: length of segment in base pairs
#t: mutation rate per site per 4N generations
#From Xue, Y. et al. Curr. Biol. (2009) doi:10.1016/j.cub.2009.07.032. Human mutation rate might be about 1/3e7. 
#With N_e of 1e5 a reasonable value of t should be around 0.001

#r: recombination rate per site per 4N generations where N is effective population size.
#In humans r varies across the genome 

#From
#Jensen-Seaman, M. I., Furey, T. S., Payseur, B. A., Lu, Y., Roskin, K. M., Chen, C.-F., ... Jacob, H. J. (2004). Comparative recombination rates in the rat, mouse, and human genomes. Genome Research, 14(4), 528-38. doi:10.1101/gr.1970304
#Averages over whole chromosomes vary from 1 to 1.6 cM/Mb with higher rates in shorter chromosomes. 
#If we assume an eff. pop. size of 1e5, this would suggest r=0.4 to 0.64. 
#However, due to recombination hot spots, the rate will be lower in most segments. 
#The same paper estimates average rate of 0.5 per 5 Mb and 0.4 per 10 Mb segment which suggests r = 0.004 to 0.0016. 

#From 
#McVean, G. A. T., Myers, S. R., Hunt, S., Deloukas, P., Bentley, D. R., & Donnelly, P. (2004). The fine-scale structure of recombination rate variation in the human genome. Science (New York, N.Y.), 304(5670), 581-4. doi:10.1126/science.1092500
#Half of all reocmbinations occur in less that 10% of the genome. In hot spots it can be as much as 4 cM/Mb and in other places can be almost zero. 
#Allowing also for variagion in N_e anything from 1e-4 to 0.15 could potentially be reasonable
#seed: seed optional


#'@title Call MACS (Chen 2009) to generate genotype data
#'@param N
#'@param length
#'@param t
#'@param r
#'@param seed
#'@param macs.loc
#'@return A matrix of genotypes
#'@export
makemacsdata<-function(N,length=15000,t=0.001, r=0.001, seed=NULL, macs.loc="~/macs/", tempdir="./"){
	f=tempfile(tmpdir=tempdir)
	if(is.null(seed)){
		system(paste(macs.loc, "macs ",2*N," ",length, " -t ", t, " -r ",  r, " 2> /dev/null | ", macs.loc, "msformatter > ", f, sep=""))
	}else{
		stopifnot(is.numeric(seed))
		system(paste(macs.loc, "macs ",2*N," ",length, " -t ", t, " -r ",  r, " -s ", seed, " 2> /dev/null | ", macs.loc, "msformatter > ", f, sep=""))
	}
	input<-readLines(f)[-(1:6)]
	unlink(f)
	haplo<-do.call(rbind,lapply(strsplit(input,""),as.integer))
	diplo<-haplo[1:N,]+haplo[(N+1):(2*N),]
	return(diplo)
}




GdsGenotypeReader <- 
function (filename, genotypeDim, genotypeVar, snpIDvar, scanIDvar, 
          ...) 
{
  if (missing(filename)) 
    stop("filename is required")
  if (missing(genotypeVar)) 
    genotypeVar <- "genotype"
  if (missing(snpIDvar)) 
    snpIDvar <- "snp.id"
  if (missing(scanIDvar)) 
    scanIDvar <- "sample.id"
  input.gds <- is(filename, "gds.class")
  tmpobj <- GdsReader(filename)
  snpDim <- getDimension(tmpobj, snpIDvar)
  scanDim <- getDimension(tmpobj, scanIDvar)
  genoDim <- getDimension(tmpobj, genotypeVar)
  if (missing(genotypeDim)) {
    if (snpDim == scanDim) {
      genotypeDim <- ""
    }
    else if (all(genoDim == c(snpDim, scanDim))) {
      genotypeDim <- "snp,scan"
    }
    else if (all(genoDim == c(scanDim, snpDim))) {
      genotypeDim <- "scan,snp"
    }
    else {
      genotypeDim <- ""
    }
  }
  tryCatch(new("GdsGenotypeReader", tmpobj, genotypeDim = genotypeDim, 
               genotypeVar = genotypeVar, snpIDvar = snpIDvar, scanIDvar = scanIDvar, 
               ...), error = function(e) {
                 if (!input.gds) 
                   close(tmpobj)
                 stop(e)
               })
}



# Functions for the core statistical analysis

# ====================================================================
# Exported Functions

# --------------------------------------------------------------------
#' Detect differentially methylated regions (DMRs) from windowed read counts from MBD-isolated genome sequencing (MiGS/MBD-seq) and similar techniques
#'
#' After the counts have been pre-processed, this function performs all the analysis. Detailed results from intermediate steps are stored in the output list object to analyze method performance and provide input for the summary, export, and plotting functions.
#' @param samp Sample data.frame from readSampleInfo()
#' @param counts Preprocessed count data from getCounts()
#' @param reads Preprocessed reads/fragments data from getReads()
#' @param poifdr False discovery rate to use during initial filtering and frequency calling. Changing this value will change the threshold used for calling the presence of methylation.
#' @param stageone.p P-value cutoff for the pairwise testing in stage one.
#' @param anodev.p P-value cutoff for the analysis of deviance (ANODEV) in stage two testing.
#' @param post.p P-value cutoff for post-tests in stage two testing.
#' @param freq Fraction of samples within groups that must agree with respect to methylation status in order for "frequent" to be "TRUE" in the output DMR list.
#' @param minsize Minimum size of DMRs to report (in bp)
#' @param joindist Extend significant windows into DMRs over non-significant stage one windows between them up to this distance large (in bp)
#' @param adjust.var Name of a column present in "samp" that will be used as a covariate adjustment in the stage two ANODEV generalized linear model (GLM)
#' @param nperms Optional, perform this number of permutations after calling DMRs. Will create a data.table called "FDR" in the output list. See also maPerm(), maPermMerge(), and maPermFdr() for manual permutation running and FDR calculation.
#' @param perm.boot If nperms > 0 and if TRUE, perform bootstrapping (sampling with replacement). Otherwise, perform permutations (sampling without permutations)
#' @param ncore Number of parallel processes to use
#' @return A list containing detailed results from each stage of the analysis.
#' @export
methylaction <- function(samp, counts, reads, poifdr=0.1, stageone.p=0.05, anodev.p=0.05, post.p=0.05, freq=2/3, minsize=150, joindist=200, adjust.var=NULL, nperms=0, perm.boot=F, ncore=1)
{
	cov <- NULL
	stagetwo.method <- "co"
	winsize <- width(counts)[1]

	ngroups <- length(unique(samp$group))
	message("Starting analysis for a ",ngroups," group comparision")
	samp <- data.table(samp)
	if(sum(samp[,length(sample),by=group]$V1 > 1)!=ngroups){stop("Each group must contain more than one replicate")}
	groups <- unique(samp$group)
	if(!is.null(adjust.var)){if(!(adjust.var %in% colnames(samp))){stop("No column with name equal to adjust.var found in samp")}}

	args <- list(samp=samp, stagetwo.method=stagetwo.method, winsize=winsize, poifdr=poifdr, stageone.p=stageone.p, freq=freq, joindist=joindist, anodev.p=anodev.p, adjust.var=adjust.var, post.p=post.p, minsize=minsize, nperms=nperms, perm.boot=perm.boot, ncore=ncore, start=Sys.time())

	# Do Initial Filtering
	fdr.filter <- methylaction:::filter(counts, samp, poifdr)

	# Get Signal Bins Only
	message("Filtering to Signal Windows")
	filter.pass <- rowSums(t(t(as.matrix(values(counts)))>fdr.filter$cuts))>0

	# Compute Size Factors based on Filtered Regions - will always use these for all subsequent tests
	message("Computing size factors")
	sizefactors <- estimateSizeFactorsForMatrix(as.matrix(values(counts[filter.pass])))

	# Make list dividing up the data
	# Normalize counts first
	message("Saving Normalized Counts")
	cds <- newCountDataSet(as.matrix(values(counts[filter.pass])),samp$group)
	sizeFactors(cds) <- sizefactors
	normcounts <- counts(cds,normalized=TRUE)

	zero <- rowSums(as.matrix(values(counts)))==0
	wingr <- GRanges(seqnames(counts),IRanges(start(counts),end(counts)))
	windows <- list()
	windows$zero <- wingr[zero]
	windows$filtered <- wingr[!filter.pass & !zero]
	windows$signal <- counts[filter.pass]
	windows$signal.norm <- wingr[filter.pass]
	values(windows$signal.norm) <- normcounts

	rm(filter.pass, cds, normcounts, zero, wingr, counts)
	gc()

	# Testing Stage 1
	test.one <- methylaction:::testOne(samp=samp,bins=windows$signal,signal.norm=windows$signal.norm,chrs=unique(as.vector(seqnames(windows$signal))),sizefactors=sizefactors,stageone.p=stageone.p,joindist=joindist,minsize=minsize,ncore=ncore)

	# Testing Stage 2 + Methylation Modeling
	test.two <- methylaction:::testTwo(samp=samp, cov=cov, reads=reads, stagetwo.method=stagetwo.method, regions=test.one$regions, sizefactors=sizefactors, fragsize=fragsize, winsize=winsize, anodev.p=anodev.p, freq=freq, adjust.var=adjust.var, post.p=post.p, fdr.filter=fdr.filter, ncore=ncore)

	fdr <- NULL
	maperm <- NULL

	ma <- list(dmr=test.two$dmrcalled, fdr=fdr, args=args, data=list(windows=windows, fdr.filter=fdr.filter, sizefactors=sizefactors, test.one=test.one, test.two=test.two, maperm=maperm))

	# Remove some things from the output to reduce the size
	ma$data$test.two$dmrcalled <- NULL
	values(ma$data$windows$filtered) <- NULL
	#ma$data$windows$signal <- NULL

	if(nperms>0)
	{
		# call to maPerm to do the actual permutations
		ma$data$maperm <- methylaction:::maPerm(ma=ma,reads=reads,nperms=nperms,perm.boot=perm.boot,save=FALSE,ncore=ncore)
		ma$fdr <- methylaction:::maPermFdr(ma=ma,maperm=ma$data$maperm,recut.p=0.05)
	}

	# Output results
	#ma <- list(fdr.filter=fdr.filter, sizefactors=sizefactors, windows=windows, test.one=test.one, test.two=test.two)
	ma$args$end <- Sys.time()
	ma$args$hours <- as.numeric(difftime(ma$args$end,ma$args$start,units="hours"))

	message("Output list size in memory=",methylaction:::sizein(ma))
	return(ma)
}
# --------------------------------------------------------------------
# ====================================================================

# ====================================================================
# Internal Functions

# --------------------------------------------------------------------
# Initial Filtering
filter <- function(counts, samp, poifdr)
{
	# Initial Filtering
	message("Computing window count frequencies")
	countints <- 0:10
	histo <- do.call(rbind,lapply(countints, FUN=function(x) {message(x);colSums(as.matrix(values(counts))==x);}))
	histo <- rbind(histo,colSums(as.matrix(values(counts))>max(countints)))
	rownames(histo) <- c(countints,paste0(">",max(countints)))
	counts.tab <- lapply(1:ncol(histo),function(x) data.frame(count=countints[-1],freq=histo[2:length(countints)-1,x]))
	names(counts.tab) <- colnames(histo)

	# This filtering equation based on code by Yaomin Xu
	dpt.poi <- function(histo, counts.tab,cutoff=poifdr)
	{
		est.pois.lambda <- function(tab)
		{
			f1 <- with(tab, freq[count==1])
			f2 <- with(tab, freq[count==2])
			2*f2/f1
		}
		est.fdr <- function(n,tb)
		{
			lambda <- est.pois.lambda(tb)
			n.total <- with(tb, freq[count==1])/dpois(1, lambda)
			n.total*dpois(n,lambda)/with(tb, freq[count == n])
		}
		fdr <- sapply(counts.tab, function(x) sapply(counts.tab[[1]]$count, est.fdr, tb=x))
		ifsatisfy <- fdr < cutoff
		cuts <- apply(ifsatisfy, MARGIN=2, FUN=function(x) max(seq(1,length(x))[!x]))
		list(histo=histo,fdr=round(fdr,5),cuts=cuts)
	}
	fdr.filter <- dpt.poi(histo, counts.tab, cutoff=poifdr)
	return(fdr.filter)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Stage one testing
testOne <- function(samp,bins,signal.norm,chrs,sizefactors,stageone.p=0.05,minsize=150,joindist=200,ncore=3)
{
	message("Begin stage one testing")

	# Perform pairwise tests
	bins.mat <- as.matrix(values(bins))

	# Compute all the pairwise comps
	comps <- methylaction:::getGroupComps(unique(samp$group))
	dotest <- function(x)
	{
		a <- x[1]
		b <- x[2]
		message("Testing ",a," vs ",b)
		out <- methylaction:::testDESeq(bins.mat,samp$group,a=a,b=b,prefix=paste0("deseq-",a," vs ",b),sizefactors=sizefactors,ncore=ncore)
		return(out)
	}
	testres <- apply(comps$combos,1,dotest)
	names(testres) <- comps$strings

	if(!all(sapply(testres,nrow)==length(bins))){stop("Ran out of memory during testing, try reducing ncore")}

	message("Calling Patterns")
	patt <- methylaction:::callPatternsN(res=testres,cutoff=stageone.p)

	# Make means columns
	colgroups <- lapply(unique(samp$group),function(x) samp[group==x,]$sample)
	names(colgroups) <- unique(samp$group)
	counts.means <- do.call(cbind, lapply(colgroups, function(i) rowMeans(as.matrix(values(signal.norm[,i])))))
	colnames(counts.means) <- paste0(levels(samp$group),".mean")

	# Combine so we get these in the output
	patt <- data.table(counts.means,patt)

	# join adjacent equivalent patterns
	message("Reduction by pattern and disjoing regions")
	bins.gr <- bins
	values(bins.gr) <- NULL
	bins.gr$pattTestOne <- patt$patt
	
	# split out by pattern
	allone <- paste(as.character(rep(1,length(unique(samp$group)))),collapse="")
	sigpatt <- bins.gr[!(bins.gr$pattTestOne %in% c("ambig","000or111",allone))]
	pattrows <- sigpatt$pattTestOne
	values(sigpatt) <- NULL
	bins.bypatt <- split(sigpatt,pattrows)

	# reduce each pattern with itself within the gap distance
	bins.red <- lapply(bins.bypatt,reduce,min.gapwidth=joindist)
	for(i in 1:length(bins.red)){bins.red[[i]]$pattTestOne <- names(bins.red)[i]}

	# only keep reduced ranges if minsize or larger
	bins.red <- lapply(bins.red,function(x) x[width(x)>=minsize])

	# deal with conflicting extensions using disjoin()
	bins.red.gr <- do.call(c,unname(bins.red))
	regions.gr <- disjoin(bins.red.gr)
	regions <- regions.gr[width(regions.gr)>=minsize]

	# output: genomic ranges to pass to stage 2

	patterns <- bins
	values(patterns) <- patt

	test.one <- list(patterns=patterns, regions=regions)

	return(test.one)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# stage two testing
testTwo <- function(samp,cov,reads,stagetwo.method="co",regions,sizefactors,fragsize,winsize,anodev.p,adjust.var,post.p, fdr.filter, freq,ncore)
{
	message("Begin stage two testing")

	message("Recounting in regions with method=",stagetwo.method)

	if(stagetwo.method=="co")
	{
		recounts <- as.matrix(values(getCounts(samp=samp,reads=reads,ranges=regions,ncore=ncore)))
	} else
	{
		stop("Bad stagetwo.method type")
	}

	# Do an ANOVA-like framework
	# First test for ANY difference with an ANOVA-style DESeq test (ANODEV since these are GLM based negative binomial tests)
	# Adjust this for multiple testing
	# Move on the significant windows to post testing, adjust within each test only then
	# Call pattern based on these values
	testDESeqANODEV<- function(recounts,groups,covar, ncore)
	{
		# chunk up data to allow parallel testing for large data sets
		chunkids <- (seq(nrow(recounts))-1) %/% 500

		dotest <- function(rows)
		{
			singleton <- FALSE
			if(nrow(recounts[rows,,drop=F])==1)
			{
				# Workaround for estimateDispersions - will not work and can not test if only one feature, so duplicate the feature we have as an estimation in the singleton case
				singleton <- TRUE
				thecounts <- rbind(recounts[rows,,drop=F],recounts[rows,,drop=F])
				cds <- newCountDataSet(thecounts,groups)
			} else
			{
				# make master object container
				cds <- newCountDataSet(recounts[rows,,drop=F],groups)
			}

			# make master object container
			#cds <- newCountDataSet(recounts[rows,,drop=F],groups)

			# estimate the size factors for the library
			# do we want to use size factors from previously for the whole library?
			sizeFactors(cds) <- sizefactors

			# pull back DEseq normalized counts - use for visualization
			# probably want to make the WIG files out of these
			#head( counts( cds, normalized=TRUE ) )

			# estimate the dispersions for each gene
			message("Estimating dispersions")
			cds <- suppressWarnings(estimateDispersions(cds,fitType="local",sharingMode="gene-est-only"))

			# perform the testing
			# compare model w/ factor versus null model
			if(!is.null(covar))
			{
				message("Fitting full GLM w/ adjustments")
				resFull <- fitNbinomGLMs( cds, modelFormula=count~groups+covar)
				message("Fitting reduced GLM w/ adjustments")
				resReduced <- fitNbinomGLMs( cds, modelFormula=count~covar)
				message("Performing test")
				res <- nbinomGLMTest(resFull, resReduced)
				if(singleton==TRUE)
				{
					res <- res[1]
				}
				return(res)
		
			} else {
				message("Fitting full GLM")
				resFull <- fitNbinomGLMs( cds, modelFormula=count~groups)
				message("Fitting reduced GLM")
				resReduced <- fitNbinomGLMs( cds, modelFormula=count~1)
				message("Performing test")
				res <- nbinomGLMTest(resFull, resReduced)
				if(singleton==TRUE)
				{
					res <- res[1]
				}
				return(res)
			}
		}
		finalres <- mclapply(0:max(chunkids),function(x) {message("Testing chunk ",x," of ",max(chunkids));dotest(chunkids==x);},mc.cores=ncore)
		finalres <- do.call(c,finalres)
		return(finalres)
	}

	if(is.null(adjust.var))
	{
		covar <- NULL
	} else {
			message("Adding adjustment for variable \"", adjust.var,"\" into the GLM")
			covar <- factor(as.data.frame(samp)[,adjust.var,drop=T])
	}
	

	anodev <- testDESeqANODEV(recounts=recounts,groups=samp$group,covar=covar,ncore=ncore)
	anodev.padj <- p.adjust(anodev,method="fdr")
	#table(anodev<0.1)
	# Filter out if it does not pass the ANODEV

	#df <- data.frame(noadj=p.adjust(anodev1,method="fdr")<0.001,adj=p.adjust(anodev,method="fdr")<0.001)
	#table(rowSums(df))

	regions$anodev.p <- anodev
	regions$anodev.padj <- anodev.padj

	# Get norm counts
	cds <- newCountDataSet(recounts,samp$group)
	sizeFactors(cds) <- sizefactors
	#values(regions) <- cbind(as.data.frame(values(regions)),counts(cds,normalized=TRUE))
	normcounts <- counts(cds,normalized=TRUE)

	anodev.keep <- (anodev.padj<anodev.p) & !(is.na(anodev.padj))

	if(sum(anodev.keep)==0){message("Found no DMRs in stage two ANODEV")}
	if(sum(anodev.keep)==0){out<-"NoDmrs"; return(out);}

	test.two <- list()
	test.two$ns <- regions[!anodev.keep]
	test.two$ns.counts <- regions[!anodev.keep]
	values(test.two$ns.counts) <- normcounts[!anodev.keep,]

	regions.sig <- regions[anodev.keep]
	recounts.sig <- recounts[anodev.keep,,drop=F]

	# a vs b
	comps <- methylaction:::getGroupComps(unique(samp$group))
	dotest <- function(x)
	{
		a <- x[1]
		b <- x[2]
		message("Testing ",a," vs ",b)
		#s <- strsplit(x,"")[[1]]
		out <- methylaction:::testDESeq(counts=recounts.sig,groups=samp$group,a=a,b=b,prefix=paste0("deseq-",a," vs ",b),sizefactors=sizefactors,ncore=ncore)
		return(out)
	}
	testres <- apply(comps$combos,1,dotest)
	names(testres) <- comps$strings

	if(!all(sapply(testres,nrow)==nrow(recounts.sig))){stop("Ran out of memory during testing, try reducing ncore")}

	patt <- methylaction:::callPatternsN(res=testres, cutoff=post.p)

	test.two$sig <- regions[anodev.keep]
	values(test.two$sig) <- patt
	test.two$sig.counts <- regions[anodev.keep]
	values(test.two$sig.counts) <- normcounts[anodev.keep,,drop=F]

	#as.data.frame(table(patt$patt))

	dmr <- data.frame(chr=seqnames(regions.sig),start=start(regions.sig),end=end(regions.sig),width=width(regions.sig),anodev.padj=regions.sig$anodev.padj,pattern=patt$patt)
	
	# Output list of sig DMRs

	# Make frequency columns
	colgroups <- lapply(unique(samp$group),function(x) samp[group==x,]$sample)
	names(colgroups) <- unique(samp$group)

	# Freq Calling via POI
	#cuts <- as.matrix(width(dmr)/50) %*% cuts
	cnt <- recounts.sig/(dmr$width/winsize)

	#meth <- t(t(counts)>cuts)
	# call yes/no for each sample based on poi cuts
	meth <- t(t(cnt)>fdr.filter$cuts)

	# call freqs
	meth.freq <- do.call(cbind, lapply(colgroups, function(i) rowSums(meth[,i,drop=F])))

	# call percentages
	meth.per <- t(t(meth.freq)/as.vector(table(samp$group)))
	colnames(meth.per) <- paste0(colnames(meth.per),".per")

	# Make means columns - for the per window counts
	perwin.means <- do.call(cbind, lapply(colgroups, function(i) rowMeans(cnt[,i,drop=F])))
	colnames(perwin.means) <- paste0(levels(samp$group),".perwin.mean")

	dmrfreq <- cbind(dmr,cnt,perwin.means,meth.freq,meth.per)
	#dmrfreq$dmrid <- 1:nrow(dmr)

	# Want to know sharpness - should be all 0 groups are no more than 1/3 meth, and all 1 groups are no less than 2/3 meth
	patts <- as.character(unique(dmrfreq$pattern))
	allone <- paste(as.character(rep(1,length(unique(samp$group)))),collapse="")
	patts <- patts[!(patts %in% c("ambig","000or111",allone))]

	sharpnessN <- function(patt,meth=freq,unmeth=(1-freq))
	{
		myfreq <- dmrfreq[dmrfreq$patt==patt,paste0(unique(samp$group),".per")]

		chars <- strsplit(patt,"")[[1]]
		checkit <- chars
		checkit[chars=="0"] <- unmeth
		checkit[chars=="1"] <- meth

		lt <- t(t(myfreq)<=checkit)
		gt <- t(t(myfreq)>=checkit)

		# Call each group as sharp or not, they must all be sharp to call the site as sharp
		sharp <- rowSums(cbind(lt[,chars=="0",drop=F],gt[,chars=="1",drop=F]))==length(unique(samp$group))
		
		ret <- dmrfreq[dmrfreq$pattern==patt,,drop=F]
		ret$frequent <- sharp
		ret
	}
	sharpness <- function(patt,meth=freq,unmeth=(1-freq))
	{
		myfreq <- dmrfreq[dmrfreq$pattern==patt,c("benign.per","low.per","high.per")]
		chars <- strsplit(patt,"")[[1]]
		checkit <- chars
		checkit[chars=="0"] <- unmeth
		checkit[chars=="1"] <- meth

		lt <- t(t(myfreq)<=checkit)
		gt <- t(t(myfreq)>=checkit)

		# Call each group as sharp or not, they must all be sharp to call the site as sharp
		sharp <- rowSums(cbind(lt[,chars=="0",drop=F],gt[,chars=="1",drop=F]))==3
		
		ret <- dmrfreq[dmrfreq$pattern==patt,]
		ret$frequent <- sharp
		ret
	}

	res <- lapply(patts,sharpnessN)
	dmrcalled <- do.call(rbind,unname(res))
	dmrcalled$dmrid <- 1:nrow(dmrcalled)

	# Make means columns
	counts.means <- do.call(cbind, lapply(colgroups, function(i) rowMeans(recounts.sig[,i,drop=F])))
	colnames(counts.means) <- paste0(unique(samp$group),".mean")

	# Final output!
	dmr <- cbind(dmr,counts.means)

	test.two$dmr <- makeGRanges(dmr)
	test.two$dmrcalled <- makeGRanges(dmrcalled)

	return(test.two)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# DEseq for 1 pairwise test
testDESeq <- function(counts,groups,a,b,prefix,sizefactors,ncore)
{
	# chunk up data to allow parallel testing for large data sets
	chunkids <- (seq(nrow(counts))-1) %/% 50000

	# need to use the same size factors for everything

	dotest <- function(rows)
	{
		singleton <- FALSE
		if(nrow(counts)==1)
		{
			# Workaround for estimateDispersions - will not work and can not test if only one feature, so duplicate the feature we have as an estimation in the singleton case
			singleton <- TRUE
			counts <- rbind(counts,counts)
			cds <- newCountDataSet(counts,groups)
		} else
		{
			# make master object container
			cds <- newCountDataSet(counts[rows,,drop=F],groups)
		}

		# estimate the size factors for the library
		#cds <- estimateSizeFactors(cds)
		sizeFactors(cds) <- sizefactors

		# estimate the dispersions for each gene
		#message("Estimating dispersions")
		cds <- suppressWarnings(estimateDispersions(cds,fitType="local",sharingMode="gene-est-only"))
		#cds <- suppressWarnings(estimateDispersions(cds,fitType="local",sharingMode="gene-est-only"))
		#cds<-estimateDispersions(cds,method="pooled-CR",sharingMode="gene-est-only",fitType="local")

		# perform the testing
		#message("Performing test")
		res <- nbinomTest(cds,a,b)

		if(singleton==TRUE)
		{
			res <- res[1,,drop=F]
		}

		return(res)
	}
	finalres <- mclapply(0:max(chunkids),function(x) {message("Testing chunk ",x," of ",max(chunkids));dotest(chunkids==x);},mc.cores=ncore)
	finalres <- rbindlist(finalres)
	finalres[,id:=NULL]
	finalres[,padj:=NULL]

	# Remove NA p-value by setting to 0 (when b mean is 0) or 1, whichever is appropriate

	return(finalres)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Function to generate all pairwise comparisions for a set of groups
# Where groups is a vector of group IDs (unique)
getGroupComps <- function(groups)
{
	ngroups <- length(groups)
	combos <- t(combn(1:ngroups,m=2))
	combos <- t(apply(combos,1,function(x) groups[x]))
	combostrings <- apply(combos,1,function(x) paste0(x[1],"_vs_",x[2]))
	return(list(combos=combos,strings=combostrings))
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Function to generate all possible pattern codes given nGroups
getGroupPatterns <- function(ngroups)
{
	ones <- 0:ngroups

	# Function to return all combinations for a given number of "1" digits in the pattern
	# We will apply over all possible numbers of "1" digits to get all the pattern codes
	digits <- function(x)
	{
		if(x==0)
		{
			return(rep(0,ngroups))
		} else if (x==ngroups)
		{
			return(rep(1,ngroups))
		} else
		{
			pool <- c(rep(0,ngroups-x),rep(1,x))
			perms <- gtools::permutations(n=ngroups,r=ngroups,v=pool,set=FALSE)
			perms <- unique(perms)
			return(perms)
		}
	}
	codes <- lapply(ones,digits)
	codes <- do.call(rbind,codes)
	return(codes)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# Decision table
getDecisionTable <- function(mygroups,human=FALSE)
{
	codes <- getGroupPatterns(length(mygroups))

	getRule <- function(x)
	{
		#message(x)
		#x <- codes[2,]
		combos <- t(combn(1:length(mygroups),m=2))
		namecombos <- t(apply(combos,1,function(y) mygroups[y]))
		combos <- t(apply(combos,1,function(y) x[y]))
		namecombos <- apply(namecombos,1,function(x) paste0(x[1],"_vs_",x[2]))

		test <- function(z)
		{
			if((z[1]==0)&(z[2]==0))
			{
				if(human==TRUE)
				{
					return("NS")
				} else
				{
					return(0)
				}
			} else if((z[1]==1)&(z[2]==0))
			{
				if(human==TRUE)
				{
					return("SIG decrease")
				} else
				{
					return(-1)
				}
			} else if((z[1]==0)&(z[2]==1))
			{
				if(human==TRUE)
				{
					return("SIG increase")
				} else
				{
					return(1)
				}
			} else if((z[1]==1)&(z[2]==1))
			{
				if(human==TRUE)
				{
					return("NS")
				} else
				{
					return(0)
				}
			}
		}
		rule <- apply(combos,1,test)
		names(rule) <- namecombos
		rule <- matrix(rule,nrow=1)
		colnames(rule) <- namecombos
		rule <- data.frame(pattern=paste(x,collapse=""),rule)
		return(rule)
	}
	codes.split <- split(codes,row(codes))
	ja <- lapply(codes.split,getRule)
	ja <- do.call(rbind,ja)
	return(ja)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
# call patterns genome wide, generalized to n-groups
callPatternsN <- function(res,cutoff=0.05)
{
	# add directions to testres
	getDirection <- function(x)
	{
		x$direction <- 0
		x[baseMeanA<baseMeanB,]$direction <- 1
		x[x$baseMeanA>x$baseMeanB,]$direction <- -1
		return(x)
	}
	res <- lapply(res,getDirection)

	# combine columns accross comparisons
	getCols <- function(x)
	{
		out <- data.table(p=x$pval,l2fc=x$log2FoldChange,sig=x$pval<cutoff,dir=x$direction)
		out[is.na(p),sig:=FALSE]
		return(out)
	}
	tests <- lapply(res,getCols)

	# filter using decision codes:
	# -1 = sig down
	# 0 = NS
	# 1 = sig up

	getCodes <- function(x)
	{
		x[,code:=dir]
		x[sig==FALSE,code:=0]
		return(x)
	}
	tests <- lapply(tests,getCodes)

	# Call pattern codes based on the decision codes
	dc <- methylaction:::getDecisionTable(unique(samp$group))

	str <- lapply(tests,function(x) x$code)
	str <- do.call(cbind,str)

	tests <- lapply(names(tests),function(x) setnames(tests[[x]],paste0(x,"_",colnames(tests[[x]]))))
	mypatt <- do.call(cbind,tests)
	mypatt[,patt:="ambig"]

	for(i in 1:nrow(dc))
	{
		message("Deciding pattern: ",dc[i,]$pattern)
		x <- as.vector(t(dc[i,-1]))
		dec <- t(t(str)==x)
		yespatt <- rowSums(dec)==ncol(dec)
		mypatt[yespatt,patt:=dc[i,]$pattern]	
	}

	mypatt[patt=="000",]$patt <- "000or111"
	mypatt$patt <- as.character(mypatt$patt)
	mypatt
}
# --------------------------------------------------------------------

# ====================================================================

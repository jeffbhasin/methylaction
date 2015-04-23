# Functions for the core statistical analysis

# --------------------------------------------------------------------
#' Detect differentially methylated regions (DMRs) from windowed read counts from MBD-isolated genome sequencing (MiGS/MBD-seq)
#'
#' Once the counts have been pre-processed, this function performs all the analysis. Detailed results from intermediate steps are stored in the output list object to analyze method performance and provide input for the summary and plotting functions.
#' @param samp Description of samples from readSampleInfo()
#' @param counts Preprocessed count data from getCounts()
#' @param bsgenome b-string genome (bsgenome) object for the genome
#' @param fragsize The average fragment length selected for in the sequencing experiment (used to extend reads when re-counting for regions in stage two testing).
#' @param winsize Size of the windows used when counting.
#' @param poifdr False discovery rate to use during initital filtering.
#' @param stageone.p P-value cutoff for stage one testing.
#' @param anodev.p P-value cutoff for the analysis of deviance (ANODEV) in stage two testing (ignored for two group comparisions).
#' @param post.p P-value cutoff for post-tests (or for the single test stage two test in the two group case).
#' @param minsize Minimum size for a reported region.
#' @param ncore Number of cores to use.
#' @return A list containing detailed results from each stage of the analysis.
#' @export
methylaction <- function(samp, counts, reads=NULL, winsize, poifdr=0.1, stageone.p=0.05, joindist=200, anodev.p=0.05, post.p=0.05, freq=2/3, adjust.var=NULL, minsize=150, nperms=0, perm.boot=F, ncore=1)
{
	cov=NULL
	stagetwo.method="co"

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

	# Testing Stage 2 + Methylation Modelling
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

# --------------------------------------------------------------------
# Given an ma object, perform n permutations
maPerm <- function(ma,reads,nperms,combos=NULL,save=T,perm.combo=F,perm.boot=F,ncore)
{
	args <- ma$args
	test.two <- ma$data$test.two
	samp <- ma$args$samp
	adjust.var <- ma$args$adjust.var
	windows <- ma$data$windows
	sizefactors <- ma$data$sizefactors
	stagetwo.method <- ma$args$stagetwo.method
	anodev.p <- ma$args$anodev.p
	post.p <- ma$args$post.p
	fdr.filter <- ma$data$fdr.filter
	fragsize <- ma$args$fragsize
	winsize <- ma$args$winsize
	dmrcalled <- ma$dmr
	freq <- ma$args$freq

	if(!is.null(combos))
	{
		nperms <- nrow(combos)
		perm.boot <- F
		perm.combo <- F
		stopifnot(ncol(combos)==nrow(samp))
	}

	doperm <- function(perm)
	{
		message("Permutation number ",perm)
		#host <- system2("hostname",stdout=TRUE)
		#write(paste0("Started permutation ",perm," on ",host," at ",date()),file=paste0(host,".txt"),append=TRUE)
		if(perm.combo==TRUE)
		{
			message("Drawing permutation order from combination space")
			ss <- as.vector(table(samp$groupcode))
			rand <- getComboSpacePerms(an=ss[1],bn=ss[2],cn=ss[3],nperms=1,ncore=ncore)[[1]]
		} else if (!is.null(combos)) {
			rand <- combos[perm,]
			message("Using pre-set combination vector:")
			print(rand)
		} else {
			rand <- sample(1:nrow(samp),nrow(samp),replace=perm.boot)
		}
		
		mysamp <- samp
		mysamp$sample <- samp$sample[rand]
		mysamp$bam <- samp$bam[rand]
		if(!is.null(adjust.var))
		{
			mysamp[[adjust.var]] <- samp[[adjust.var]][rand]
		}

		mybins <- windows$signal
		values(mybins) <- values(mybins)[,rand]
		mysig <- windows$signal.norm
		values(mysig) <- values(mysig)[,rand]
		mysizes <- sizefactors[rand]

		# rename samples so there are no identical names in the bootstrap case
		mysamp$sample <- paste0("BS",1:nrow(samp),mysamp$sample)
		names(mysizes) <- mysamp$sample
		colnames(values(mybins)) <- mysamp$sample
		colnames(values(mysig)) <- mysamp$sample

		# need to check for and randomize cov/reads here too
		mycov <- NULL
		myreads <- NULL
			
		if(stagetwo.method=="co")
		{
			# shuffle reads
			myreads <- reads[rand]
			names(myreads) <- mysamp$sample
		} else if(stagetwo.method=="ac")
		{
			# shuffle cov
			mycov <- cov[rand]
		}
		#message("T1")
		test.one <- methylaction:::testOne(samp=mysamp,bins=mybins,signal.norm=mysig,chrs=unique(as.vector(seqnames(mybins))),sizefactors=mysizes,stageone.p=ma$args$stageone.p,joindist=ma$args$joindist,minsize=ma$args$minsize,ncore=ncore)
		#message("T2")
		test.two <- methylaction:::testTwo(samp=mysamp, cov=mycov, reads=myreads, stagetwo.method=stagetwo.method, regions=test.one$regions, freq=freq, sizefactors=mysizes, fragsize=fragsize, winsize=winsize, anodev.p=anodev.p, post.p=post.p, adjust.var=adjust.var, fdr.filter=fdr.filter,ncore=ncore)
		if(is.list(test.two))
		{
			ret <- test.two$dmrcalled
		} else
		{
			ret <- test.two
		}
			
		rm(test.one,test.two,mybins,mysig)
		gc()
		#maperm2[[perm]] <<- ret
		#save(rand,file="permdump.rd",compress=T)
		return(ret)
	}

	deathproof <- function(x)
	{
		ret <- tryCatch(
			doperm(x),
			error = function(e) 
			{
				message("Permutation ",x," had error: ",e$message)
				character()
			}
		)
		return(ret)
	}

	maperm <- lapply(1:nperms, deathproof)

	#ma$fdr <- fdr
	#ma$data$maperm <- maperm

	if(save==T)
	{
		# save Rd to disk that can be loaded using maPermMerge
		file <- paste0("maperm_",R.utils::getHostname.System(),"_",as.integer(as.POSIXct(Sys.time())),".rd")
		save(maperm,file=file,compress=T)
		message("Saved maperm object in ",file)
	} else
	{
		return(maperm)
	}
}
# --------------------------------------------------------------------

# Recompute FDRs based on permutation data in ma$data$maperm
maPermFdr <- function(ma,maperm,recut.p=0.05)
{
	samp <- ma$args$samp
	#maperm <- ma$data$maperm
	nperms <- length(maperm)
	dmrcalled <- ma$dmr

	checkdead <- sapply(maperm,class)
	dead <- (checkdead != "GRanges")&(maperm!="NoDmrs")
	message("Out of ",nperms," perms, ",sum(dead)," died and were ignored")
	maperm <- maperm[!dead]

	nodmr <- maperm=="NoDmrs"
	message("Out of ",nperms," perms, ",sum(nodmr)," found no DMRs of any pattern")

	ma$args$permerrors <- sum(dead)
	ma$args$permnodmrs <- sum(nodmr)

	# make observed event table
	gettab <- function(dmrcalled)
	{
		dmrcalled <- dmrcalled[dmrcalled$anodev.padj<recut.p]
		dmrcalled$frequent <- factor(dmrcalled$frequent,levels=c("FALSE","TRUE"))
		matab <- as.matrix(table(dmrcalled$pattern,dmrcalled$frequent))
		matab <- cbind(matab,rowSums(matab))
		colnames(matab) <- c("other","frequent","all")
		matab <- rbind(matab,colSums(matab))
		rownames(matab)[nrow(matab)] <- "all"
		return(matab)
	}

	# make expected event table - mean of all permutations
	realtab <- gettab(dmrcalled)
	allone <- paste(as.character(rep(1,length(unique(samp$group)))),collapse="")
	realtab <- realtab[!rownames(realtab) %in% c("000or111","ambig","1111",allone),]

	# If no DMR, put out a blank table
	#mapermout <- maperm #this one won't be subsetted
	maperm <- maperm[!nodmr]

	if(sum(nodmr)>0)
	{
		notabs <- lapply(1:sum(nodmr),function(x) {k<-realtab; apply(k,c(1,2),function(x) 0);})
	}

	# Make perm tables
	# make sure the perms report missing levels so table() gives same output for everything
	maperm <- lapply(maperm,function(x){
		x$pattern <- factor(x$pattern,levels=levels(dmrcalled$pattern));
		return(x);
	})

	permtab <- lapply(maperm,gettab)
	allone <- paste(as.character(rep(1,length(unique(samp$group)))),collapse="")
	permtab <- lapply(permtab,function(x) x[!(rownames(x) %in% c("000or111","ambig","1111",allone)),])

	# Need to make sure they are ordered right before doing the division
	permtab <- lapply(permtab,function(x) x[rownames(realtab),])

	# Add the blanks back in if they're there
	if(sum(nodmr)>0)
	{
		permtab <- c(permtab, notabs)
	}

	permmeans <- apply(simplify2array(permtab), c(1,2), mean)
	permmeans <- permmeans[rownames(realtab),]
	permsds <- apply(simplify2array(permtab), c(1,2), sd)
	permsds <- permsds[rownames(realtab),]
	permcv <- permsds/permmeans
	permcv <- permcv[rownames(realtab),]


	# calculate FDR by pattern and overall
	#expected <- Reduce("+",permtab2)/nperms
	fdrpercents <- (permmeans/realtab)*100

	# Make neat summary tables
	real.m <- reshape::melt.matrix(realtab)
	colnames(real.m) <- c("pattern","type","value")
	real.m$var <- "nDMRs"
	mean.m <- reshape::melt.matrix(permmeans)
	colnames(mean.m) <- c("pattern","type","value")
	mean.m$var <- "permMean"
	sd.m <- reshape::melt.matrix(permsds)
	colnames(sd.m) <- c("pattern","type","value")
	sd.m$var <- "permSD"
	cv.m <- reshape::melt.matrix(permcv)
	colnames(cv.m) <- c("pattern","type","value")
	cv.m$var <- "permCV"
	fdr.m <- reshape::melt.matrix(fdrpercents)
	colnames(fdr.m) <- c("pattern","type","value")
	fdr.m$var <- "FDRpercent"
	longdf <- rbind(real.m,mean.m,sd.m,cv.m,fdr.m)
	longdf$var <- factor(longdf$var,levels=unique(longdf$var))
	fdr <- reshape::cast(longdf,formula="pattern+type~var",value="value")
	#ma$perms$dmrs <- maperm
	#ma$perms$fdr <- fdr
	# Remove count columns from maperm to save space and not duplicate this data
	#maperm2 <- lapply(maperm,function(x){
		#values(x) <- as.data.frame(values(x))[,!(colnames(values(x)) %in% samp$sample)]
		#return(x);
	#})

	return(fdr)
}

# --------------------------------------------------------------------
# Draws permtuation orders out of combination space by running a series of chooses
# This may be preferable to simply shuffling the sample label list because it ensures each permtation ends up with distinct configurations of groups
# a: size of group 1
# b: size of group 2
# c: size of group 3
# nperms: how many permutations to provide
getComboSpacePerms <- function(an,bn,cn,nperms,ncore=1)
{
	nn <- an+bn+cn
	totalperms <- choose(nn,an)*choose(nn-an,bn)
	#message("There are ",formatC(totalperms,format="d",big.mark=",")," total possible combinations")
	samples <- 1:nn

	# Generate all possible ways to draw the first group
	co1 <- gtools::combinations(n=nn,r=an,v=samples)

	# Want to pick one at random, then pick one of the next possible choices at random up to the number of permtuations
	rand <- sample(1:nrow(co1),nperms,replace=F)

	dorow <- function(x)
	{
		#message(x)
		# Get one draw of the first round
		a <- co1[x,]

		# Get everything that's left after the first draw
		nota <- samples[!(samples %in% a)]

		# Get all combos for the second draw
		b <- gtools::combinations(n=22-7,r=6,v=nota)

		# Get the remainder for group 3
		c <- t(apply(b,1,function(x) samples[!(samples %in% a) & !(samples %in% x)]))

		# Build Output table
		tab <- cbind(t(replicate(nrow(b),a)),b,c)
		colnames(tab) <- c(rep("a",7),rep("b",6),rep("c",9))

		# Report one at random
		ret <- tab[sample(1:nrow(tab),1),]
		return(ret)
	}
	permcos <- mclapply(rand,dorow,mc.cores=ncore)
	#out <- do.call(rbind,permcos)
	#return(out)
	return(permcos)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
binnedAverage <- function(bins, numvar)
{
	stopifnot(is(bins, "GRanges"))
	stopifnot(is(numvar, "RleList"))
	stopifnot(identical(seqlevels(bins), names(numvar)))
	bins_per_chrom <- split(ranges(bins), seqnames(bins))
	means_list <- lapply(names(numvar),
	function(seqname)
	{
		views <- Views(numvar[[seqname]],bins_per_chrom[[seqname]])
		viewMeans(views)
	})
	new_mcol <- as.integer(round(unsplit(means_list, as.factor(seqnames(bins)))))
	#mcols(bins)[[mcolname]] <- new_mcol
	#bins
	return(new_mcol)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------

# Count for any GRanges object
regionCounts <- function(cov, regions)
{
	ba <- mclapply(cov, function(x) binnedAverage(bins=regions, numvar=x),mc.cores=5)
	ba2 <- mclapply(ba, Rle)
	ba <- do.call(cbind,ba)
	#ba <- round(ba,digits=0)
	return(ba)
}

# Count in windows spanning the genome
windowCounts <- function(reads, bsgenome, chrs, winsize)
{
	gb <- Repitools::genomeBlocks(genome=bsgenome, chrs=chrs, width=winsize)
	cov <- mclapply(reads,coverage,mc.cores=ncore)
	countmat <- regionCounts(cov=cov, regions=gb)
	values(gb) <- countmat
	return(gb)
}
# --------------------------------------------------------------------

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

		# pull back DEseq normalized counts - use for visualization
		# probably want to make the WIG files out of these
		#head( counts( cds, normalized=TRUE ) )

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
	# call patterns genome wide based on ab, ac, and bc test results
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

		# Trying to make this very easy - we just make the decision codes into strings, then all we have to do is match these strings with the decision table and get the pattern as a table lookup
		# These strings will be ordered the same as comparisions are ordered when returned by the get comps function
		str <- lapply(tests,function(x) x$code)
		str <- do.call(cbind,str)
		#str <- data.table(str)
		#str$id <- 1:nrow(str)
		#str <- str[,list(toString(.SD)),by=id]$V1
		# This is a potential bottleneck
		# The alternative would be to iterate over the rows of dc, check each against all of str (as  matrix), assign patterns that way
		#str <- apply(str,1,toString)
		# Now make the DecTable into codestrings
		#decide <- data.frame(pattern=dc$pattern,codestring=apply(dc[,-1,drop=F],1,toString))

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

		#doDecide <- function(x)
		#{
		#	dec <- t(t(str)==as.vector(t(x)))
		#	yespatt
		#}
		#dc2 <- split(dc[,-1,drop=F],1:nrow(dc))
		#decs <- lapply(dc2,doDecide)
		#mypatt$codestring <- str
		#mypatt$patt <- decide[match(mypatt$codestring,decide$codestring),]$pattern
		#mypatt[is.na(patt),]$patt <- "ambig"
		mypatt[patt=="000",]$patt <- "000or111"
		mypatt$patt <- as.character(mypatt$patt)
		mypatt
	}


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
		#s <- strsplit(x,"")[[1]]
		out <- methylaction:::testDESeq(bins.mat,samp$group,a=a,b=b,prefix=paste0("deseq-",a," vs ",b),sizefactors=sizefactors,ncore=ncore)
		return(out)
	}
	testres <- apply(comps$combos,1,dotest)
	names(testres) <- comps$strings

	if(!all(sapply(testres,nrow)==length(bins))){stop("Ran out of memory during testing, try reducing ncore")}

	# call patterns genome wide based on ab, ac, and bc test results
	callPatterns <- function(ab,ac,bc,cutoff=0.05)
	{
		# beacuse the bins all line up, we can work in column space
		ab$direction <- 0
		ab[ab$baseMeanA<ab$baseMeanB,]$direction <- 1
		ab[ab$baseMeanA>ab$baseMeanB,]$direction <- -1
		ac$direction <- 0
		ac[ac$baseMeanA<ac$baseMeanB,]$direction <- 1
		ac[ac$baseMeanA>ac$baseMeanB,]$direction <- -1
		bc$direction <- 0
		bc[bc$baseMeanA<bc$baseMeanB,]$direction <- 1
		bc[bc$baseMeanA>bc$baseMeanB,]$direction <- -1

		tests <- data.table(ab.p=ab$pval,ab.l2fc=ab$log2FoldChange,ab.sig=ab$pval<cutoff,ab.dir=ab$direction,ac.p=ac$pval,ac.l2fc=ac$log2FoldChange,ac.sig=ac$pval<cutoff,ac.dir=ac$direction,bc.p=bc$pval,bc.l2fc=bc$log2FoldChange,bc.sig=bc$pval<cutoff,bc.dir=bc$direction)

		# Fix for NA p-values
		# According to S. Anders: "Only in case 1 (zero counts in _all_ samples that are involved in the comparison), the p values is NA. This makes sense because if you do not observe anything from a gene you cannot say anything about it."
		# I will set these to be p-value 1, so they never come up as sig differences
		tests[is.na(ab.p),ab.sig:=FALSE]
		tests[is.na(ac.p),ac.sig:=FALSE]
		tests[is.na(bc.p),bc.sig:=FALSE]

		# filter using codes:
		# -1 = sig down
		# 0 = NS
		# 1 = sig up
		# could also just dispense with testing and join where directions agree, or use a much higher p-value cutoff like 0.5

		tests[,ab.code:=ab.dir]
		tests[ab.sig==FALSE,ab.code:=0]
		tests[,ac.code:=ac.dir]
		tests[ac.sig==FALSE,ac.code:=0]
		tests[,bc.code:=bc.dir]
		tests[bc.sig==FALSE,bc.code:=0]

		# table of all results that came out of the test
		tests[,codestring:=paste(ab.code,ac.code,bc.code,sep=",")]
		#table(tests$codestring)

		tests[,patt:="ambig"]
		tests[(ab.code==0)&(ac.code==1)&(bc.code==1),patt:="001"]
		tests[(ab.code==1)&(ac.code==1)&(bc.code==0),patt:="011"]
		tests[(ab.code==-1)&(ac.code==-1)&(bc.code==0),patt:="100"]
		tests[(ab.code==0)&(ac.code==-1)&(bc.code==-1),patt:="110"]
		tests[(ab.code==1)&(ac.code==0)&(bc.code==-1),patt:="010"]
		tests[(ab.code==-1)&(ac.code==0)&(bc.code==1),patt:="101"]
		tests[(ab.code==0)&(ac.code==0)&(bc.code==0),patt:="000or111"]
		return(tests)
	}

	# compare new vs old	
	#patt1 <- callPatterns(ab=testres[[1]],ac=testres[[2]],bc=testres[[3]],cutoff=stageone.p)
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
	# regions.gr

	# save to an Rd
	#save(bins.gr,testres,tests,regions.gr,file="testone.rd",compress=T)
	patterns <- bins
	values(patterns) <- patt

	test.one <- list(patterns=patterns, regions=regions)

	return(test.one)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
testTwo <- function(samp,cov,reads,stagetwo.method,regions,sizefactors,fragsize,winsize,anodev.p,adjust.var,post.p, fdr.filter, freq,ncore)
{
	message("Begin stage two testing")

	message("Recounting in regions with method=",stagetwo.method)

	if(stagetwo.method=="co")
	{
		recounts <- as.matrix(values(getCounts(samp=samp,reads=reads,ranges=regions,ncore=ncore)))
	} else if(stagetwo.method=="fc")
	{
		gdf <- data.frame(GeneID=1:length(regions),Chr=seqnames(regions),Start=start(regions),End=end(regions),Strand="+")
		recounts <- methylaction:::featureCountsDt(files=samp$bam, annot.ext=gdf, useMetaFeatures=F, allowMultiOverlap=T, read2pos="5", readExtension3=fragsize, strandSpecific="0", nthreads=ncore)$counts
		colnames(recounts) <- samp$sample
	} else if(stagetwo.method=="ac")
	{
		ba <- mclapply(cov, function(x) binnedAverage(bins=regions, numvar=x),mc.cores=1)
		recounts <- do.call(cbind,ba)
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

	#patt1 <- callPatterns2(testres[[1]], testres[[2]], testres[[3]], cutoff=post.p)
	#patt2 <- callPatterns2N(res=testres, cutoff=post.p)
	#table(patt1$patt==patt2$patt)
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

	#s1 <- sharpness("010")
	#s2 <- sharpnessN("010")
	#table(s1$frequent==s2$frequent)

	res <- lapply(patts,sharpnessN)
	dmrcalled <- do.call(rbind,unname(res))
	dmrcalled$dmrid <- 1:nrow(dmrcalled)

	# Make means columns
	counts.means <- do.call(cbind, lapply(colgroups, function(i) rowMeans(recounts.sig[,i,drop=F])))
	colnames(counts.means) <- paste0(unique(samp$group),".mean")

	# Final output!
	dmr <- cbind(dmr,counts.means)
	#dmr$f <- paste(dmr$a,dmr$b,dmr$c,sep=",")

	#test.two <- list()
	# test.two$ns
	#test.two$sig <- # GRanges with the normalized DEseq counts + deviances used from the ANODEV stage - also have these available in the $ns object
	test.two$dmr <- makeGRanges(dmr)
	test.two$dmrcalled <- makeGRanges(dmrcalled)

	return(test.two)
}

# given a directory with RData files saved by maPerm(), find them and merge the perms together, re-calculating the FDRs
maPermMerge <- function(dir=".")
{
	rds <- dir(dir,pattern="maperm_",full.names=T)

	mas <- list()
	idx <- 1
	for(i in rds)
	{
		message("Loading: ",i)
		load(i)
		mas[[idx]] <- maperm
		message("Found ",length(maperm)," perms")
		idx <- idx + 1
		rm(maperm)
	}

	# lapply had memory errors, sticking with one by one on the loop
	# isn't so bad in terms of runtime
	#mas <- mclapply(rds,function(x) {message("Loading: ",x);load(x);return(ma)},mc.cores=15)


	# now just grab all the perms
	mas <- do.call(c,mas)
	message("Loaded ",length(mas)," permutations from ",length(rds)," maPerm() .rd files")

	#load(rds[1])
	#ma$data$maperm <- mas
	#summary(ma$data$maperm)
	#message("Loaded ",length(ma$data$maperm)," perms")
	#tab <- methylaction:::maPermFdr(ma,recut.p=0.05)
	#ma$fdr <- tab
	return(mas)
}

# --------------------------------------------------------------------

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
methylaction <- function(samp, counts, bsgenome, fragsize, winsize, poifdr, stageone.p, joindist, anodev.p, post.p, minsize=150, ncore=1)
{
	# Assign groups from samp to a, b, c
	# Validate that we have 3 groups and each is replicated
	ngroups <- length(unique(samp$group))
	if(ngroups!=3){stop("samp must contain 3 groups - detected ",ngroups)}
	samp <- data.table(samp)
	if(sum(samp[,length(sample),by=group]$V1 > 1)!=3){stop("Each group must contain more than one replicate")}
	groupcodes <- c("a","b","c")
	names(groupcodes) <- unique(samp$group)
	samp$groupcode <- groupcodes[match(samp$group,names(groupcodes))]

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
	windows$filtered <- counts[!filter.pass & !zero]
	windows$signal <- counts[filter.pass]
	windows$signal.norm <- wingr[filter.pass]
	values(windows$signal.norm) <- normcounts

	# Testing Stage 1
	test.one <- methylaction:::testOne(samp=samp,bins=windows$signal,signal.norm=windows$signal.norm,chrs=unique(as.vector(seqnames(counts))),sizefactors=sizefactors,stageone.p=stageone.p,joindist=joindist,minsize=minsize,ncore=ncore)

	# Testing Stage 2 + Methylation Modelling
	test.two <- methylaction:::testTwo(samp=samp, regions=test.one$regions, sizefactors=sizefactors, bsgenome=bsgenome, fragsize=fragsize, anodev.p=anodev.p, post.p=post.p)

	# Output results
	ma <- list(opts=list(samp=samp,fragsize=fragsize,poifdr=poifdr,stageone.p=stageone.p,joindist=joindist,winsize=winsize,anodev.p=anodev.p,post.p=post.p,minsize=minsize,ncore=ncore),fdr.filter=fdr.filter, sizefactors=sizefactors, windows=windows, test.one=test.one, test.two=test.two)

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
	chunkids <- (seq(nrow(counts))-1) %/% 100000

	# need to use the same size factors for everything

	dotest <- function(rows)
	{
		# make master object container
		cds <- newCountDataSet(counts[rows,],groups)

		# estimate the size factors for the library
		#cds <- estimateSizeFactors(cds)
		sizeFactors(cds) <- sizefactors

		# pull back DEseq normalized counts - use for visualization
		# probably want to make the WIG files out of these
		#head( counts( cds, normalized=TRUE ) )

		# estimate the dispersions for each gene
		#message("Estimating dispersions")
		cds <- estimateDispersions(cds,fitType="local",sharingMode="gene-est-only")

		# perform the testing
		#message("Performing test")
		res <- nbinomTest(cds,a,b)
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
# Stage one testing
testOne <- function(samp,bins,signal.norm,chrs,sizefactors,stageone.p=0.05,minsize=150,joindist=200,ncore=3)
{
	message("Begin stage one testing")

	# Perform pairwise tests
	bins.mat <- as.matrix(values(bins))

	# a vs b
	todo <- c("ab","ac","bc")
	dotest <- function(x)
	{
		message("Testing ",x)
		s <- strsplit(x,"")[[1]]
		out <- methylaction:::testDESeq(bins.mat,samp$groupcode,a=s[1],b=s[2],prefix=paste0("deseq-",x),sizefactors=sizefactors,ncore=ncore)
		return(out)
	}
	testres <- lapply(todo,dotest)
	names(testres) <- todo

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
		tests[is.na(ab.p),ab.p:=1]
		tests[is.na(ac.p),ac.p:=1]
		tests[is.na(bc.p),bc.p:=1]

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
	
	patt <- callPatterns(testres$ab, testres$ac, testres$bc, cutoff=stageone.p)

	# Make means columns
	colgroups <- list(a=samp[samp$groupcode=="a",]$sample,b=samp[samp$groupcode=="b",]$sample,c=samp[samp$groupcode=="c",]$sample)
	counts.means <- do.call(cbind, lapply(colgroups, function(i) rowMeans(as.matrix(values(signal.norm[,i])))))
	colnames(counts.means) <- paste0(unique(samp$group),".mean")

	# Combine so we get these in the output
	patt <- data.table(counts.means,patt)

	# join adjacent equivalent patterns
	message("Reduction by pattern and disjoing regions")
	bins.gr <- bins
	values(bins.gr) <- NULL
	bins.gr$pattTestOne <- patt$patt
	
	# split out by pattern
	sigpatt <- bins.gr[!(bins.gr$pattTestOne %in% c("ambig","000or111"))]
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
testTwo <- function(samp,regions,bsgenome,sizefactors,fragsize,anodev.p, post.p)
{
	message("Begin stage two testing")

	# Recount from BAMs inside these regions (try something like easyRNAseq - see DESeq vingette)
	recounts <- Repitools::annotationCounts(x=samp$bam,anno=regions,seq.len=fragsize,up=0,down=0)
	colnames(recounts) <- samp$sample

	# Do an ANOVA-like framework
	# First test for ANY difference with an ANOVA-style DESeq test (ANODEV since these are GLM based negative binomial tests)
	# Adjust this for multiple testing
	# Move on the significant windows to post testing, adjust within each test only then
	# Call pattern based on these values
	testDESeqANODEV<- function(recounts,groups,ncore)
	{
		# chunk up data to allow parallel testing for large data sets
		chunkids <- (seq(nrow(recounts))-1) %/% 500

		dotest <- function(rows)
		{
			# make master object container
			cds <- newCountDataSet(recounts[rows,],groups)

			# estimate the size factors for the library
			# do we want to use size factors from previously for the whole library?
			sizeFactors(cds) <- sizefactors

			# pull back DEseq normalized counts - use for visualization
			# probably want to make the WIG files out of these
			#head( counts( cds, normalized=TRUE ) )

			# estimate the dispersions for each gene
			message("Estimating dispersions")
			cds <- estimateDispersions(cds,fitType="local",sharingMode="gene-est-only")

			# perform the testing
			# compare model w/ factor versus null model
			message("Fitting full GLM")
			resFull <- fitNbinomGLMs( cds, modelFormula=count~groups)
			message("Fitting reduced GLM")
			resReduced <- fitNbinomGLMs( cds, modelFormula=count~1)
			message("Performing test")
			res <- nbinomGLMTest(resFull, resReduced)
			return(res)
		}
		finalres <- mclapply(0:max(chunkids),function(x) {message("Testing chunk ",x," of ",max(chunkids));dotest(chunkids==x);},mc.cores=ncore)
		finalres <- do.call(c,finalres)
		return(finalres)
	}
	anodev <- testDESeqANODEV(recounts=recounts,groups=samp$group,ncore=ncore)
	anodev.padj <- p.adjust(anodev,method="fdr")
	#table(anodev<0.1)
	# Filter out if it does not pass the ANODEV

	regions$anodev.p <- anodev
	regions$anodev.padj <- anodev.padj

	# Get norm counts
	cds <- newCountDataSet(recounts,samp$group)
	sizeFactors(cds) <- sizefactors
	#values(regions) <- cbind(as.data.frame(values(regions)),counts(cds,normalized=TRUE))
	normcounts <- counts(cds,normalized=TRUE)

	anodev.keep <- (anodev.padj<anodev.p) & !(is.na(anodev.padj))

	test.two <- list()
	test.two$ns <- regions[!anodev.keep]
	test.two$ns.counts <- regions[!anodev.keep]
	values(test.two$ns.counts) <- normcounts[!anodev.keep,]

	regions.sig <- regions[anodev.keep]
	recounts.sig <- recounts[anodev.keep,]

	# a vs b
	todo <- c("ab","ac","bc")
	dotest <- function(x)
	{
		message("Testing ",x)
		s <- strsplit(x,"")[[1]]
		out <- methylaction:::testDESeq(recounts.sig,samp$groupcode,a=s[1],b=s[2],prefix=paste0("deseq-test2-",x),sizefactors,ncore)
		return(out)
	}
	testres <- mclapply(todo,dotest,mc.cores=1)
	names(testres) <- todo


	callPatterns2 <- function(ab,ac,bc,cutoff)
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

		# row by row p adjustment
		#p.adjust(c(ab$pval,ac$pval,bc$pval),method="fdr")

		# you set cutoff to 0.017 or 0.05/3 then you bonferroni adjust everything by row

		tests <- data.table(ab.p=ab$pval,ab.l2fc=ab$log2FoldChange,ab.sig=ab$pval<cutoff,ab.dir=ab$direction,ac.p=ac$pval,ac.l2fc=ac$log2FoldChange,ac.sig=ac$pval<cutoff,ac.dir=ac$direction,bc.p=bc$pval,bc.l2fc=bc$log2FoldChange,bc.sig=bc$pval<cutoff,bc.dir=bc$direction)


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

	patt <- callPatterns2(testres$ab, testres$ac, testres$bc, cutoff=post.p)
	
	test.two$sig <- regions[anodev.keep]
	values(test.two$sig) <- patt
	test.two$sig.counts <- regions[anodev.keep]
	values(test.two$sig.counts) <- normcounts[anodev.keep,]

	#as.data.frame(table(patt$patt))

	dmr <- data.frame(chr=seqnames(regions.sig),start=start(regions.sig),end=end(regions.sig),width=width(regions.sig),anodev.padj=regions.sig$anodev.padj,pattern=patt$patt)
	#write.table(dmr,file="tmp.bed",row.names=F,col.names=F,quote=F,sep="\t")
	
	#dpt <- fread("~/lustre/EurBLHsp0.1/methyl.20140701185731.sites.csv")
	#dmr$dpt <- makeGRanges(dmr) %over% makeGRanges(dpt)
	#dmr$url <- goldmine:::getBrowserURLs(makeGRanges(dmr),"hg19")
	#write.csv(dmr,file="tmp.csv",row.names=F)

	# want to look at DPT sites this one does not pick up
	#dpt$madmr <- makeGRanges(dpt) %over% makeGRanges(dmr)
	#dpt$url <- goldmine:::getBrowserURLs(makeGRanges(dpt),"hg19")
	#write.csv(dpt,file="tmp2.csv",row.names=F)

	#table(makeGRanges(dpt) %over% makeGRanges(dmr))
	#table(makeGRanges(dmr) %over% makeGRanges(dpt))

	#table(makeGRanges(dpt[pattern=="001",]) %over% makeGRanges(dmr[dmr$patt=="001",]))
	#table(makeGRanges(dmr) %over% makeGRanges(dpt))

	# Output list of sig DMRs
	
	# Run baymeth on these regions to add freq columns
	#message("Calculating CpG density for BayMeth")
	#gbA <- resize(regions.sig, 1, fix="center")
	#cpgdens <- cpgDensityCalc(gbA, organism=bsgenome, w.function="linear", window=700)

	# Make BayMethList object
	#message("Running BayMeth")
	#bml <- BayMethList(window=regions.sig,control=matrix(),sampleInterest=recounts.sig,cpgDens=cpgdens)

	# Set 0-offets because we have no M.SssI control
	#fOffset(bml) <- t(matrix(rep(0,ncol(recounts))))

	# For estimating parameters
	#bml <- empBayes(bml, ngroups = 100, ncomp = 1, maxBins = 50000, method="beta", ncpu=ncore, verbose=F)

	# Get methylation estimates
	#bml <- methylEst(bml, verbose=F, controlCI = list(compute = FALSE))

	# Obtain Final Values
	#me <- methEst(bml)$mean

	# Make frequency columns
	colgroups <- list(a=samp[samp$groupcode=="a",]$sample,b=samp[samp$groupcode=="b",]$sample,c=samp[samp$groupcode=="c",]$sample)

	#test.two$baymeth <- me

	#me.call <- me>0.75
	#meth.freq <- do.call(cbind, lapply(colgroups, function(i) rowSums(me.call[,i]>0.75)))
	#colnames(meth.freq) <- paste0(unique(samp$group),".freq")

	# Make means columns
	counts.means <- do.call(cbind, lapply(colgroups, function(i) rowMeans(recounts.sig[,i])))
	colnames(counts.means) <- paste0(unique(samp$group),".mean")

	# Final output!
	dmr <- cbind(dmr,counts.means)
	#dmr$f <- paste(dmr$a,dmr$b,dmr$c,sep=",")

	#test.two <- list()
	# test.two$ns
	#test.two$sig <- # GRanges with the normalized DEseq counts + deviances used from the ANODEV stage - also have these available in the $ns object
	test.two$dmr <- makeGRanges(dmr)

	return(test.two)
}
# --------------------------------------------------------------------

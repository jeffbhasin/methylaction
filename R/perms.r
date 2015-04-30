# Functions for performing permutations/bootstraps

# ====================================================================
# Exported Functions

# --------------------------------------------------------------------
#' Permute or Bootstrap DMR Detection
#'
#' Will perform permutations or bootstraps after methylaction() has been called. See also maPermMerge() and maPermFdr().
#' @param ma Output object from methylaction()
#' @param reads Preprocessed reads/fragments data from getReads()
#' @param nperms Number of permutation (or bootstrap) iterations to perform
#' @param combos A matrix of pre-set combinations to use. Useful for smaller sample sizes where there are only a limited number of possible combinations. The row should have ncols = number of samples, and each row represents a re-ordering of samples into groups (where the groups will be set in the order they appear in the "samp" given to methylaction())
#' @param save If true, save an RData file of the permutations. These can be merged using maPermMerge(). Useful for running permutations across multiple computers or in a cluster environment. If FALSE, will return the permutation results.
#' @param perm.boot If nperms > 0 and if TRUE, perform bootstrapping (sampling with replacement). Otherwise, perform permutations (sampling without permutations)
#' @param ncore Number of parallel processes to use
#' @return A list of DMRs arising from each permutation or bootstrap. Will save as an RData if save==TRUE.
#' @export
maPerm <- function(ma,reads,nperms,combos=NULL,save=T,perm.boot=F,ncore=1)
{
	perm.combo <- F
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

# --------------------------------------------------------------------
#' Merge permutations generated and saved by maPerm()
#'
#' The function maPerm() with save=TRUE will save permutations into RData files. Move all these RData files into the same directory, and indicate this directory using "dir". The output from maPermMerge() will be a list of DMRs detected in all permutations, and can be given to maPermFdr() to compute false discovery rates (FDRs).
#' @param dir Directory containing RData files saved by maPerm() where save was equal to TRUE
#' @return A list of DMRs arising from each permutation or bootstrap, merged across all the detected RData files.
#' @export
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

	mas <- do.call(c,mas)
	message("Loaded ",length(mas)," permutations from ",length(rds)," maPerm() .rd files")

	return(mas)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Compute false discovery rates (FDRs) from permuted DMR lists
#'
#' Using DMRs arising from permutations saved using maPerm(), compute FDRs for each pattern.
#' @param ma Output object from methylaction()
#' @param maperm Output object from maPerm()
#' @param recut.p ANODEV p-value cutoff to use
#' @param wide Produce a more readable "wide" format summary table, sorted by "frequent" FDR
#' @return A table of FDRs for each pattern, stratified by "frequent" status
#' @export
maPermFdr <- function(ma,maperm,recut.p=0.05,wide=FALSE)
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

	if(wide==TRUE)
	{
		tab <- fdr
		tab <- tab[tab$type %in% c("frequent","other"),c("pattern","type","nDMRs","FDRpercent")]
		tab <- tab[tab$pattern!="all",]
		c1 <- reshape::cast(tab,formula=pattern~type,value="nDMRs")
		c2 <- reshape::cast(tab,formula=pattern~type,value="FDRpercent")
		stopifnot(c1$pattern==c2$pattern)
		df <- data.frame(frequent_dmrs=c1$frequent,frequent_fdr=c2$frequent,other_dmrs=c1$other,other_fdr=c2$other)
		df <- cbind(pattern=c1$pattern,df)
		out <- df[order(df$frequent_fdr,decreasing=F),]
		out$frequent_fdr <- round(out$frequent_fdr,digits=1)
		out$other_fdr <- round(out$other_fdr,digits=1)
		fdr <- out
	}
	return(fdr)
}
# --------------------------------------------------------------------
# ====================================================================

# ====================================================================
# Internal Functions

# --------------------------------------------------------------------
# Draws permutation orders out of combination space by running a series of chooses
# This may be preferable to simply shuffling the sample label list because it ensures each permutation ends up with distinct configurations of groups
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

	# Want to pick one at random, then pick one of the next possible choices at random up to the number of permutations
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
# ====================================================================

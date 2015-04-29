# Functions for performing permutations/bootstraps

# ====================================================================
# Exported Functions

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

# --------------------------------------------------------------------
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
# ====================================================================

# ====================================================================
# Internal Functions

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
# ====================================================================

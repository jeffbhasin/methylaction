# Functions for summary and plotting after analysis has been run

# --------------------------------------------------------------------
#' Summary stats for a run of methylaction()
#'
#' Will return information about number of windows/regions that pass cutoffs at each stage of the analysis. Useful for paramater tuning.
#' @param ma Output object from methylaction()
#' @return A data.frame with the summary statistics.
#' @export
maSummary <- function(ma)
{
	per <- function(x,digits=2){round(x*100,digits)}

	# Initial Filtering
	df <- data.frame(stat="Window Size",count=ma$args$winsize,percent="",stringsAsFactors=F)
	wins <- length(ma$data$windows$zero) + length(ma$data$windows$filtered) + length(ma$data$windows$signal.norm)
	df <- rbind(df,c("Total Windows",wins,""))
	zero <- length(ma$data$windows$zero)
	filt <- length(ma$data$windows$filtered)
	signal <- length(ma$data$windows$signal.norm)
	df <- rbind(df,c("All Zero Windows (filtered)",zero,per(zero/wins)))
	df <- rbind(df,c("All Below FDR Windows (filtered)",filt,per(filt/wins)))
	df <- rbind(df,c("Signal Windows (move on to stage one)",signal,per(signal/wins)))

	# Stage One Testing
	owins <- length(ma$data$test.one$patterns)
	osig <- length(ma$data$test.one$patterns[!(ma$data$test.one$patterns$patt %in% c("000or111","ambig"))])
	ons <- length(ma$data$test.one$patterns[(ma$data$test.one$patterns$patt %in% c("000or111"))])
	oamb <- length(ma$data$test.one$patterns[(ma$data$test.one$patterns$patt %in% c("ambig"))])
	df <- rbind(df,c("Windows Tested in Stage One",owins,""))
	df <- rbind(df,c("Sig Pattern in Stage One",osig,per(osig/owins)))
	df <- rbind(df,c("Non-Sig Pattern in Stage One",ons,per(ons/owins)))
	df <- rbind(df,c("Ambig Pattern in Stage One",oamb,per(oamb/owins)))
	df <- rbind(df,c("Regions Formed By Joining Adjacent Patterns",length(ma$data$test.one$regions),""))

	# Stage Two Testing
	treg <- length(ma$data$test.two$sig) + length(ma$data$test.two$ns)
	df <- rbind(df,c("Regions Tested in Stage Two",treg,""))
	tsig <- length(ma$data$test.two$sig)
	df <- rbind(df,c("Regions That Pass ANODEV",tsig,per(tsig/treg)))

	tpsig <- length(ma$data$test.two$dmr[!(ma$data$test.two$dmr$pattern %in% c("000or111","ambig"))])
	tpns <- length(ma$data$test.two$dmr[(ma$data$test.two$dmr$pattern %in% c("000or111"))])
	tpamb <- length(ma$data$test.two$dmr[(ma$data$test.two$dmr$pattern %in% c("ambig"))])

	df <- rbind(df,c("ANODEV Sig with Sig Pattern",tpsig,per(tpsig/tsig)))
	df <- rbind(df,c("ANODEV Sig with Non-sig Pattern",tpns,per(tpns/tsig)))
	df <- rbind(df,c("ANODEV Sig with Ambig Pattern",tpamb,per(tpamb/tsig)))

	# DMRs
	df <- rbind(df,c("Total DMRs",tpsig,""))
	return(df)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Heatmap of the differentially methylated regions (DMRs) found by a run of methylaction()
#'
#' Will plot a heatmap of the DMRs based on the normalized read counts.
#' @param samp Description of samples from readSampleInfo()
#' @param ma Output object from methylaction()
#' @param pdf PDF file to output.
#' @return Saves plot to disk.
#' @export
maHeatmap <- function(ma,frequentonly=TRUE,bias=2,file)
{
	sites <- ma$dmr

	# Restrict to frequent only if requested
	if(frequentonly==TRUE)
	{
		sites <- sites[sites$frequent==TRUE]
	}

	# Extract counts matrix
	mat <- as.matrix(values(sites))
	samp <- ma$args$samp
	mat <- mat[,colnames(mat) %in% samp$sample]
	stopifnot(colnames(mat)==samp$sample)

	# Set scale inflection point based on POI cutoffs
	inflect <- sqrt(mean(ma$data$fdr.filter$cuts))

	# Get upper bound for saturation
	#upper <- sqrt(round(quantile(mat,0.9999),0))
	upper <- sqrt(max(mat))

	# Adjust by size
	#mat <- sqrt(mat/(width(sharpsites)/winsize))

	# set saturation point (hard in matrix)
	#mat[mat>upper] <- upper

	# Generate colors
	ncolors <- 100
	per <- round((inflect/upper)*100,0)/100
	pal1size <- round(per*ncolors)
	#pal1 <- colorRampPalette(c("#fff5f0","#fff5f0","#fff5f0","#fee0d2","#fc9272","#fb6a4a","#ef3b2c"),bias=1)(pal1size)
	#pal2 <- colorRampPalette(c("#ef3b2c","#cb181d","#a50f15","#67000d"),bias=1)(ncolors-pal1size)
	#pal1 <- colorRampPalette(c("#313695","#fee090","#f46d43"),bias=1)(pal1size)
	#pal2 <- colorRampPalette(c("#f46d43","#d73027","#a50026"),bias=1)(ncolors-pal1size)
	#pal1 <- colorRampPalette(c("#313695","#74add1","#fee090","#f46d43"),bias=1)(pal1size)
	#pal2 <- colorRampPalette(c("#f46d43","#d73027","#a50026"),bias=1)(ncolors-pal1size)
	pal1 <- colorRampPalette(rev(c("#313695","#4575b4","#fee090")),bias=bias)(pal1size)
	pal2 <- colorRampPalette(c("#fee090","#f46d43","#d73027","#a50026"),bias=bias)(ncolors-pal1size)
	cols <- c(rev(pal1),pal2)

	# Do sorting
	#ord <- order(ma$test.two$dmr$pattern)
	#cnt <- cnt[ord,]
	#allpatts <- apply(methylaction2:::getGroupPatterns(length(unique(ma$args$samp$group))),1,paste0,collapse="")
	#fac <- factor(sites$pattern,levels=allpatts)
	#cnt <- mat[order(sites$pattern,sites$anodev.padj),]
	cnt <- mat
	#cnt <- mat
	# Do transformations
	cnt <- sqrt(cnt/ma$args$winsize)

	# Do plotting
	#pdf(file=pdf,width=8,height=10.5)
	png(filename=file,width=2550,height=3300,res=300)
	sc <- c(benign="#4daf4a",low="#377eb8",high="#e41a1c")
	samp <- ma$args$samp
	sc <- unique(samp$color)
	names(sc) <- unique(samp$group)
	csc <- sc[match(samp$group,names(sc))]

	cs <- numeric(0)
	last <- ""
	for(i in 1:nrow(samp))
	{
		if(samp[i,]$group != last)
		{
			cs <- c(cs,i)
		}
		last <- samp[i,]$group
	}

	cs <- (cs-1)[-1] 

	rs <- match(unique(sites$pattern),sites$pattern)
	rs <- (rs-1)[-1] 
	gplots::heatmap.2(cnt,Colv=F,Rowv=F,trace="none",labRow=F,col=cols,ColSideColors=csc,colsep=cs, sepwidth=c(0.15,5),rowsep=rs)
	dev.off()

	message("Plot saved to ",file)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Clustering of samples based on read counts
#'
#' Will plot clustering based on some subset of read counts.
#' @param samp Description of samples from readSampleInfo()
#' @param ma Output object from methylaction()
#' @param mincv Minimum cv (cv=mean/sd) for a window to be included in the clustering.
#' @param type One of "mds" or "hca"
#' @param pdf PDF file to output.
#' @return Saves plot to disk.
#' @export
maClustering <- function(samp,ma,mincv,type,pdf)
{

}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Table of differentially methylated regions (DMRs) by pattern
#'
#' Will return a table of which patterns were detected and the number of DMRs in each.
#' @param samp Description of samples from readSampleInfo()
#' @param ma Output object from methylaction()
#' @return A data.frame with the summary statistics.
#' @export
maTable <- function(ma,recut.p=NULL,use.perms=NULL)
{
	library(reshape)

	# Add code to recompute FDR here for different ANODEV.P cutoffs
	# make observed event table
	if(!is.null(recut.p))
	{
		gettab <- function(dmrcalled)
		{
			dmrcalled <- dmrcalled[dmrcalled$anodev.padj<recut.p]
			dmrcalled$frequent <- factor(dmrcalled$frequent,levels=c(FALSE,TRUE))
			matab <- as.matrix(table(dmrcalled$pattern,dmrcalled$frequent))
			matab <- cbind(matab,rowSums(matab))
			colnames(matab) <- c("other","frequent","all")
			matab <- rbind(matab,colSums(matab))
			rownames(matab)[nrow(matab)] <- "all"
			return(matab)
		}

		# make expected event table - mean of all permutations
		realtab <- gettab(ma$dmr)
		realtab <- realtab[!rownames(realtab) %in% c("000or111","ambig"),]

		# Do we only want to do for a subset of perms?
		if(is.null(use.perms))
		{
			touse <- length(ma$data$maperm)
		} else {
			touse <- use.perms
			if(use.perms>length(ma$data$maperm)){stop("Asked to use more perms than are present in ma object")}
		}

		permtab <- lapply(1:touse,function(x){
			#message(x);
			gettab(ma$data$maperm[[x]]);
		})
		permtab <- lapply(permtab,function(x) x[!(rownames(x) %in% c("000or111","ambig")),])

		# Need to make sure they are ordered right before doing the division
		permtab <- lapply(permtab,function(x) x[rownames(realtab),])

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
		tab <- fdr
	} else
	{
			tab <- ma$fdr
	}

	#tab <- tab[tab$type %in% c("frequent","other"),c("pattern","type","nDMRs","FDRpercent")]
	#tab <- tab[tab$pattern!="all",]

	#c1 <- cast(tab,formula=pattern~type,value="nDMRs")
	#c2 <- cast(tab,formula=pattern~type,value="FDRpercent")
	#stopifnot(c1$pattern==c2$pattern)

	#df <- data.frame(frequent_dmrs=c1$frequent,frequent_fdr=c2$frequent,other_dmrs=c1$other,other_fdr=c2$other)

	#abc <- do.call(rbind,str_split(c1$pattern,""))[,-1]
	#colnames(abc) <- c("a","b","c")

	#df <- cbind(pattern=c1$pattern,abc,df)
	#df$pattern <- factor(df$pattern,levels=c("001","110","011","100","010","101"))
	#out <- df[order(df$pattern),]
	#return(out)

	return(tab)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Compare between various runs of methylaction()
#'
#' Reports shared and non-shared regions detected by different runs of methylaction(). Useful for comparing between studies/cohorts and paramater tuning.
#' @param malist List of output objects from methylaction(). Set the names() attribute of this list to unique descriptive names for each run.
#' @return A list containing both site by site comparisions and a data.frame summary table of total and shared sites for each pairwise comparision between runs.
#' @export
maCompare <- function(malist)
{

}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Write BED and BIGWIG files for normalized, filter-passed window count values
#'
#' Creates a BED file suitable for uploading as a custom track to the UCSC genome browser.
#' @param ma Output list from a run of methylaction()
#' @param path Folder to save the files in (defulat: current working directory)
#' @param bigwig Convert to BIGWIG files, requires wigToBigWig in $PATH (default: FALSE)
#' @param chrs
#' @param bsgenome
#' @param ncore
#' @return Writes BED file to disk.
#' @export
maTracks <- function(ma, path=".", bigwig=FALSE, chrs=NULL, bsgenome=NULL, ncore=NULL)
{
	if((bigwig==TRUE)&((is.null(chrs))|(is.null(bsgenome)))){stop("Must give chrs and bsgenome if bigwig=TRUE")}
	#if(is.null(counts)){stop("Must give list of data as counts argument")}

	# Make one GRanges for all bin coordinates
	#bins.gr <- suppressWarnings(do.call(c,lapply(names(counts$bins), function(x) GRanges(x,counts$bins[[x]]))))

	chrlens <- seqlengths(bsgenome)[chrs]

	signal <- ma$data$windows$signal.norm

	writeBed <- function(x)
	{
		message(paste(x,": Creating BedGraph",sep=""))
		filename <- file.path(path,paste(x,".bed",sep=""))
		values <- values(signal)[,x]
		bed <- data.frame(chr=seqnames(signal), start=as.integer(start(signal)-1), end=as.integer(end(signal)), value=values)
		bed <- bed[bed$value>0,]
		#trackheader <- paste("track","type=wiggle_0",paste("name=",x,sep=""),sep=" ")
		#write(trackheader,file=filename)
		write.table(bed, file=filename, append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

		if(bigwig==TRUE)
		{
			message(paste(x,": Converting to BigWig",sep=""))
			cs <- tempfile(x)
			write.table(data.frame(chrs,chrlens),col.names=F,row.names=F,file=cs, quote=F)
			cmd <- paste("wigToBigWig",filename,cs,file.path(path,paste(x,".bw",sep="")),sep=" ")
			message(cmd)
			system(cmd)
			file.remove(cs)
		}
	}
	# want to print all series of data in the file expect the bin metadata
	mclapply(ma$args$samp$sample, writeBed, mc.cores=ncore)
	NULL
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Write BED file of DMR regions
#'
#' Creates a BED file suitable for uploading as a custom track to the UCSC genome browser.
#' @param ma Output list from a run of methylaction()
#' @param file Name of BED file to create
#' @param dmr.only Don't report regions without a significant pattern (show up listed as NS, default: FALSE)
#' @return Writes BED file to disk.
#' @export
maBed <- function(ma, file, dmr.only=F)
{
	if(dmr.only==T)
	{
		# Only give ranges where we have sig pattern
		call.gr <- ma$dmr
		values(call.gr) <- NULL
		call.gr$call <- as.character(ma$dmr$pattern)
		call.gr$call[ma$dmr$sharp==F] <- paste0("111_",as.character(call.gr$call[ma$dmr$sharp==F]))
	} else
	{
		# Give classification for the entire genome

		# A: All zero bins
		#a <- ma$data$windows$zero
		#a <- reduce(a)
		#a$call <- "zero"

		# B: All below POI bins
		#b <- ma$data$windows$filtered
		#b <- reduce(b)
		#b$call <- "filtered"

		# C: NS in stage one
		c <- ma$data$test.one$patterns[!(ma$data$test.one$patterns %over% ma$data$test.one$regions)]
		values(c) <- NULL
		c <- reduce(c)
		c$call <- "NS1"

		# D: NS in stage two ANODEV
		d <- ma$data$test.two$ns
		values(d) <- NULL
		d <- reduce(d)
		d$call <- "NS2"

		# E: NS in stage two pairwise
		e <- ma$data$test.two$sig[!(ma$data$test.two$sig %over% ma$dmr)]
		values(e) <- NULL
		e <- reduce(e)
		e$call <- "NS3"

		# F: Sig DMR called w/ pattern
		f <- ma$dmr
		values(f) <- NULL
		f$call <- as.character(ma$dmr$pattern)
		f$call[ma$dmr$sharp==F] <- paste0("111_",as.character(f$call[ma$dmr$sharp==F]))

		# Merge
		#call.gr <- c(a,b,c,d,e,f) # obvious that anything not in there is zero or filtered
		call.gr <- c(d,e,f)
	}

	# Make BED
	bed <- data.frame(chr=seqnames(call.gr),start=as.integer(start(call.gr))-1,end=as.integer(end(call.gr)),name=call.gr$call, score=0, strand="+", thickStart=as.integer(start(call.gr))-1, thickEnd=as.integer(end(call.gr)))

	patts <- unique(bed$name)
	pal <- colorRampPalette(RColorBrewer::brewer.pal(8,"Set2"))(length(patts))
	#pal <- sample(pal)

	mycols <- data.frame(name=patts,color=pal)
	mycols$itemRgb <- apply(col2rgb(mycols$color),2,function(x) paste(x,collapse=","))
	bed$itemRgb <- mycols[match(bed$name,mycols$name),]$itemRgb

	# Write BED
	if(file.exists(file)){file.remove(file)}
	header <- "track name=\"methylaction\" itemRgb=\"On\" visibility=\"pack\""
	cat(header, '\n',file=file)
	write.table(bed,file=file, col.names=F, row.names=F, sep=" ", quote=F,append=T)

	NULL
}
# --------------------------------------------------------------------

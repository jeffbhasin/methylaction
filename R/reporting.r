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
	df <- data.frame(stat="Window Size",count=ma$opts$winsize,percent="",stringsAsFactors=F)
	wins <- length(ma$windows$zero) + length(ma$windows$filtered) + length(ma$windows$signal)
	df <- rbind(df,c("Total Windows",wins,""))
	zero <- length(ma$windows$zero)
	filt <- length(ma$windows$filtered)
	signal <- length(ma$windows$signal)
	df <- rbind(df,c("All Zero Windows (filtered)",zero,per(zero/wins)))
	df <- rbind(df,c("All Below FDR Windows (filtered)",filt,per(filt/wins)))
	df <- rbind(df,c("Signal Windows (move on to stage one)",signal,per(signal/wins)))

	# Stage One Testing
	owins <- length(ma$test.one$patterns)
	osig <- length(ma$test.one$patterns[!(ma$test.one$patterns$patt %in% c("000or111","ambig"))])
	ons <- length(ma$test.one$patterns[(ma$test.one$patterns$patt %in% c("000or111"))])
	oamb <- length(ma$test.one$patterns[(ma$test.one$patterns$patt %in% c("ambig"))])
	df <- rbind(df,c("Windows Tested in Stage One",owins,""))
	df <- rbind(df,c("Sig Pattern in Stage One",osig,per(osig/owins)))
	df <- rbind(df,c("Non-Sig Pattern in Stage One",ons,per(ons/owins)))
	df <- rbind(df,c("Ambig Pattern in Stage One",oamb,per(oamb/owins)))
	df <- rbind(df,c("Regions Formed By Joining Adjacent Patterns",length(ma$test.one$regions),""))

	# Stage Two Testing
	treg <- length(ma$test.two$sig) + length(ma$test.two$ns)
	df <- rbind(df,c("Regions Tested in Stage Two",treg,""))
	tsig <- length(ma$test.two$sig)
	df <- rbind(df,c("Regions That Pass ANODEV",tsig,per(tsig/treg)))

	tpsig <- length(ma$test.two$dmr[!(ma$test.two$dmr$pattern %in% c("000or111","ambig"))])
	tpns <- length(ma$test.two$dmr[(ma$test.two$dmr$pattern %in% c("000or111"))])
	tpamb <- length(ma$test.two$dmr[(ma$test.two$dmr$pattern %in% c("ambig"))])

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
maHeatmap <- function(samp,ma,pdf)
{

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
maTable <- function(samp, ma)
{

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


maTracks <- function(ma, path="", bigwig=FALSE, chrs=NULL, bsgenome=NULL, ncore=NULL)
{
	if((bigwig==TRUE)&((is.null(chrs))|(is.null(bsgenome)))){stop("Must give chrs and bsgenome if bigwig=TRUE")}
	#if(is.null(counts)){stop("Must give list of data as counts argument")}

	# Make one GRanges for all bin coordinates
	#bins.gr <- suppressWarnings(do.call(c,lapply(names(counts$bins), function(x) GRanges(x,counts$bins[[x]]))))

	signal <- ma$windows$signal

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
	mclapply(samp$sample, writeBed, mc.cores=ncore)
	NULL
}
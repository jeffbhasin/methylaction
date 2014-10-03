# Functions for data input/output/preprocessing

# --------------------------------------------------------------------
#' Load a CSV containing required information about each sample
#'
#' The CSV file must contain the following columns: "sample" - unique sample IDs, "group" - group IDs, "bam" - path to BAM file containing aligned reads for the sample. Columns with other names will be ignored. Note that in subsequent reporting of pattern strings (where each digit represents a group), digits for each group will be ordered in the order they first appear in this samplesheet.
#' @param file Path to the CSV samplesheet to open. Must contain the columns described above.
#' @param colors Vector of colors (one for each group) in same order as groups appear in the sample file. These will be uniform colors used in the plotting functions for these groups. Give colors as hex codes. If none provided they will be auto-selected with RColorBrewer. If there is a column named "color" in the CSV, then this will always be used.
#' @return A data.frame of the samplesheet that will be valid input to the other function's "samp" arguments.
#' @export
readSampleInfo <- function(file=NULL,colors=NULL)
{
	if(is.null(file)){stop("Must provide a path to a file to open.")}
	if(!file.exists(file)){stop(paste("Could not open: ",file,sep=""))}
	samplesheet <- read.csv(file, header=TRUE, stringsAsFactors=FALSE)
	if(sum(c("sample","group","bam") %in% colnames(samplesheet))!=3){stop("CSV must contain columns: sample, group, bam")}
	nsamp <- nrow(samplesheet)
	ngroup <- length(unique(samplesheet$group))
	message(paste("Found ", nsamp, " samples in ", ngroup, " groups",sep=""))
	grouporder <- unique(samplesheet$group)
	samplesheet$group <- factor(samplesheet$group,levels=grouporder,ordered=TRUE)
	message(paste("Group order for pattern and contrast digits will be: ", paste(grouporder,collapse=", ", sep=""), sep=""))
	message("Number of samples in each group:")
	print(summary(samplesheet$group))
	samplesheet <- samplesheet[order(samplesheet$group,samplesheet$sample),]

	if(is.null(colors) & sum(colnames(samplesheet)=="color")==0)
	{
		message("Auto-picking group colors from RColorBrewer")
		if(ngroup<=9)
		{
			colors <- RColorBrewer::brewer.pal(ngroup,"Set1")
		} else
		{
			colors <- colorRampPalette(brewer.pal(9,"Set1"))(ngroup)
		}
		samplesheet$color <- colors[as.numeric(samplesheet$group)]
	} else if(!is.null(colors))
	{
		if(length(colors)!=ngroup){stop("colors vector must equal number of groups in file")}
		samplesheet$color <- colors[as.numeric(samplesheet$group)]
	}

	samplesheet	
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Produce a GRanges with counts of overlapping reads for a set of ranges
#'
#' Will return a GenomicRanges object for non-overlapping windows genome-wide for the genome given as a bsgenome object for the chromosomes given in chrs. The values() of the GRanges will contain a table of counts for each sample at each window.
#' @param samp Sample data.frame from readSampleInfo()
#' @param bsgenome
#' @param ranges
#' @param fragsize Average fragment length from the sequencing experiment. Reads will be extended up to this size when computing coverage.
#' @param winsize Size of the non-overlapping windows.
#' @return A GenomicRanges object with values() containing a table of counts for each sample at each window.
#' @export
getCounts <- function(samp, reads, bsgenome=NULL, ranges=NULL, chrs=NULL, winsize=50, ncore=1)
{
	# Either give bsgenome to automatically tile the genome, or give a pre-defined ranges set
	if(is.null(bsgenome)&is.null(ranges)){stop("Need to give either bsgenome (count in genome-wide non-overlapping windows for the given genome) or ranges (count in pre-defined GenomicRanges object).")}
	if(!is.null(bsgenome)&!is.null(ranges)){stop("Please give either bsgenome or ranges, not both.")}

	# Make GRanges of non-overlapping windows accross the genome
	if(!is.null(bsgenome))
	{
		message("Generating Window Positions")
		gb <- Repitools::genomeBlocks(Hsapiens,chrs=chrs,width=winsize)
	} else
	{
		gb <- ranges
	}

	# Do countOverlaps
	docount <- function(x)
	{
		message("Counting Overlaps for ",x)
		return(countOverlaps(gb,reads[[x]],type="any",ignore.strand=F))
	}
	if(ncore>5){ncore <- 5}
	counts <- mclapply(names(reads),docount,mc.cores=ncore)
	names(counts) <- samp$sample
	counts <- do.call(cbind,counts)

	# Build GRanges object to return so counts are never separated from their coordinates
	message("Building counts GRanges")
	values(gb) <- counts
	message("Counts GRanges Size in Memory=",format(object.size(gb),units="auto"))
	gc()
	return(gb)
}
# --------------------------------------------------------------------

# --------------------------------------------------------------------
#' Store GenomicRanges of BAM reads in a list
#'
#' Will return a list of GenomicRanges objects
#' @param samp Sample data.frame from readSampleInfo()
#' @param bsgenome
#' @param chrs
#' @param fragsize Average fragment length from the sequencing experiment. Reads will be extended up to this size when computing coverage.
#' @return A list of GenomicRanges objects.
#' @export
getReads <- function(samp, bsgenome, chrs, fragsize, ncore)
{
	message("Reading BAMs from disk with ",ncore," concurrent processes")
	# Validate BAM existence
	# Check that BAM files exist
	check <- sapply(samp$bam,file.exists)
	if(sum(check)!=nrow(samp)){stop(paste0("BAM file(s) could not be found: ",toString(samp$bam[check==F])))}

	# Validate asked for chrs are in the bsgenome
	# Check that BAM has the chrs asked for
	if(sum(chrs %in% seqlevels(bsgenome)) != length(chrs)){stop(paste0("Could not find chrs: ",toString(chrs[!(chrs %in% seqlevels(bsgenome))]), " in given bsgenome"))}

	# Use the index so we don't bother reading in from the unaligned chrs
	which.gr <- GRanges(chrs,IRanges(1,seqlengths(bsgenome)[chrs]))

	bam2gr <- function(bampath)
    {
    	message("Reading BAM file: ",bampath)
    	# What: fields to read in (column filtering)
    	# Flag: records to read in (row filtering)
    	# Which: what sequences must be overlaped (chr/pos filtering)
    	param <- Rsamtools::ScanBamParam(what=character(), which=which.gr, flag=Rsamtools::scanBamFlag(isUnmappedQuery=FALSE))
    	bam.ga <- GenomicAlignments::readGAlignments(bampath, param = param)
    	# Not using which
	    #bam.ga3 <- bam.ga2[seqnames(bam.ga2) %in% chrs]
    	#param <- Rsamtools::ScanBamParam(what=character(), flag=scanBamFlag(isUnmappedQuery=FALSE))
    	#bam.ga2 <- GenomicAlignments::readGAlignments(bampath, param = param)
		bam.gr <- as(bam.ga, "GRanges")
		message("Extending read length to ",fragsize)
		bam.gr <- resize(bam.gr, fragsize)
		seqlevels(bam.gr) <- chrs
		return(bam.gr)
	}
	reads <- mclapply(samp$bam,bam2gr,mc.cores=ncore)
	names(reads) <- samp$sample
	message("Read GRanges Size in Memory=",format(object.size(reads),units="auto"))
	return(reads)
}
# --------------------------------------------------------------------

# How much memory does any R object x take up?
sizein <- function(x)
{
	format(object.size(x),units="auto")
}

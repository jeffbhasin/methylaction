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
#' Produce a GRanges of windows and their read counts
#'
#' Will return a GenomicRanges object for non-overlapping windows genome-wide for the genome given as a bsgenome object for the chromosomes given in chrs. The values() of the GRanges will contain a table of counts for each sample at each window.
#' @param samp Sample data.frame from readSampleInfo()
#' @param chrs Which chromosomes to consider.
#' @param fragsize Average fragment length from the sequencing experiment. Reads will be extended up to this size when computing coverage.
#' @param winsize Size of the non-overlapping windows.
#' @return A GenomicRanges object with values() containing a table of counts for each sample at each window.
#' @export
getCounts <- function(samp, chrs, fragsize, winsize, regions=NULL, ncore)
{
	# Check that BAM files exist
	check <- sapply(samp$bam,file.exists)
	if(sum(check)!=nrow(samp)){stop(paste0("BAM file(s) could not be found: ",toString(samp$bam[check==F])))}

	# Read in BAMs as GAlignments
	# If we want to provide a unified filtering when reading in from BAM, this is the place
	message("Reading data from BAM files using ",ncore," processes")
	bams.ga <- mclapply(samp$bam,GenomicAlignments::readGAlignments,use.names=F,mc.cores=ncore)
	names(bams.ga) <- samp$sample

	# Check that BAM has the chrs asked for
	if(sum(chrs %in% seqlevels(bams.ga[[1]])) != length(chrs)){stop(paste0("Could not find chrs: ",toString(chrs[!(chrs %in% seqlevels(bams.ga[[1]]))]), " in BAM header"))}

	# Filter to only requested chrs
	#message("Filtering reads to only those on chrs: ",toString(chrs))
	#bams.ga <- mclapply(bams.ga,function(x) x[seqnames(x) %in% chrs],mc.cores=ncore)

	# Get chromosome lengths (were extracted from the BAM header by readGAlignments)
	chrlens <- seqlengths(bams.ga[[1]])[chrs]

	# Convert to GRanges
	bams.gr <- mclapply(bams.ga,function(x) as(x, "GRanges"),mc.cores=ncore)
	rm(bams.ga)
	gc()

	# Extend reads
	message("Extending reads to length ",fragsize,"bp")
	bams.gr2 <- mclapply(bams.gr, function(x) resize(x, fragsize),mc.cores=ncore)
	rm(bams.gr)
	gc()

	# Generate windows
	# Otherwise use a provided anno so we can use this function for the stageTwo re-counting
	if(is.null(regions))
	{
		message("Calculating Windows")
		gb <- Repitools::genomeBlocks(chrlens,chrs=chrs,width=winsize)
	} else
	{
		gb <- regions
	}

	# Count inside windows
	message("Counting")
    counts <- mclapply(bams.gr2, function(x) matrix(countOverlaps(gb, x)),mc.cores=ncore)
    counts2 <- do.call(cbind,counts)
    colnames(counts2) <- names(counts)

    # Return GRanges Result
    counts.gr <- gb
    values(counts.gr) <- counts2
    gc()
    return(counts.gr)
}
# --------------------------------------------------------------------
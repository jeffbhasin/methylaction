# Functions for data input/output/preprocessing

# --------------------------------------------------------------------
#' Load a CSV containing required information about each sample
#'
#' The CSV file must contain the following columns: "sample" - unique sample IDs, "group" - group IDs, "bam" - path to BAM file containing aligned reads for the sample. Columns with other names will be ignored. Note that in subsequent reporting of pattern strings (where each digit represents a group), digits for each group will be ordered in the order they first appear in this samplesheet.
#' @param file Path to the CSV samplesheet to open. Must contain the columns described above.
#' @param colors Vector of colors (one for each group) in same order as groups appear in the sample file. These will be uniform colors used in the plotting functions for these groups. Give colors as hex codes. If none provided they will be auto-selected with RColorBrewer.
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

	if(is.null(colors))
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
	} else
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
#' @param begenome A b-string genome (bsgenome) object from Bioconductor for the genome of interest.
#' @param chrs Which chromosomes to consider.
#' @param fragsize Average fragment length from the sequencing experiment. Reads will be extended up to this size when computing coverage.
#' @param winsize Size of the non-overlapping windows.
#' @return A GenomicRanges object with values() containing a table of counts for each sample at each window.
#' @export
getCounts <- function(samp, bsgenome, chrs, fragsize, winsize)
{
	bins.gr <- genomeBlocks(bsgenome,chrs=chrs,width=winsize)
	message("Generated ",prettyNum(length(bins.gr),big.mark=",",scientific=F)," non-overlapping windows of width ", winsize,"bp genome-wide")
	counts.genome <- annotationCounts(x=samp$bam,anno=bins.gr,seq.len=fragsize,up=0,down=0)
	colnames(counts.genome) <- samp$sample
	values(bins.gr) <- counts.genome
	return(bins.gr)
}
# --------------------------------------------------------------------

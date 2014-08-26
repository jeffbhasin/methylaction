# Functions for summary and plotting after analysis has been run

# --------------------------------------------------------------------
#' Summary stats for a run of methylaction()
#'
#' Will return information about number of windows/regions that pass cutoffs at each stage of the analysis. Useful for paramater tuning.
#' @param samp Description of samples from readSampleInfo()
#' @param ma Output object from methylaction()
#' @return A data.frame with the summary statistics.
#' @export
maSummary <- function(samp, ma)
{

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

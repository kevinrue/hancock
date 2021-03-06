% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/learn-methods.R
\name{learnMarkersByPositiveProportionDifference}
\alias{learnMarkersByPositiveProportionDifference}
\title{Identify Markers by Largest Difference of Detection Rate in Clusters}
\usage{
learnMarkersByPositiveProportionDifference(se, cluster.col,
  assay.type = "counts", threshold = 0, n = Inf, min.diff = 0.1,
  min.prop = 0.1, diff.method = c("min", "mean", "median", "max"))
}
\arguments{
\item{se}{An object of class inheriting from "\code{\link{SummarizedExperiment}}".}

\item{cluster.col}{Name of a column in \code{colData(se)} that contains
a factor indicating cluster membership for each column (i.e. sample) in \code{se}.}

\item{assay.type}{A string specifying which assay values to use, e.g., \code{"counts"} or \code{"logcounts"}.}

\item{threshold}{Value \emph{above which} the marker is considered detected.}

\item{n}{Maximal number of markers allowed for each signature.}

\item{min.diff}{Minimal difference in detection rate between the target cluster
and the summarized detection rate in any other cluster (in the range 0-1).
See argument \code{diff.method}.}

\item{min.prop}{Minimal proportion of samples in the target cluster where the combined set of markers is detected.}

\item{diff.method}{Method to contrast the detection rate in the target cluster to that of all other clusters.
See Details section.}
}
\value{
A collection of signatures as a "\code{\link{Sets}}" object.
}
\description{
This function computes the detection rate of each feature in each cluster.
For each cluster, it ranks all the features by decreasing difference between
the detection rate in the target cluster, and the detection rate in all other clusters.
The function can limit results up to \code{n} markers for each cluster.
}
\details{
\code{diff.method} controls how the detection rate in all clusters \emph{other than the target one} are summarized before comparison with the detection in the target cluster.
It is possible to rank features using the minimal (\code{"min"}), \code{"mean"}, \code{"median"}, or maximal (\code{"max"}) difference between the detection rate in the target cluster and those of all other clusters.
}
\examples{
# Example data ----
library(SummarizedExperiment)
nsamples <- 100
u <- matrix(rpois(20000, 1), ncol=nsamples)
rownames(u) <- paste0("Gene", sprintf("\%03d", seq_len(nrow(u))))
colnames(u) <- paste0("Cell", sprintf("\%03d", seq_len(ncol(u))))
se <- SummarizedExperiment(assays=list(counts=u))

colData(se)[, "cluster"] <- factor(sample(head(LETTERS, 3), ncol(se), replace=TRUE))

# Example usage ----

baseset <- learnMarkersByPositiveProportionDifference(se, cluster.col="cluster")

relations(baseset)
}
\seealso{
\code{\link{learnSignatures}}.
}
\author{
Kevin Rue-Albrecht
}

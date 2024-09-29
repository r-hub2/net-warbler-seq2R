\name{find.points}
\alias{find.points}

\title{
Simple method to detect compositional changes in genomic sequences
}
\description{
\code{find} is used to detect changes at genomic sequences composition. The method is based on fitting nonparametric models by using local linear kernel smoothers.
}
\usage{
find.points(x,kbin= 300, p= 3, bandwidth=-1, weights= 1, nboot=100, kernel="gaussian",
n.bandwidths= 20, seed = NULL, ...)
}

\arguments{
  \item{x}{
Sequences in binary system (by using \code{change.binary} function previously) are to be analyzed from.
}
  \item{kbin}{
The number of binning nodes over which the function is to be estimated.
}
  \item{p}{
Degree of the polynomial. By default \code{p=3}.
}
  \item{bandwidth}{
The kernel bandwidth or smoothing parameter. Large values of bandwidth make smoother estimates, smaller values of bandwidth make less smooth estimates. The default \code{h=-1} is a bandwidth compute by cross validation.
}
  \item{weights}{
Weights.
}
 \item{nboot}{
Number of bootstrap repeats.
}
\item{kernel}{
Character which denotes the kernel function (a symmetric density). By default \code{kernel = "gaussian"}, this is, the Gaussian density function. Also, other types of kernel functions can be used: Epanechnikov and triangle, \code{kernel="Epanech"} and \code{kernel="triang"}, respectively.
}
\item{n.bandwidths}{Number that it will be used to calculate the grid of bandwidths in a range between 0 and 1. In this grid, it will be selected the optimum bandwidth by cross-validation.If the optimum bandwidth value is close to 0, we will obtain rough estimates; when it is close to 1, we will obtain smooth estimates.
}
\item{seed}{Seed to be used in the bootstrap procedure.}
\item{\ldots}{Other options.}

}
\details{
For each genomic sequence the AT and CG skews profiles were calculated as \eqn{A vs. T = (A-T)/(A+T)} and \eqn{C vs. G = (C-G)/(C+G)}.
}
\value{
The function computes and returns a list of short information for a fitted \code{change.points} object.
  \item{Number of A-T base pairs}{The returned value is the total nucleotide  (adenine and thymine) contained at the sequence analyzed.}
 \item{Number of C-G base pairs }{In this case, the returned value is the sum of cytosine and guanine contained at the sequence.}
  \item{Number of binning nodes}{The number of binning nodes over which the function is to be estimated.}
  \item{Number of bootstrap repeats}{Number of bootstrap repeats.}
  \item{Bandwidth}{Value of the kernel bandwidth or smoothing parameter used in the fitting for A vs. T and C vs. G.}
   \item{Exists any critical point}{Emphasize if there is or not any critical.}
}

\author{
Nora M. Villanueva and Marta Sestelo.
}

\references{N. M. Villanueva, M. Sestelo, M. M. Fonseca and J. Roca-Pardinas (2023).
seq2R: An R package to detect change points in DNA sequences. Mathematics, 11 (10), 2299.}


\examples{
library(seq2R)


mtDNAhuman <- read.genbank("NC_012920")
DNA<- transform(mtDNAhuman)
seq1<-find.points(DNA)
seq1
}

\keyword{change points}


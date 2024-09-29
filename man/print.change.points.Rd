\name{print.change.points}
\alias{print.change.points}

\title{
Short find.points summary
}
\description{
\code{\link{find.points}} summary.}
\usage{
\method{print}{change.points}(x=model,...)
}

\arguments{
  \item{x}{\code{change.points} object.}
\item{\ldots}{Other options.}
}

\value{
The function computes and returns a list of short information for a fitted \code{change.points} object.
  \item{Number of A-T base pairs}{The returned value is the total nucleotide  (adenine and thymine) contained in the sequence analyzed.}
 \item{Number of C-G base pairs }{In this case, the returned value is the sum of cytosine and guanine contained at the sequence.}
  \item{Number of binning nodes}{The number of binning nodes over which the function is to be estimated.}
  \item{Number of bootstrap repeats}{Number of bootstrap repeats.}
  \item{Bandwidth}{Value of the Kernel bandwidth or smoothing parameter used in the fiiting for A vs. T and C vs. G.}
   \item{Exists any critical point}{Emphasize if there is or not any critical.}
}

\author{
Nora M. Villanueva and Marta Sestelo.
}
\note{
See details in \code{\link{find.points}}.
}


\examples{
library(seq2R)
mtDNAhuman <- read.genbank("NC_012920")
DNA <- transform(mtDNAhuman)
seq1 <- find.points(DNA)
seq1
}

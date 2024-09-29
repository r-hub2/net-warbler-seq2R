\name{critical}
\alias{critical}

\title{
Critical points (maxima and minima)
}
\description{
Function that maximizes or minimizes the first derivative of the model obtained with \code{\link{find.points}} function. Also, it is included their 95\% confidence intervals.
}
\usage{
critical(model, base.pairs = NULL)
}

\arguments{
  \item{model}{\code{change.points} object.}
  \item{base.pairs}{Character string for A vs. T or C vs G.}
}

\value{
The returned list has two component (\code{$AT}, \code{$CG}). Both of them containing  a matrix  with values about their critical points (maxima and minima), lower  and upper 95\% confidence intervals.
  \item{AT }{Critical points for AT.}
  \item{CG }{Critical points for CG.}
}

\author{
Nora M. Villanueva and Marta Sestelo.
}


\references{N. M. Villanueva, M. Sestelo, M. M. Fonseca and J. Roca-Pardinas (2023).
seq2R: An R package to detect change points in DNA sequences. Mathematics, 11 (10), 2299.}

\examples{
library(seq2R)

mthumanDNA <- read.genbank("NC_012920")
DNA <- transform(mthumanDNA)
seq1 <- find.points(DNA, nboot = 10)


critical(seq1,base.pairs="CG")


critical(seq1,base.pairs="AT")
}

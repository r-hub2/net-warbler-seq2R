\name{transform}
\alias{transform}

\title{
Convert biological sequences into binary code
}
\description{
Biological sequences are categorical variables. With this function the four nucleotides are coded with two bits, 0 and 1 (binary numeral system) for being used by almost all modern computers.
}
\usage{
transform(x)
}

\arguments{
  \item{x}{
 The object obtained with \code{\link{read.all}} or \code{\link{read.genbank}} functions is the argument required for \code{transform}. The nature of the sequence is DNA. Sequences are returned as a vector of single characters.
}
}

\value{
The returned list has two component (\code{$AT}, \code{$CG}). Both of them containing  a matrix  with values about their critical points (maximum and minimum), and their lower  and upper  95\% confidence intervals.
  \item{AT}{Variable A and T with binary system.}
  \item{CG}{Variable C and G with binary system. }
}

\author{
Nora M. Villanueva and Marta Sestelo.
}

\examples{
library(seq2R)

mthumanDNA <- read.genbank("NC_012920")
DNA <- transform(mthumanDNA)
DNA
}


\keyword{binary}
\keyword{AT}
\keyword{GC}

\name{plot.change.points}
\alias{plot.change.points}

\title{Visualization of change.points objects
}
\description{
Useful for drawing the estimation and first derivative of the skew profile.
}
\usage{

\method{plot}{change.points}(x = model, y = NULL, base.pairs = NULL, der = NULL,
xlab = "x", ylab = "y", col = "black", CIcol = "black", main = NULL, type = "l",
CItype = "l", critical = FALSE, CIcritical = FALSE,ylim=NULL,...)
}

\arguments{
  \item{x}{
 \code{change.points} object.
}
  \item{y}{\code{NULL}
}
\item{base.pairs}{Character string about the skew profile for A vs. T or C vs. G.
}
  \item{der}{ Number which determines inference process to be drawing into the plot. By default \code{der} is \code{NULL}. If it is \code{0}, the plot represents the initial estimate. If \code{der} is \code{1}, the first derivative is plotted.
}
  \item{xlab}{
Title for x axis.
}
  \item{ylab}{
Title for y axis.
}
  \item{col}{
A specification for the default plotting color.
}
  \item{CIcol}{
A specification for the default confidence intervals plotting color.
}
  \item{main}{
An overall title for the plot.
}
  \item{type}{
 Type of plot should be drawn. Possible types are, \code{p} for points, \code{l} for lines, \code{o} for overplotted, etc. See details in \code{?par}.
}
  \item{CItype}{ Type of plot should be drawn for confidence intervals. Possible types are, \code{p} for points, \code{l} for lines, \code{o} for overplotted, etc. See details in \code{?par}.
}
\item{critical}{ A logical value. If \code{TRUE} (not by default), the critical points are drawn into the plot.
}
\item{CIcritical}{ A logical value. If \code{TRUE} (not by default), the 95\% confidence intervals for the critical points are drawn into the plot.
}
\item{ylim}{The y limits of the plot.
}
  \item{\dots}{
Other options.
}
}

\value{
  Simply produce a plot.
}

\author{
Nora M. Villanueva and Marta Sestelo.
}

\examples{
library(seq2R)


#mtDNAhum <- read.genbank("NC_012920")
data(mtDNAhum)
DNA <- transform(mtDNAhum)
seq1 <- find.points(DNA)



plot(seq1,der=0,base.pairs="CG",CIcritical=TRUE,ylim=c(0.08,0.67))
plot(seq1,der=1,base.pairs="CG",CIcritical=TRUE,ylim=c(-0.0005,0.00045))
abline(h=0)



plot(seq1,critical=TRUE, CIcritical=TRUE)
}

\keyword{plot}

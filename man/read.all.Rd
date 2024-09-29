\name{read.all}
\alias{read.all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Read FASTA and GBK formatted files
}
\description{
Read nucleic acid sequences from a file in FASTA or GBK format.
}
\usage{
read.all(file = system.file(""), seqtype = "DNA")
}

\arguments{
  \item{file}{ The name of the file which the sequences in FASTA or GBK format are to be
  read from.
}
  \item{seqtype}{The nature of the sequence. Nowadays only \code{DNA}, in further updates it will be possible to use for  different type of sequences.}
}
\details{
Fasta is a widely used format in molecular biology. Sequence in FASTA format starts with a single-line description, distinguished by a greater-than `>' symbol,   followed by sequence data on the next lines.

'GenBank' format files have the extension GBK, by convention. Files contain fields with different types of information well-labeled. The header of the file has information describing the sequence, such as its  type, shape, length  and source. The features of the genome sequence follow the header, and include protein translations. The DNA sequence is the last element of the file, which ends with (and must include) a soluble slash. Complete genomes in this format are available at the \url{https://ftp.ncbi.nlm.nih.gov/genbank/}.
}
\value{

  \item{Sequence }{The returned list has a component  \code{Sequence} containing the DNA sequence taken from the field ``ORIGIN'' in GenBank. The sequence is a vector of single characters.}
  \item{Locus or accession}{the returned list has a component \code{Locus/Accession} containing the names of the locus or accession number taken from the field ``LOCUS'' or ``ACCESSION'' in 'GenBank'. Also, return sequence length.}
}

\author{Nora M. Villanueva and Marta Sestelo.
}

\examples{
library(seq2R)
data(mtDNAhum)
\dontrun{
data<-read.all("file.fasta")
data<-read.all("file.gbk")
}
}

\keyword{fasta}
\keyword{gbk}
\keyword{sequence}

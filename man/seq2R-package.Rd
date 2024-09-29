\name{seq2R-package}
\alias{seq2R-package}
\alias{seq2R}
\docType{package}
\title{
Simple method to detect compositional changes in genomic sequences
}
\description{\code{seq2R} is just a shortcut for "sequence to \code{R}". This software is useful for loading .fasta or .gbk files, and for recovering sequences from GenBank database.
This package allows  to detect differences or asymmetries based on nucleotide composition by using local linear kernel smoothers. Also, it is possible to draw inference about critical points (i. e. maximum or minimum points) related with the derivative curves. Additionally, bootstrap methods have been used  for estimating confidence intervals and speed computational techniques (binning techniques) have been implemented in "seq2R".
}

\author{
Nora M. Villanueva and Marta Sestelo.

Maintainer: Nora M. Villanueva <nmvillanueva@uvigo.es>
}
\references{
Efron, B. (1979). Bootstrap methods: another look at the jackknife. Annals of Statistics, 7:1-26.

Efron, E. and Tibshirani, R. J. (1993). An introduction to the Bootstrap. Chapman and Hall, London.

Fan, J. and Marron, J.S. (1994). Fast implementation of nonparametric curve estimators. Journal of Computational and Graphical Statistics, 3:35-56.

Hastie, T. and  Tibshirani, R. (1990). Generalized additive models. Chapman and Hall, London.

Stone, C. J. (1977). Consistent nonparametric regression. Annals of Statistics, 5: 595-620.

Villanueva, N. M., Sestelo, M.,  Fonseca, M. M. and Roca-Pardinas, J.  (2023).
seq2R: An R package to detect change points in DNA sequences. Mathematics, 11 (10), 2299.

Wand, M. P. and Jones, M. C. (1995). Kernel Smoothing. Chapman and Hall, London.
}


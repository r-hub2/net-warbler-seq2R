#-------------------------------------------------------
# demo package "seq2R", realease = 1.0 data= 2012-07-31
#-------------------------------------------------------


# 1. Example Human mitochondrial DNA
#------------------------------------------------------
library(seq2R)
mtDNAhuman <- read.genbank("NC_012920")
mtDNAhuman

DNA <- transform(mtDNAhuman)
#DNA

seq1 <- find(DNA)
seq1


# Specifying the base pairs CG,
#the estimates, their first derivative and
# 95% confidence intervals of the critical points

plot(
  seq1,
  der = 0,
  base.pairs = "CG",
  CIcritical = TRUE,
  ylim = c(0.08, 0.67)
)
plot(
  seq1,
  der = 1,
  base.pairs = "CG",
  CIcritical = TRUE,
  ylim = c(-0.0005, 0.00045)
)
abline(h = 0)

# Critical points CG
critical(seq1, base.pairs = "CG")


# Specifying the base pairs AT,
#the estimates, their first derivative and
# 95% confidence intervals of the critical points

plot(seq1, der = 0, base.pairs = "AT", CIcritical = TRUE)
plot(seq1, der = 1, base.pairs = "AT", CIcritical = TRUE)
abline(h = 0)


# Critical points AT
critical(seq1,base.pairs="AT")



# Visualization of change.points objects
plot(seq1, critical = TRUE, CIcritical = TRUE)

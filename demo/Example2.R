#-------------------------------------------------------
# demo package "seq2R", realease = 1.0 data= 2012-07-31
#-------------------------------------------------------


# 2. Example Nematode DNA
#------------------------------------------------------
library(seq2R)
nematode <- read.genbank("NC_013253")
nematode

nem <- transform(nematode)
#nem

seq2 <- find(nem, kbin = 450, nh = 10)
seq2

# Specifying the base pairs AT,
#the estimates, their first derivative and
# 95% confidence intervals of the critical points

plot(seq2, der = 0, base.pairs = "AT", CIcritical = TRUE)
plot(seq2, der = 1, base.pairs = "AT", CIcritical = TRUE,
     ylim = c(-0.0002,0.0002))
abline(h=0)


# Critical points AT
critical(seq2, base.pairs = "AT")



# Visualization of change.points objects
plot(seq2, critical = TRUE, CIcritical = TRUE)

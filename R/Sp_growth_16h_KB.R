growth <- read.csv("Sp_growth_16h_KB.csv", header=T)

seqe <- seq(1, 95)
seqe <- rep(seqe, each=8)

library(dplyr)

growth <- mutate(growth,
                 Time = seqe)

summary(growth)

library(ggplot2)

par(mfrow=c(3,3))

growthPC <- ggplot(growth, aes(y=PC, x=Time)) +
  geom_point(aes(col=Well)) +
  xlab(label = "Read number/10 mins") +
  ylab(label = "OD of Pseudomonas corrugata")

growthPC

growthSR <- ggplot(growth, aes(y=SR, x=Time)) +
  geom_point(aes(col=Well)) +
  xlab(label = "Read number/10 mins") +
  ylab(label = "OD of Stenotropomonas rhisophilia")

growthSR

growthAA <- ggplot(growth, aes(y=AA, x=Time)) +
  geom_point(aes(col=Well)) + 
  xlab(label = "Read number/10 mins") +
  ylab(label = "OD of Achromobacter agilis")

growthAA

growthVG <- ggplot(growth, aes(y=VG, x=Time)) +
  geom_point(aes(col=Well)) +
  xlab(label = "Read number/10 mins") +
  ylab(label = "OD of Variovorax guangxiensis")

growthVG

growthOD <- ggplot(growth, aes(y=OD, x=Time)) +
  geom_point(aes(col=Well)) +
  xlab(label = "Read number/10 mins") +
  ylab(label = "OD of Ochrobactrum daejonense")

growthOD

install.packages("mtcars")
library(mtcars)

growth2 <- 

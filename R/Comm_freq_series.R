#Load packages

library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(ggridges)
library(stringr)
library(lme4)

#Read in data

Commdata <- read.csv("Comm_freq_series.csv", stringsAsFactors = FALSE)

#Convert colony counts to CFUs/mL. 10^-5, 60ul plated
#Density_CFU = CFU per mL = Count*Dilution*ML

# ML = dilf fac to convert from plated vol to 1 mL (1000 µl / plated_vol)
# Dilution = serial dilution plated
# Count = number counted on plate
# Plated_vol = µl plated

#Formula- counts(x) * 10^5 * (1000/60)

Commdata <- mutate(Commdata,
                   CFU.mL = Count * Dilution * Dil_fact)

# All data now converted to CFUs/mL

#For graph purposes, may be easier to visualise as relative proportions

#sum of total CFUs/mL in each community
Comm_total_CFUs <- aggregate(Commdata$CFU.mL, by = list(Commdata$Comm, Commdata$Time), sum)

colnames(Comm_total_CFUs) <- c("Comm", "Time", "sum_cfu")

Comm_total_CFUs <- mutate(Comm_total_CFUs,
                          Comm =  substr(Comm, 5, 6))

Commdata <- mutate(Commdata,
                   Comm =  substr(Comm, 5, 6))

#merge

Commdata2 <- merge(Commdata, Comm_total_CFUs, all.x=T)

#SUCCESS- that only took 2 hours...
#tidy dataframe

Commdata_seq <- Commdata2[c("Comm", "Time", "Sp", "CFU.mL", "sum_cfu")]

Commdata_seq <- mutate(Commdata2, 
                       Rel_freq = CFU.mL / sum_cfu)

#plot species abundance over time for each community

Commplot <- ggplot(Commdata_seq, aes(x= forcats::fct_reorder(Time, Sp), y= Rel_freq, group= Sp), na.rm=T) + 
  geom_point(position = position_dodge(0.1), aes(colour=Sp)) +
  geom_line(position = position_dodge(0.1), aes(colour=Sp)) +
  scale_color_discrete("Sp") +
  facet_wrap(~Comm) +
  xlab("Time point") +
  ylab("Species relative frequency")

Commplot

#variation between communities by species frequency

Sp_plot <- ggplot(Commdata_seq, aes(x= Comm, y= Rel_freq), na.rm=T) +
  geom_col(aes(fill=Sp)) +
  facet_wrap(~Time) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Community replicate") +
  ylab("Species relative frequency")

Sp_plot

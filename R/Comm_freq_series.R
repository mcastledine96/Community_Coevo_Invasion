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

names(Commdata)

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
  ylab("Species proportion in community")

Commplot

Commplot_CFU <- ggplot(Commdata_seq, aes(x= forcats::fct_reorder(Time, Sp), y= CFU.mL, group= Sp), na.rm=T) + 
  geom_point(position = position_dodge(0.1), aes(colour=Sp)) +
  geom_line(position = position_dodge(0.1), aes(colour=Sp)) +
  scale_color_discrete("Sp") +
  facet_wrap(~Comm) +
  xlab("Time point") +
  ylab("Species abundance (CFU/mL)") +
  scale_y_continuous(limits = c(0, 5.0e+08))

Commplot_CFU

#variation between communities by species frequency

#mosaic plot
Sp_plot <- ggplot(Commdata_seq, aes(x= Comm, y= Rel_freq), na.rm=T) +
  geom_col(aes(fill=Sp)) +
  facet_wrap(~Time) +
  scale_fill_brewer(palette = "Set3") +
  xlab("Community replicate") +
  ylab("Species proportion in community")

Sp_plot

#bar plot

Sp_barplot <- ggplot(Commdata_seq, aes(x=Comm, y= Rel_freq), na.rm=T) +
  geom_col(aes(fill=factor(Comm))) +
  facet_wrap(~Sp + Time) + 
  xlab("Community replicate") +
  ylab("Species proportion in community") +
  theme(axis.text.x = element_text(size=7))

Sp_barplot

#boxplot

Sp_boxplot <- ggplot(Commdata_seq, aes(x=Sp, y=Rel_freq), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Time) +
  xlab("Species") +
  ylab("Species proportion in Community")

Sp_boxplot

Sp_boxplot_CFU <- ggplot(Commdata_seq, aes(x=Sp, y=CFU.mL), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Time) +
  xlab("Species") +
  ylab("Species abundance (CFU/mL)") +
  scale_y_continuous(limits = c(0, 5.0e+08))

Sp_boxplot_CFU #boxplot with CFU/mL instead of relative frequency

#Community productivity over time

prodplot <- ggplot(Comm_total_CFUs, aes(x=forcats::fct_reorder(Time, Comm), y= sum_cfu, group=Comm), na.rm=T) + 
  geom_point(position = position_dodge(0.1), aes(colour=Comm)) +
  geom_line(position = position_dodge(0.1), aes(colour=Comm)) +
  scale_color_discrete("Comm") +
  xlab("Time point") +
  ylab("Community productivity (total CFU/mL)") +
  scale_y_continuous(limits = c(0, 1.1e+09))

prodplot

prodplot2 <- ggplot(Comm_total_CFUs, aes(x=Time, y= sum_cfu, group=Comm), na.rm=T) + 
  geom_point(position = position_dodge(0.1), aes(colour=Comm)) +
  geom_line(position = position_dodge(0.1), aes(colour=Comm)) +
  facet_wrap(~Comm) +
  xlab("Time point") +
  ylab("Community productivity (total CFU/mL)") +
  scale_y_continuous(limits = c(0, 1.1e+09))

prodplot2

prodplot3 <- ggplot(Comm_total_CFUs, aes(x=Time, y= sum_cfu), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  scale_color_discrete("Comm") +
  xlab("Time point") +
  ylab("Community productivity (total CFU/mL)") +
  scale_y_continuous(limits = c(0, 9e+08))

prodplot3

#set scale to remove weird outlier at T1 which squidged the other lines. 

mod3 <- lm(CFU.mL ~ Sp + as.factor(Comm), data = Commdata_seq)
mod4 <- lm(CFU.mL ~ Sp * as.factor(Comm), data = Commdata_seq)
anova(mod3, mod4)
#the difference between species CFU/mL varies between communities significantly- we'd expect this considering the species interact within communities to determine their CFU/mL

#try fitting a mixed effect model to see whether species relative abundance is significantly different to one another with community ID and time as random effects (groups non-independent)

mod <- lmer(CFU.mL ~ Sp + (1 | Time) + (1 | Comm), data=Commdata_seq)
mod2 <- lmer(CFU.mL ~ 1 + (1 | Time) + (1 | Comm), data=Commdata_seq)
anova(mod, mod2)
#significant effect of Sp

lsmeans::lsmeans(mod, pairwise ~ Sp)
#significant comparisons between all species

Species <- cbind(Commdata_seq$CFU.mL[Commdata_seq$Sp == "Pc"], Commdata_seq$CFU.mL[Commdata_seq$Sp == "Od"], Commdata_seq$CFU.mL[Commdata_seq$Sp == "Vg"], Commdata_seq$CFU.mL[Commdata_seq$Sp == "Sr"], Commdata_seq$CFU.mL[Commdata_seq$Sp == "Aa"])

summary(manova(Species ~ Time + Comm, data=Commdata_seq)) #need to wrangle dataframe to make col lengths same and informative

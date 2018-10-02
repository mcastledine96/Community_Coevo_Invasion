library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(ggridges)
library(stringr)
library(lme4)
install.packages("lsmeans")
library(lsmeans)

#load un data
la.dat <- read.csv("Exp1_data_300918.csv", stringsAsFactors = FALSE)
names(la.dat)

#calculate CFU/mL of resident community and invading species
la.dat <- mutate(la.dat,
                 T1.res.cfu = Res.count * Dilution * ML * Froz.fact,
                 T1.inv.cfu = Inv.count * Dilution * ML * Froz.fact)

#Calculate relative fitness based on Malthusian fitness parameters
la.dat <- mutate(la.dat,
                 rel.fitness = log(T1.inv.cfu/T0.inv)/log(T1.res.cfu/T0.res))

#plot relative fitness of coevolved and coevolved invaders, separating out individual species
ggplot(la.dat, aes(x=coev, y=rel.fitness), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Inv.sp) +
  scale_y_continuous(limits = c(0.7, 2)) +
  ylab(label = "Relative invader fitness") +
  xlab(label = "Invader coevolutionary history to community") +
  scale_x_discrete(labels = c("Non-coev.", "Coev."))

#Same as above, ignoring species ID
ggplot(la.dat, aes(x=coev, y=rel.fitness), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  scale_y_continuous(limits = c(0.7, 2)) +
  ylab(label = "Relative invader fitness") +
  xlab(label = "Invader coevolutionary history to community") +
  scale_x_discrete(labels = c("Non-coev.", "Coev."))

#axis labels

labels <- c(expression(atop(Achromobacter,
                       agilis)),
            expression(atop(Ochrobactrum,
                       daejonense)),
            expression(atop(Pseudomonas,
                            corrugata)),
            expression(atop(Stenotropomonas,
                            rhisophilia)),
            expression(atop(Variovorax,
                            guangxiensis)))

#plot invasion rel fitness based on species ID
ggplot(la.dat, aes(x=forcats::fct_reorder(Inv.sp, coev), y=rel.fitness), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_point(position = position_jitter(0.2), aes(shape = coev), size =1.5) +
  scale_shape_manual(values=c(15,16), name = "Invader type", labels = c("Non-coev.", "Coev.")) +
  scale_y_continuous(limits = c(0.7, 2)) +
  scale_x_discrete(labels = labels) +
  xlab(label = "Invading species") +
  ylab(label = "Relative invader fitness")
  
?scale_shape_manual

#model data to test significance using lmer- generalised linear mixed effect model

#remove problem data
la.dat2 <- filter(la.dat, ! Treatment == 15)

mod1 <- lmer(rel.fitness ~ Inv.sp * coev + (1|Res.com.ID), la.dat2)
mod2 <- lmer(rel.fitness ~ Inv.sp + coev + (1|Res.com.ID), la.dat2)
anova(mod1, mod2)
#no interaction between species ID and coevolutionary history

mod3 <- lmer(rel.fitness ~ Inv.sp + (1|Res.com.ID), la.dat2)
mod4 <- lmer(rel.fitness ~ 1 + (1|Res.com.ID), la.dat2)
anova(mod3, mod2)
#Coevolutionary history not significant 
anova(mod3, mod4)
#Invasive species ID significant

lsmeans::lsmeans(mod3, pairwise ~ Inv.sp)

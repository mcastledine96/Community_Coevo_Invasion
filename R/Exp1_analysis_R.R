library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(ggridges)
library(stringr)
library(lme4)
library(emmeans)

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
  xlab(label = "Invader Type") +
  scale_x_discrete(labels = c("Ancestral", "Random", "Coev."))

#ignoring species ID
ggplot(la.dat, aes(x=coev, y=rel.fitness), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  scale_y_continuous(limits = c(0.7, 2)) +
  ylab(label = "Relative invader fitness") +
  xlab(label = "Invader type") + 
  scale_x_discrete(labels = c("Ancestral", "Random", "Coev."))

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
  scale_shape_manual(values=c(1,15,16), name = "Invader type", labels = c("Ancestral","Random", "Coev.")) +
  scale_x_discrete(labels = labels) +
  xlab(label = "Invading species") +
  ylab(label = "Relative invader fitness")
  

#model data to test significance using lmer- generalised linear mixed effect model accounting for experimental design- invasion of species in resident community within block

#remove problem data
la.dat2 <- filter(la.dat, ! Treatment == 15)
la.dat2 <- filter(la.dat2, ! Treatment == 40)


## Does invasive species ID and coev affect community productivity?

la.dat4 <- mutate(la.dat, 
                  Comm.summ.T1 = T1.res.cfu + T1.inv.cfu,
                  Comm.summ.T0 = T0.res + T0.inv)

la.dat4 <- mutate(la.dat4,
                          Comm.prod = log(Comm.summ.T1/Comm.summ.T0))

ggplot(la.dat4, aes(x=coev, y=Comm.prod), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Inv.sp) +
  ylab(label = "Community productivity log(CFU/mL)") +
  xlab(label = "Invader type") +
  scale_x_discrete(labels = c("Ancestral","Random", "Coev."))

la.dat4 <- filter(la.dat4, ! Treatment == 10)

prod.mod1 <- lmer(Comm.prod ~ Inv.sp * coev + (1|Block), la.dat4)
prod.mod2 <- lmer(Comm.prod ~ Inv.sp + coev + (1|Block), la.dat4)
anova(prod.mod1, prod.mod2)
#no interaction
prod.mod3 <- lmer(Comm.prod ~ Inv.sp + (1|Block), la.dat4)
anova(prod.mod2, prod.mod3)
#effect of coevolution
prod.mod4 <- lmer(Comm.prod ~ coev + (1|Block), la.dat4)
anova(prod.mod2, prod.mod4)
#no effect of invasive species ID
prod.mod5 <- lmer(Comm.prod ~ 1 + (1|Block), la.dat4)
anova(prod.mod4, prod.mod5)
#sig effect of coev

emmeans(prod.mod4, pairwise ~ coev)

##rerun rel fitness models

#Linear mixed effect model of relative fitness (invader - community) against coevolutionary history of invader and invader species identity with random effect of resident community within block

la.max <- lmer(rel.fitness ~ Inv.sp * coev + (1|Block), data = la.dat2)

#model checks
par(mfrow = c(2,2))
plot(la.max)
plot(la.max, resid(., scaled=TRUE) ~ fitted(.) | Block, abline = 0)

#looks homoscedastic
par(mfrow = c(1,1))
qqnorm(resid(la.max))
#some deviations but overall looks normal enough

#Does invasion sp ID and coev interact?
la.max.noint <- lmer(rel.fitness ~ Inv.sp + coev + (1|Block), data = la.dat2)
anova(la.max, la.max.noint) #no sig effect of including interaction (P = 0.9885)

summary(la.max.noint)

#Effect of coev
la.max.nc <- update(la.max.noint,~. - coev)
anova(la.max.nc, la.max.noint)
#Coevo sig. P < 0.001

#Effect of sp ID
la.nsp <- update(la.max.noint,~. - Inv.sp)
anova(la.max.noint, la.nsp)
#Sp ID sig. Drop coev
summary(la.max.noint)

la.null <- lmer(rel.fitness ~ 1 + (1|Block), data = la.dat2)
la.nulls <- lmer(rel.fitness ~ Inv.sp + (1|Block), data = la.dat2)
anova(la.null, la.nulls) #effect of species sig- p<0.001
la.nulls3 <- lmer(rel.fitness ~ coev + (1|Block), data = la.dat2)
anova(la.null, la.nulls3) #effect of coev non-sig on own. P = 0.056

#multi-model inference
la.max <- lmer(rel.fitness ~ Inv.sp * coev + (1|Block), data = la.dat2, na.action = "na.fail")

d.mod <- dredge(la.max)
d.mod #multi-model inference- most informative model contains Inv.sp only


#relevant packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(ggridges)
library(stringr)
library(lme4)
library(lsmeans)
library(nlme)

#load data
la.ca <- read.csv("LA_AN_MO_20.11.18.CSV", header = T)

#calculate CFU/mL of replacement species- *2 to account for plating from frozen. 

la.ca <- mutate(la.ca,
                sp.cfu = Count * ML * Dil * 2)

#calculate fitness of replacement species

la.ca <- mutate(la.ca,
                sp.fit = log(sp.cfu/T0))

##redo null- based on 1.1 null done at same time. 

null <- read.csv("Null_trial_2.csv", header = T)

null <- mutate(null, 
               null.cfu = Count * Dil * ML * 2)

null <- mutate(null,
               null.fit = log(null.cfu/T0))

dat2 <- mutate(la.ca,
                sp.rel.fit2 = sp.fit/N.fit2)

#
write.csv(dat2, "Exp1_1_1.FINAL.csv", row.names = F)
#plot this

ggplot(dat2, aes(x = Sp, y = sp.rel.fit2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Treat) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = 
          element_blank())+
  theme_bw() +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  ylab("Species fitness relative to same species local to community") +
  xlab("Species") +
  scale_y_continuous(limits = c(0.6, 1.75))

##take a look at fitness itself

sepdat <- read.csv("Exp1_1.1.FINAL.fitsep.csv", header = T)

ggplot(sepdat, aes(x = Treat, y = Fit)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Sp, scales = "free") +
  scale_x_discrete(labels = c("Ancestral", "Random", "Coev.")) +
  ylab("Species fitness log(CFU/mL)") +
  xlab("Treatment")

##different way of calculating relative fitness for local adaptation experiment- rel fitness compared to self when in own community versus in diff community

self.fit <- read.csv("LA_self_fit.csv", header = T)

self.fit <- mutate(self.fit,
                   rel.fit = sp.fit/self.null)

ggplot(self.fit, aes(x = Sp, y = rel.fit)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = 
          element_blank())+
  theme_bw() +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  ylab("Fitness of random species relative to when in sympatric community") +
  xlab("Species")

t.test(self.fit$rel.fit[self.fit$Sp == "OD"], mu = 1) #sig
t.test(self.fit$rel.fit[self.fit$Sp == "SR"], mu = 1) #non-sig
t.test(self.fit$rel.fit[self.fit$Sp == "AA"], mu = 1) #non-sig
t.test(self.fit$rel.fit[self.fit$Sp == "PC"], mu = 1) #sig
t.test(self.fit$rel.fit[self.fit$Sp == "VG"], mu = 1) #non-sig

##trial analyses

com.mod <- lmer(Fit ~ Treat * Sp + (1|Block), data = sepdat)
com.mod.noin <- lmer(Fit ~ Treat + Sp + (1|Block), data = sepdat)


anova(com.mod, com.mod.noin) #no sig interaction

#test effect of species
com.mod.sp <- lmer(Fit ~ Treat + (1|Block), data = sepdat)
anova(com.mod.sp, com.mod.noin) #sig effect of species

#test effect of Treatment
com.mod.tr <- lmer(Fit ~ Sp + (1|Block), data = sepdat)
anova(com.mod.noin, com.mod.tr) #sig effect of treatment on species fitness

emmeans(com.mod.noin, pairwise ~ Treat+Sp)



## Try linear mixed effect model with just local adaptation and null data

la.dat <- read.csv("LA_dat.csv", header = T)

la.mod <- lmer(sp.fit ~ Treat + Sp + (1|Pair) + (1|Block), data = la.dat)
summary(la.mod)

la.modin <- lmer(sp.fit ~ Treat * Sp + (1|Pair) + (1|Block), data = la.dat)

anova(la.mod, la.modin) #no interaction

la.mod.nt <- lmer(sp.fit ~ Sp + (1|Pair) + (1|Block), data = la.dat)

anova(la.mod, la.mod.nt) #sig effect of Treatment

la.mod.ns <- lmer(sp.fit ~ Treat + (1|Pair) + (1|Block), data = la.dat)

anova(la.mod, la.mod.ns) #sig effect of species

emmeans::emmeans(la.mod, pairwise ~ Sp + Treat)
#something not quite right here

ggplot(la.dat, aes(x = Treat, y = sp.fit, group = Pair), na.rm = T) +
  geom_point(position = position_dodge(0.2), aes(colour=Pair)) +
  geom_line(position = position_dodge(0.2), aes(colour=Pair)) +
  facet_wrap(~Sp, scales = "free")


##              junk code probably               ##

#calculate fitness of replacement species relative to same species sympatric to community - not sure if this is the best way of doing this. Null = rel fitness = 1

la.ca <- mutate(la.ca,
                sp.rel.fit = sp.fit/N.fit)

#plot this

ggplot(la.ca, aes(x = Sp, y = sp.rel.fit)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Treat) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = 
          element_blank())+
  theme_bw() +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  ylab("Species fitness relative to same species local to community") +
  xlab("Species")

##take a look at fitness itself

ggplot(ca.sep, aes(x = Treat, y = sp.fit)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Sp, scales = "free") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = 
          element_blank())+
  theme_bw()


##model- will have to use fitness lines as separate- relative fitness doesn't allow for comparison between fitness between being in sympatric versus allopatric community. 
#Will have to chat to smart people about the unbalanced design- take means of the ancestral and la replicates for comparison? Until then:

ca.sep <- read.csv("1_1.modelling.csv", header = T)

com.mod <- lmer(sp.fit ~ Treat * Sp + (1|Block), data = ca.sep)
com.mod.noin <- lmer(sp.fit ~ Treat + Sp + (1|Block), data = ca.sep)

anova(com.mod, com.mod.noin) #no sig interaction

#test effect of species
com.mod.sp <- lmer(sp.fit ~ Treat + (1|Block), data = ca.sep)
anova(com.mod.sp, com.mod.noin) #sig effect of species

#test effect of Treatment
com.mod.tr <- lmer(sp.fit ~ Sp + (1|Block), data = ca.sep)
anova(com.mod.noin, com.mod.tr) #sig effect of treatment on species fitness

emmeans(com.mod.noin, pairwise ~ Treat) 
#sig differences in species fitness between ancestral and null and ancestral and allopatric species treatment. No sig difference between allopatric species and null. Suggests species are adapted to being in a community but this is not localised. 
#Focus on same-species pairwise comparisons between treatments: This effect is being driven by PC and AA- sig comparisons between AN and NU. 

## effects on resident species

res <- read.csv("Exp1_1_1_res.fit.var.csv", header = T)

res <- mutate(res,
              res.cfu = Count * ML * Dil * 2)

res <- mutate(res,
              res.fit = log(res.cfu/T0))

ggplot(res, aes(x = Sp.rep, y = res.fit)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Sp * Treat, scales="free")




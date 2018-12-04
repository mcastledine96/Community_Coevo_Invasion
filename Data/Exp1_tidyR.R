#Tidy up analyses from experiment 1

#Load in packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(ggridges)
library(stringr)
library(lme4)
library(lsmeans)
library(nlme)

#load in wrangled data from experiment 1

la.dat.expx <- read.csv("Exp_wrangled.csv", header = T)

#For now, remove some problem data with NAs or negative fitness
la.dat.exp <- filter(la.dat.expx, ! Treatment == 15)
la.dat.exp <- filter(la.dat.exp, ! Treatment == 40)

#Linear mixed effect model of relative fitness (invader - community) against coevolutionary history of invader and invader species identity with random effect of resident community within block

la.max <- lmer(rel.fitness ~ Inv.sp * coev + (1|Res.com.ID/Block), data = la.dat.exp)

#model checks
par(mfrow = c(2,2))
plot(la.max)
plot(la.max, resid(., scaled=TRUE) ~ fitted(.) | Block, abline = 0)
plot(la.max, resid(., scaled=TRUE) ~ fitted(.) | Res.com.ID, abline = 0)

#looks homoscedastic
par(mfrow = c(1,1))
qqnorm(resid(la.max))
#some deviations but overall looks normal enough

#Does invasion sp ID and coev interact?
la.max.noint <- lmer(rel.fitness ~ Inv.sp + coev + (1|Block), data = la.dat.exp)
anova(la.max, la.max.noint) #no sig effect of including interaction (P = 0.968)

summary(la.max.noint)

#Effect of coev
la.max.nc <- update(la.max.noint,~. - coev)
anova(la.max.nc, la.max.noint)
#Coevo non-sig. P = 0.3395

#Effect of sp ID
la.nsp <- update(la.max.noint,~. - Inv.sp)
anova(la.max.noint, la.nsp)
#Sp ID sig. Drop coev
summary(la.max.noint)

la.null <- update(la.max.nc,~. - Inv.sp)
anova(la.null, la.max.nc)
#Significant effect of species ID on relative invader fitness P<0.001

#Post-hoc analysis - Tukey multiple comparisons
emmeans::emmeans(la.max.nc, pairwise ~ Inv.sp)

#Significant pairwise comparisons between AA and other species but not PC. Negative estimates- AA lower rel fitness
#OD sig against AA, PC and SR. Lower rel fitness than SR but greater than PC and AA. Non-sig against VG
#PC sig against all except AA. Lower rel fitness than SR, VG and OD.
#SR sig greater rel fitness than all other species. 
#VG- sig greater against AA, PC but sig lower than SR. Nons-sig to OD

#Plot showing lack of coev effect (ignoring species)

ggplot(la.dat.expx, aes(x=coev, y=rel.fitness), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  scale_y_continuous(limits = c(0.7, 2)) +
  ylab(label = "Relative invader fitness") +
  xlab(label = "Invader coevolutionary history to community") +
  scale_x_discrete(labels = c("Non-coev.", "Coev."))

#Separate out species and coev- something interesting perhaps with variation in invasion success between coev and non-coev species

ggplot(la.dat.expx, aes(x=coev, y=rel.fitness), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Inv.sp) +
  scale_y_continuous(limits = c(0.7, 2)) +
  ylab(label = "Relative invader fitness") +
  xlab(label = "Invader coevolutionary history to community") +
  scale_x_discrete(labels = c("Non-coev.", "Coev."))

#Plot of this relative species invasion sucess ignoring coevolution:

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


#for some reason it no longer likes the wrangled dataframe

la.dat.expx$Inv.sp <- as.character(la.dat.expx$Inv.sp)
la.dat.expx$coev <- as.character(la.dat.expx$coev)

ggplot(la.dat.expx, aes(x=forcats::fct_reorder(Inv.sp, coev), y=rel.fitness), na.rm=T) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_point(position = position_jitter(0.2), aes(shape = coev), size =1.5) +
  scale_shape_manual(values=c(15,16), name = "Invader type", labels = c("Non-coev.", "Coev.")) +
  scale_y_continuous(limits = c(0.7, 2)) +
  scale_x_discrete(labels = labels) +
  xlab(label = "Invading species") +
  ylab(label = "Relative invader fitness")


#naive GLM reaches same conclusion even without block included. Formal analyses will include block effect anyway to account for experimental design
lol <- glm(rel.fitness ~ Inv.sp * coev, data = la.dat.exp)
lol2 <- glm(rel.fitness ~ Inv.sp + coev, data = la.dat.exp)
anova(lol, lol2, test = "F") #still non-sig
lol3 <- glm(rel.fitness ~ Inv.sp, data = la.dat.exp)
anova(lol2, lol3, test = "F") #still non-sig
lol4 <- glm(rel.fitness ~ Inv.sp, data = la.dat.exp) 
anova(lol2, lol4, test = "F") #still non-sig
lol5 <- glm(rel.fitness ~ 1, data = la.dat.exp)
anova(lol4, lol5, test = "F") #still sig

## Effect of invader on community productivity? 
exp1.prod <- read.csv("Exp1_inv_prod.csv", header = T)

prod.mod1 <- lmer(Comm.prod ~ Inv.sp * coev + (1|Res.com.ID/Block), exp1.prod)
prod.mod2 <- lmer(Comm.prod ~ Inv.sp + coev + (1|Res.com.ID/Block), exp1.prod)
anova(prod.mod1, prod.mod2)
#no interaction. P = 0.959

prod.mod3 <- lmer(Comm.prod ~ Inv.sp + (1|Res.com.ID/Block), exp1.prod)
anova(prod.mod2, prod.mod3)
#no effect of coevolution. P = 0.761

prod.mod4 <- lmer(Comm.prod ~ coev + (1|Res.com.ID/Block), exp1.prod)
anova(prod.mod2, prod.mod4)
#no effect of invasive species ID. P = 0.697. Drop ID as coev is less sig

prod.mod5 <- lmer(Comm.prod ~ 1 + (1|Res.com.ID/Block), exp1.prod)
anova(prod.mod3, prod.mod5) #no sig effect of coev. 
#Conclusion: invasive species ID and coev does not effect total comm productivity

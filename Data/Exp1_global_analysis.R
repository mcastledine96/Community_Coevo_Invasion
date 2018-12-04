#Try global model

glob <- read.csv("Exp1_global2.csv", header = T)

glob <- mutate(glob,
               Res.Fit = log(T1.res/T0.res),
               Inv.Fit = log(T1.inv/T0.inv))

ggplot(glob, aes(x = Inv.sp, y = Res.Fit), na.rm = T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2), aes(shape = coev)) +
  facet_wrap(~Res.ID, scales="free")

#Try linear mixed effect model with multi-model inference

glob2 <- filter(glob, ! Treatment == 15) #remove problem datapoints
glob2 <- filter(glob2, ! Treatment == 40)
glob2 <- filter(glob2, ! Treatment == 82)

mod <- lmer(Res.Fit ~ Inv.sp*Inv.Fit + Res.ID + coev + (1|Res.com.ID/Block), data = glob2, na.action = "na.fail")

library(MuMIn)
mod.dred <- dredge(mod)

mod.dred #Simplest model only contains Resident species ID

#model with two factors interacting to create broad invader id
glob2$int <- interaction(glob2$coev, glob2$Inv.sp)

ggplot(glob2, aes(x = int, y = Res.Fit), na.rm = T) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(0.2)) +
  facet_wrap(~Res.ID, scales="free")

mod3 <- lmer(Res.Fit ~ int + Inv.Fit + Res.ID + (1|Block), data = glob2, na.action = "na.fail")

mod3.dredge <- dredge(mod3)
mod3.dredge

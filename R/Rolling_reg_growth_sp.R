library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(ggridges)
library(stringr)
install.packages("zoo")

files <- list.files("C:/Users/Meaghan/Documents/4th Year/Research project/Data/MSci_analysesR/Growth_curves", pattern = ".txt", full.names = T)

filer <- "C:/Users/Meaghan/Documents/4th Year/Research project/Data/MSci_analysesR/Growth_curves/SR_ct_mb.txt"

files <- files[! files %in% filer]

d.growth <- purrr::map_df(files, bind_biolog_sheet, .id = "sp")

file.data <- data.frame(files = files, stringsAsFactors = F) %>%
  mutate(., species = substr(basename(files), 1, 2),
         sp = 1:5) %>%
  select(- files)

d.growth <- merge(d.growth, file.data, by = "sp") %>%
  separate(., well, c('row', 'column'), sep = '_', remove = FALSE)

# specific cases to bind into dataframe
d.special <- filter(d.growth, well == "A_8" & species == "SR")


d.growth2 <- filter(d.growth, row%in%c("B", "G"))
d.growth2[which(d.growth2$well == "B_10" & d.growth2$species == "SR"), ]$OD <- d.special$OD



d.growth2 <- mutate(d.growth2, 
                    t = (id*20)/60,
                    log10_od = log10(OD))


ggplot(d.growth2, aes(x= t, y= OD)) +
  geom_point(aes(col = row)) +
  facet_wrap(~species) 

# run a rolling regression ####

# make a function for output
output <- function(x){
  fm <- lm(as.data.frame(x))
  slope <- coef(fm)[[2]]
  ci <- confint(fm)[2, ]
  c(slope = slope, conf_lower = ci[[1]], conf_upper = ci[[2]])
}
rolling_coefs <- function(x, window_size) zoo::rollapplyr(x, window_size, output, by.column = FALSE, fill = NA)

models <- d.growth2 %>%
  group_by(species, well, row, column) %>%
  arrange(., t) %>%
  do(cbind(., select(., log10_od, t) %>% rolling_coefs(window_size = 5)))

d_growth_rate <- mutate(models, slope = ifelse(is.infinite(slope) | is.nan(slope), NA, slope)) %>%
  group_by(species, well, row, column) %>%
  filter(., slope == max(slope, na.rm = TRUE)) %>%
  ungroup()

min(d_growth_rate$slope)

ggplot(d_growth_rate, aes(species, slope, col = row)) +
  geom_point(position = position_dodge(0.5)) 

##

d_growth_rate <- janitor::clean_names(d_growth_rate) %>%
  mutate(.,
         K = temp + 273.15,
         curve_id = group_indices(., growth_temp, rep, phage),
         growth_temp2 = paste('gt_', growth_temp, sep = ''))

# model using standard approach in nls.multstart ####
sharpeschoolhigh_1981 <- function(temp_k, rtref, e, eh, th, tref){
  tref <- 273.15 + tref
  k <- 8.62e-05
  boltzmann.term <- rtref*exp(e/k * (1/tref - 1/temp_k))
  inactivation.term <- 1/(1 + exp(eh/k * (1/th - 1/temp_k)))
  return(boltzmann.term * inactivation.term)
}


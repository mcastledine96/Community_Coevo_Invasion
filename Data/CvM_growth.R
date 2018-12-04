library(dplyr)
library(ggplot2)

bind_biolog_sheet <- function(file){
  
  # read in file
  dat <- readr::read_tsv(file, col_names = F, skip = 3)
  
  # process file
  d <- dplyr::select(dat, 2:13) %>%
    tidyr::drop_na(.) %>%
    dplyr::mutate(., sum_row = rowSums(.)) %>%
    dplyr::filter(., sum_row != sum(1:12)) %>%
    dplyr::select(., -sum_row)
  
  # work out how many plates per sheet
  num_plates <- nrow(d)/8
  
  d <- dplyr::mutate(d, id = rep(seq(1, num_plates, 1), each = 8),
                     well = rep(LETTERS[1:8],num_plates)) %>%
    tidyr::gather(., 'col_plate', 'OD', starts_with('X')) %>%
    mutate(well = paste(well, readr::parse_number(col_plate) - 1, sep = '_')) %>%
    dplyr::select(., -col_plate)
  return(d)
  
}


##SR 2

file <- "C:/Users/Meaghan/Documents/4th Year/Research project/Data/MSci_analysesR/Growth_curves/SR_ct_mb2.txt"

x <- bind_biolog_sheet(file)

dat <- filter(x, well == "B_1" | well == "B_2" | well == "B_3" | well == "B_4" | well == "B_5" | well == "B_6" | well == "B_7"| well == "B_8" | well == "B_9" | well == "B_10" | well == "B_11" | well == "B_12" | well == "G_1" | well == "G_2" | well == "G_3" | well == "G_4" | well == "G_5" | well == "G_6" | well == "G_7" | well == "G_8" | well == "G_9" | well == "G_10" | well == "G_11" | well == "G_12")

new_well <- filter(x, well == "A_8")
tr2 <- rep(1, length.out = 49)
tr2 <- as.factor(tr2)
new_well <- mutate(new_well,
                   Tr = tr2)

tr <- seq(1, 0)
tr <- rep(tr, length.out = 1176)
tr <- as.factor(tr)

dat <- mutate(dat,
              Tr = tr)

dat <- merge(dat, new_well, all.x = T, all.y = T)


dat <- filter(dat, well == "A_8" | well == "B_1" | well == "B_2" | well == "B_3" | well == "B_4" | well == "B_5" | well == "B_6" | well == "B_7"| well == "B_8" | well == "B_9" | well == "B_11" | well == "B_12" | well == "G_1" | well == "G_2" | well == "G_3" | well == "G_4" | well == "G_5" | well == "G_6" | well == "G_7" | well == "G_8" | well == "G_9" | well == "G_10" | well == "G_11" | well == "G_12")

ggplot(dat, aes(x = id, y = OD)) +
  geom_point(aes(col = Tr))

##VG growth

file <- "C:/Users/Meaghan/Documents/4th Year/Research project/Data/MSci_analysesR/Growth_curves/VG_ct_mb.txt"

x <- bind_biolog_sheet(file)

dat <- filter(x, well == "B_1" | well == "B_2" | well == "B_3" | well == "B_4" | well == "B_5" | well == "B_6" | well == "B_7"| well == "B_8" | well == "B_9" | well == "B_10" | well == "B_11" | well == "B_12" | well == "G_1" | well == "G_2" | well == "G_3" | well == "G_4" | well == "G_5" | well == "G_6" | well == "G_7" | well == "G_8" | well == "G_9" | well == "G_10" | well == "G_11" | well == "G_12")

tr <- seq(1, 0)
tr <- rep(tr, length.out = 1152)
tr <- as.factor(tr)

dat <- mutate(dat,
              Tr = tr)

dat <- filter(dat, well == "B_1" | well == "B_2" | well == "B_3" | well == "B_4" | well == "B_5" | well == "B_6" | well == "B_9" | well == "B_10" | well == "B_11" | well == "B_12" | well == "G_1" | well == "G_2" | well == "G_3" | well == "G_4" | well == "G_5" | well == "G_6" | well == "G_7" | well == "G_8" | well == "G_9" | well == "G_10" | well == "G_11" | well == "G_12")

ggplot(dat, aes(x = id, y = OD)) +
  geom_point(aes(col = Tr)) +
  scale_y_continuous(limits = c(0.05, 0.25)) +
  scale_x_continuous(limits = c(20, 50))

##OD growth 

file <- "C:/Users/Meaghan/Documents/4th Year/Research project/Data/MSci_analysesR/Growth_curves/OD_ct_mb.txt"

x <- bind_biolog_sheet(file)

dat <- filter(x, well == "B_1" | well == "B_2" | well == "B_3" | well == "B_4" | well == "B_5" | well == "B_6" | well == "B_7"| well == "B_8" | well == "B_9" | well == "B_10" | well == "B_11" | well == "B_12" | well == "G_1" | well == "G_2" | well == "G_3" | well == "G_4" | well == "G_5" | well == "G_6" | well == "G_7" | well == "G_8" | well == "G_9" | well == "G_10" | well == "G_11" | well == "G_12")

tr <- seq(1, 0)
tr <- rep(tr, length.out = 1176)
tr <- as.factor(tr)

dat <- mutate(dat,
              Tr = tr)

dat <- filter(dat, well == "B_1" | well == "B_2" | well == "B_3" | well == "B_4" | well == "B_5" | well == "B_6" | well == "B_9" | well == "B_10" | well == "B_11" | well == "B_12" | well == "G_1" | well == "G_2" | well == "G_3" | well == "G_4" | well == "G_5" | well == "G_6" | well == "G_7" | well == "G_8" | well == "G_9" | well == "G_10" | well == "G_11" | well == "G_12")

ggplot(dat, aes(x = id, y = OD)) +
  geom_point(aes(col = Tr))
  

##AA growth

file <- "C:/Users/Meaghan/Documents/4th Year/Research project/Data/MSci_analysesR/Growth_curves/AA_ct_mb.txt"

x <- bind_biolog_sheet(file)

dat <- filter(x, well == "B_1" | well == "B_2" | well == "B_3" | well == "B_4" | well == "B_5" | well == "B_6" | well == "B_7"| well == "B_8" | well == "B_9" | well == "B_10" | well == "B_11" | well == "B_12" | well == "G_1" | well == "G_2" | well == "G_3" | well == "G_4" | well == "G_5" | well == "G_6" | well == "G_7" | well == "G_8" | well == "G_9" | well == "G_10" | well == "G_11" | well == "G_12")

tr <- seq(1, 0)
tr <- rep(tr, length.out = 1176)
tr <- as.factor(tr)

dat <- mutate(dat,
              Tr = tr)

ggplot(dat, aes(x = id, y = OD)) +
  geom_point(aes(col = Tr))


##PC growth

file <- "C:/Users/Meaghan/Documents/4th Year/Research project/Data/MSci_analysesR/Growth_curves/PC_ct_mb.txt"

x <- bind_biolog_sheet(file)

dat <- filter(x, well == "B_1" | well == "B_2" | well == "B_3" | well == "B_4" | well == "B_5" | well == "B_6" | well == "B_7"| well == "B_8" | well == "B_9" | well == "B_10" | well == "B_11" | well == "B_12" | well == "G_1" | well == "G_2" | well == "G_3" | well == "G_4" | well == "G_5" | well == "G_6" | well == "G_7" | well == "G_8" | well == "G_9" | well == "G_10" | well == "G_11" | well == "G_12")

tr <- seq(1, 0)
tr <- rep(tr, length.out = 1176)
tr <- as.factor(tr)

dat <- mutate(dat,
              Tr = tr)

ggplot(dat, aes(x = id, y = OD)) +
  geom_point(aes(col = Tr)) #looks like mono lines have a greater carrying capacity

#just community lines
com.pc <- filter(dat, Tr == "1")

ggplot(com.pc, aes(x = id, y = OD)) +
  geom_point()

#just mono lines
mono.pc <- filter(dat, Tr == "0")

ggplot(mono.pc, aes(x = id, y = OD)) +
  geom_point()

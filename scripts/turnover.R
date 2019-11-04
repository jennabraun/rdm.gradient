#analytical steps
#aggregrate samples to shrub/open by sites
#transform into relative abundances because we lost some samples
#figure out what the spatial gradient is
#calculate spatial turnover

library(tidyverse)
env <- read.csv("clean_data/cov.csv")
comm <- read.csv("clean_data/comm.csv")
#remove larrea
env <- filter(env, Microsite != 'larrea')
#subset non-matching rows
comm <- comm[comm$X %in% env$X == TRUE,]
comm <- separate(comm, X, c('site', 'microsite', 'rep'), sep = " ")

comm <- comm %>% select(-rep) 




%>% group_by(site,microsite) %>% summarise_at(vars(c(3)), 


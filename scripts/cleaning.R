#data cleaning R script
library(tidyverse)
library(XLConnect)
rdm <- readWorksheetFromFile("rdm.gradient.2019.datasheet.xlsx", "data")
str(rdm)
rdm$Microsite <- as.factor(rdm$Microsite)
rdm$max.veg.height <- as.numeric(rdm$max.veg.height)

ggplot(rdm, aes(Microsite, RDM)) + geom_boxplot() + facet_wrap(~Region)
ggplot(rdm, aes(Microsite, rdm.cov)) + geom_boxplot() + facet_wrap(~Region)
ggplot(rdm, aes(Microsite, max.veg.height)) + geom_boxplot() + facet_wrap(~Region)


arth <- readWorksheetFromFile("arthropods.xlsx", "Main ID")

str(arth)
arth$Microsite <- as.factor(arth$Microsite)
arth$Microsite <- gsub(" ", "", arth$Microsite)

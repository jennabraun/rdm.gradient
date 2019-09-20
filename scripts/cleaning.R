#data cleaning R script
library(tidyverse)
library(XLConnect)
library(vegan)
rdm <- readWorksheetFromFile("rdm.gradient.2019.datasheet.xlsx", "data")
str(rdm)
rdm$Microsite <- as.factor(rdm$Microsite)
rdm$max.veg.height <- as.numeric(rdm$max.veg.height)

ggplot(rdm, aes(Microsite, RDM)) + geom_boxplot() + facet_wrap(~Region)
ggplot(rdm, aes(Microsite, rdm.cov)) + geom_boxplot() + facet_wrap(~Region)
ggplot(rdm, aes(Microsite, max.veg.height)) + geom_boxplot() + facet_wrap(~Region)

#arthopod data cleaning

arth <- readWorksheetFromFile("arthropods.xlsx", "Main ID")
str(arth)
arth$Microsite <- as.factor(arth$Microsite)
arth$Microsite <- gsub(" ", "", arth$Microsite)
arth$Microsite <- gsub("shrub", "ephedra", arth$Microsite)
arth$uniID <- paste(arth$Site, arth$Microsite, arth$Rep)
arth$morph <- paste(arth$highest.rtu, arth$Morpho_code)
#aggregrate observations
arth.ag <- dplyr::select(arth, uniID, morph, Quantity)
arth.ag <- arth.ag %>% group_by(uniID, morph) %>% summarise(Quantity = sum(Quantity)) 
sum(arth.ag$Quantity)

arth.ag$morph <- gsub(" ","", arth.ag$morph)

wide <- arth.ag %>% spread(morph, Quantity)

wide <- select(wide, -"ignoreNA", -"damagedNA", -"NA1", -"NA2", -"NA3", -"NANA")
wide[is.na(wide)] <- 0


##make metadata dataframe
metadata <- rdm
metadata$uniID <- paste(metadata$site, metadata$Microsite, metadata$ID)


row.names(metadata) <- metadata$uniID
row.names(wide) <- wide$uniID
metadata <- metadata[match(wide$uniID, metadata$uniID),]
missing <- anti_join(wide, metadata, by = "uniID")
wide <- wide %>% ungroup(uniID) %>% select(-uniID)

all.equal(rownames(wide), rownames(metadata))


metadata$abun <- apply(wide, 1, sum)
#check for total
sum(metadata$abun)
H <- diversity(wide)
simp <- diversity(wide, "simpson")
S <- specnumber(wide)
J <- H/log(S)
metadata$H <- H
metadata$Simpson <- simp
metadata$Species <- S
metadata$Even <- J



ggplot(metadata, aes(Microsite, Species)) + geom_boxplot() + facet_wrap(~Region)
ggplot(metadata, aes(Microsite, Simpson)) + geom_boxplot() + facet_wrap(~Region)

ggplot(metadata, aes(Microsite, abun)) + geom_boxplot() + facet_wrap(~Region)

ggplot(metadata, aes(RDM, Species, color = Microsite)) + geom_smooth() + facet_wrap(~Region)


ggplot(rdm, aes(RDM)) + geom_density()
qqnorm(rdm$RDM)
sq <- rdm$RDM^2
shapiro.test(sq)

plot(sq)

m1 <- lm(RDM~Microsite + Region, data = rdm)
shapiro.test(residuals(m1))
summary(m1)

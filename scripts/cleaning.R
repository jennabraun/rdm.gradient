#data cleaning R script
library(tidyverse)
library(XLConnect)
library(vegan)
rdm <- readWorksheetFromFile("raw_data/rdm.gradient.2019.datasheet.xlsx", "data")
str(rdm)
rdm$Microsite <- as.factor(rdm$Microsite)
rdm$max.veg.height <- as.numeric(rdm$max.veg.height)
write.csv(rdm, "clean_data/veg_covariates.csv")

#arthopod data cleaning
arth <- readWorksheetFromFile("raw_data/arthropods.xlsx", "Main ID")
str(arth)
arth$Microsite <- as.factor(arth$Microsite)
arth$Microsite <- gsub(" ", "", arth$Microsite)
arth$Microsite <- gsub("shrub", "ephedra", arth$Microsite)
arth$uniID <- paste(arth$Site, arth$Microsite, arth$Rep)
arth$morph <- paste(arth$highest.rtu, arth$Morpho_code)


#filter out ignores
arth <- filter(arth, highest.rtu != "alate" & highest.rtu != "ignore" & highest.rtu != "damaged")

sum(arth$Quantity)

write.csv(arth, "clean_data/arth_long.csv")
arth <- read.csv("clean_data/arth_long.csv")
#aggregrate observations
arth <- filter(arth, Order != "Lepidoptera")
arth.ag <- dplyr::select(arth, uniID, morph, Quantity)
arth.ag <- arth.ag %>% group_by(uniID, morph) %>% summarise(Quantity = sum(Quantity)) 
sum(arth.ag$Quantity)

arth.ag$morph <- gsub(" ","", arth.ag$morph)


wide <- arth.ag %>% spread(morph, Quantity)
wide[is.na(wide)] <- 0
row.names(wide) <- wide$uniID
##make metadata dataframe
metadata <- rdm
metadata$uniID <- paste(metadata$site, metadata$Microsite, metadata$ID)
row.names(metadata) <- metadata$uniID

#row.names(metadata) <- metadata$uniID

metadata <- metadata[match(wide$uniID, metadata$uniID),]

missing <- anti_join(wide, metadata, by = "uniID")
all.equal(rownames(wide), rownames(metadata))
wide <- wide %>% ungroup(uniID) %>% select(-uniID)


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

#checking to ensure rep names are correct
rdm$uniID <- paste(rdm$site, rdm$Microsite, rdm$ID)
test <- left_join(rdm, metadata, by = "uniID")
test <- select(test, PF.x, abun, uniID)


t <- filter(test, PF.x == "Y" & is.na(abun) == TRUE)
t <- filter(test, is.na(PF.x) == TRUE & is.na(abun) == FALSE)

count(test, PF.x == "N")

write.csv(metadata, "clean_data/cov.csv")
write.csv(wide, "clean_data/comm.csv")


#I want the same community and env without ephedra
arth.ag <- dplyr::select(arth, uniID, morph, Microsite, Quantity)
arth.ag <- filter(arth.ag, Microsite != "larrea")
#arth.ag <- filter(arth.ag, uniID != "PAN2 open 11")
arth.ag <- arth.ag %>% group_by(uniID, morph) %>% summarise(Quantity = sum(Quantity)) 
sum(arth.ag$Quantity)

arth.ag$morph <- gsub(" ","", arth.ag$morph)
wide.eph <- arth.ag %>% spread(morph, Quantity)
wide.eph[is.na(wide.eph)] <- 0
row.names(metadata) <- metadata$uniID
row.names(wide.eph) <- wide.eph$uniID
metadata <- metadata[match(wide.eph$uniID, metadata$uniID),]
missing <- anti_join(wide.eph, metadata, by = "uniID")
wide.eph <- wide.eph %>% ungroup(uniID) %>% select(-uniID)

all.equal(rownames(wide.eph), rownames(metadata))
metadata$abun <- apply(wide.eph, 1, sum)
#check for total
sum(metadata$abun)
H <- diversity(wide.eph)
simp <- diversity(wide.eph, "simpson")
S <- specnumber(wide.eph)
J <- H/log(S)
metadata$H <- H
metadata$Simpson <- simp
metadata$Species <- S
metadata$Even <- J



write.csv(metadata, "clean_data/cov_eph.csv")
write.csv(wide.eph, "clean_data/comm_eph.csv")
write.csv(wide, "clean_data/comm.csv")
#make a data frame of the regions each morphospecies is found
arth$pres <- "Y"
arth <- select(arth, morph, Region, pres)
arth <- distinct(arth)
regions <- arth %>% pivot_wider(., morph, names_from = Region, values_from = pres)
write.csv(regions, "clean_data/regions.csv")
morphs <- read.csv("clean_data/morphospecies_descriptions.csv")


morphs$morph <- paste(morphs$highest.taxanomic.resolution, morphs$Morphotype)
morphs$morph <- gsub(" ", "", morphs$morph)
regions$morph <- gsub(" ", "", regions$morph)
morphs <- left_join(morphs, regions, by = "morph")
a <- filter(morphs, is.na(Cuyama) == TRUE & is.na(Mojave) == TRUE & is.na(Panoche == "TRUE"))
write.csv(morphs, "clean_data/morphospecies_descriptions.csv")

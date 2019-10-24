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


#aggregrate observations
arth.ag <- dplyr::select(arth, uniID, morph, Quantity)
arth.ag <- arth.ag %>% group_by(uniID, morph) %>% summarise(Quantity = sum(Quantity)) 
sum(arth.ag$Quantity)

arth.ag$morph <- gsub(" ","", arth.ag$morph)

wide <- arth.ag %>% spread(morph, Quantity)
wide <- select(wide, -"ignoreNA", -"damagedNA")
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

#checking to ensure rep names are correct
rdm$uniID <- paste(rdm$site, rdm$Microsite, rdm$ID)
test <- left_join(rdm, metadata, by = "uniID")
test <- select(test, PF.x, abun, uniID)


t <- filter(test, PF.x == "Y" & is.na(abun) == TRUE)
t <- filter(test, is.na(PF.x) == TRUE & is.na(abun) == FALSE)

count(test, PF.x == "N")

write.csv(metadata, "clean_data/cov.csv")
write.csv(wide, "clean_data/comm.csv")

#potential prey items only
codes <- read.csv("clean_data/morphospecies_descriptions.csv")
codes$morpho <- paste(codes$highest.rtu, codes$Morphotype)
codes$morpho <- gsub(" ", "", codes$morpho)
wide$uniID <- row.names(wide)
prey <- gather(wide, morpho, count, 1:207)
prey <- right_join(prey, codes, by = "morpho")
prey <- filter(prey, BNLL.prey.item. == "Y")
prey <- select(prey, 1:3)
prey <- spread(prey, morpho, count)
prey <- filter(prey, uniID != "NA")
row.names(prey) <- prey$uniID
prey <- select(prey, -uniID)
all.equal(row.names(prey), row.names(metadata))
str(prey)

#get rid of columns containing NAs becasue they are unassigned morphospecies
prey <- prey %>% select_if(~ !any(is.na(.)))

metadata$prey.abun <- apply(prey, 1, sum)
metadata$prey.Species <- specnumber(wide)

write.csv(prey, "clean_data/bnll_prey.csv")
write.csv(metadata, "clean_data/cov.csv")

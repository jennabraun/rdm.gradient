#data cleaning R script


#make a data frame of the regions each morphospecies is found
# arth$pres <- "Y"
library(dplyr)
library(tidyr)
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


#filter out ignores
arth <- filter(arth, highest.rtu != "alate" & highest.rtu != "ignore" & highest.rtu != "damaged" & highest.rtu != "missing" & highest.rtu != "immature")
arth <- filter(arth, Order != "Lepidoptera")
sum(arth$Quantity)

write.csv(arth, "clean_data/arth_long.csv")
arth <- read.csv("clean_data/arth_long.csv")
#aggregrate observations





arth.ag <- dplyr::select(arth, uniID, highest.rtu, Quantity)
arth.ag <- arth.ag %>% group_by(uniID, highest.rtu) %>% summarise(Quantity = sum(Quantity)) 
sum(arth.ag$Quantity)

arth.ag$highest.rtu <- gsub(" ","", arth.ag$highest.rtu)

chk <- arth.ag %>% ungroup() %>% select(highest.rtu) %>% count(highest.rtu)

wide <- arth.ag %>% spread(highest.rtu, Quantity)
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


#with no singletons 
total <- arth.ag %>% group_by(highest.rtu) %>% summarise(sum = sum(Quantity))
unique(total$highest.rtu)
singles <- total %>% filter(sum == 1)
singles <- singles$highest.rtu
'%notin%' <- Negate('%in%')
nosingles <- arth.ag %>% filter(!highest.rtu %in% singles)
sum(nosingles$Quantity)

wide.nosingles <- nosingles%>% spread(highest.rtu, Quantity)
wide.nosingles[is.na(wide.nosingles)] <- 0
row.names(wide.nosingles) <- wide.nosingles$uniID




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

#add singletonless calculations to metadata

metadata <- metadata[match(wide.nosingles$uniID, metadata$uniID),]

missing <- anti_join(wide.nosingles, metadata, by = "uniID")
all.equal(rownames(wide.nosingles), rownames(metadata))
wide.nosingles <- wide.nosingles %>% ungroup(uniID) %>% select(-uniID)


metadata$abun.nos <- apply(wide.nosingles, 1, sum)
#check for total
sum(metadata$abun.nos)

S <- specnumber(wide.nosingles)

metadata$Species.nos <- S


#checking to ensure rep names are correct
rdm$uniID <- paste(rdm$site, rdm$Microsite, rdm$ID)
test <- left_join(rdm, metadata, by = "uniID")
test <- select(test, PF.x, abun, uniID)


t <- filter(test, PF.x == "Y" & is.na(abun) == TRUE)
t <- filter(test, is.na(PF.x) == TRUE & is.na(abun) == FALSE)

count(test, PF.x == "N")

write.csv(metadata, "clean_data/cov.csv")
write.csv(wide, "clean_data/comm.csv")
write.csv(wide.nosingles, "clean_data/comm_nos.csv")




# arth <- select(arth, morph, Region, pres)
# arth <- distinct(arth)
# regions <- arth %>% pivot_wider(., morph, names_from = Region, values_from = pres)
# write.csv(regions, "clean_data/regions.csv")
# morphs <- read.csv("clean_data/morphospecies_descriptions.csv")
# 
# 
# morphs$morph <- paste(morphs$highest.taxanomic.resolution, morphs$Morphotype)
# morphs$morph <- gsub(" ", "", morphs$morph)
# regions$morph <- gsub(" ", "", regions$morph)
# morphs <- left_join(morphs, regions, by = "morph")
# a <- filter(morphs, is.na(Cuyama) == TRUE & is.na(Mojave) == TRUE & is.na(Panoche == "TRUE"))
# write.csv(morphs, "clean_data/morphospecies_descriptions.csv")

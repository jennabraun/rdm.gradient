old bits



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
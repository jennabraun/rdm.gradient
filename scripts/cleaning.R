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
filter(arth, Family == "Formicidae") %>% write.csv(.,"desertants.csv")

str(arth)
arth$Microsite <- as.factor(arth$Microsite)
arth$Microsite <- gsub(" ", "", arth$Microsite)
arth$Microsite <- gsub("shrub", "ephedra", arth$Microsite)
arth$uniID <- paste(arth$Site, arth$Microsite, arth$Rep)
arth$morph <- paste(arth$highest.rtu, arth$Morpho_code)

arth <- filter(arth, highest.rtu != "ignore" & highest.rtu != "damaged")



#aggregrate observations
arth.ag <- dplyr::select(arth, uniID, morph, Quantity)
arth.ag <- arth.ag %>% group_by(uniID, morph) %>% summarise(Quantity = sum(Quantity)) 
sum(arth.ag$Quantity)

arth.ag$morph <- gsub(" ","", arth.ag$morph)

wide <- arth.ag %>% spread(morph, Quantity)

wide <- select(wide, -"ignoreNA", -"damagedNA", -"NANA")
wide[is.na(wide)] <- 0

#remove destroyed reps


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


filter()

#checking to ensure rep names are correct
rdm$uniID <- paste(rdm$site, rdm$Microsite, rdm$ID)
test <- left_join(rdm, metadata, by = "uniID")
test <- select(test, PF.x, abun, uniID)


t <- filter(test, PF.x == "Y" & is.na(abun) == TRUE)
t <- filter(test, is.na(PF.x) == TRUE & is.na(abun) == FALSE)

count(test, PF.x == "N")


#remove missing reps.. maybe we don't need to because they are na and not 0. 

#EDA for data checking and general variable distributions, correlations etc

#density plots
write.csv(metadata, "cov.csv")
write.csv(wide, "comm.csv")
ggplot(rdm, aes(rdm)) + geom_density()
ggplot(rdm, aes(max.veg.height)) + geom_density()
ggplot(rdm, aes(bur.drip)) + geom_histogram()
ggplot(rdm, aes(rdm.cov)) + geom_density()
ggplot(rdm, aes(green.cov)) + geom_density()


ggplot(metadata, aes(abun)) + geom_histogram()

#without the two really big samples
metadata %>% filter(abun < 400) %>% ggplot(aes(abun)) + geom_histogram()

metadata %>% filter(abun < 150) %>% ggplot(aes(abun)) + geom_histogram()


ggplot(metadata, aes(Species)) + geom_histogram()

ggplot(metadata, aes(H)) + geom_density()
shapiro.test(metadata$H)
#holy moly p = 0.06155 its  normal


ggplot(metadata, aes(Even)) + geom_density()
shapiro.test(metadata$Even)


m1 <- glm(Species ~ RDM  + Microsite + site, family = "poisson", data = eph)
summary(m1)
car::Anova(m1, type = 3)


library(lsmeans)
lsmeans(m1, pairwise~Microsite, adjust = "Tukey")
lsmeans(m1, pairwise~RDM|Region, adjust = "Tukey")

ggplot(metadata, aes(RDM, Species)) + geom_smooth(method = lm) + facet_wrap(~Region)


m2 <- lm(H ~ Microsite, data = metadata)
m3 <- glmmTMB::glmmTMB(Species ~ Microsite + site+(1|Region), family = gaussian, data = eph)
summary(m3)


eph <- filter(metadata, Microsite != "larrea")

shapiro.test(residuals(m2))
summary(m2)
car::Anova(m2, type = 2)

m3 <- lm(Even ~ RDM + Microsite + site, data = eph)
shapiro.test(residuals(m3))
summary(m3)
car::Anova(m3, type = 2)

m3 <- glm(abun ~ RDM + Region * Microsite, family = "poisson", data = metadata)
summary(m3)

#remove that huge sample

m4 <- metadata %>% filter(abun < 350) %>% glm(abun ~ RDM + site + Microsite, family = "quasipoisson", data = .)
#metadata$abun.scl <- scale(metadata$abun)
summary(m4)



m5 <- glm(abun.scl ~ RDM + site + Microsite, family = "quasipoisson", data = metadata)
summary(m5)

metadata %>% filter(abun < 350) %>% ggplot(aes(abun)) + geom_histogram(binwidth = 1)

metadata %>% group_by(Region,Microsite) %>% summarise(mean(Simpson))

#not sure how to deal with dispersion issues at the moment, but it appears that outlier is really driving the negative influence of RDM

summary(m4)

plot(residuals(m4))

  
ggplot(metadata, aes(Microsite, Species)) + geom_boxplot() + facet_wrap(~Region)
ggplot(metadata, aes(Microsite, Simpson)) + geom_boxplot() + facet_wrap(~Region)
ggplot(metadata, aes(Microsite, abun)) + geom_boxplot() + facet_wrap(~Region)

ggplot(metadata, aes(x, RDM, color = Microsite)) + geom_smooth() + facet_wrap(~Region)


ggplot(rdm, aes(RDM)) + geom_density()
qqnorm(rdm$RDM)
sq <- rdm$RDM^2
shapiro.test(sq)

plot(sq)

m1 <- lm(Species~Microsite * Region + RDM, data = metadata)
shapiro.test(residuals(m1))
summary(m1)

eph <- 

env <- metadata %>% select(Region, Microsite, RDM, max.veg.height)
str(env)
env$Region <- as.factor(env$Region)
wide <- decostand(wide, method = "hellinger")
#quick ordination
c1 <- cca(wide ~ Region + Microsite + RDM, data = env)
c1
summary(c1)
car::Anova(c1)
RsquareAdj(c1)

#try a perMANOVA
p1 <- adonis(wide ~ Region * Microsite + RDM, env)
summary(p1)
p1


library(indicspecies)
indval <- multipatt(wide, env$Region)
indval
summary(indval)


pan <- 

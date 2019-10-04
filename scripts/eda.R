library(tidyverse)
library(purrr)
library(ggthemes)

env <- read.csv("cov.csv")
comm <- read.csv("comm.csv")
rdm <- read.csv("rdm.csv")

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


ggplot(rdm, aes(RDM)) + geom_density() + facet_grid(~site)
ggplot(env, aes(Species, fill = Microsite)) + geom_density(alpha = 0.7) + facet_grid(~site)
env %>% filter(abun < 350) %>% ggplot(aes(abun, fill = Microsite)) + geom_density(alpha = 0.7) + facet_grid(~site)
ggplot(rdm, aes(x, RDM)) + geom_smooth() + facet_grid(~site) + xlab("Shrub Width")
ggplot(rdm, aes(y, RDM)) + geom_smooth() + facet_grid(~site) + xlab("Shrub Height")
ggplot(rdm, aes(Microsite, RDM)) + geom_boxplot() + facet_grid(~site)
ggplot(env, aes(Microsite, Species)) + geom_boxplot() + facet_grid(~site)
ggplot(env, aes(Microsite, Even)) + geom_boxplot() + facet_grid(~site)
env %>% filter(abun < 350) %>% ggplot(aes(Microsite, abun)) + geom_boxplot() + facet_grid(~site)


ggplot(rdm, aes(Microsite, RDM)) + geom_boxplot() + facet_grid(~Region)
ggplot(env, aes(Microsite, Species)) + geom_boxplot() + facet_grid(~Region)
ggplot(env, aes(Microsite, Even)) + geom_boxplot() + facet_grid(~Region)
env %>% filter(abun < 350) %>% ggplot(aes(Microsite, abun)) + geom_boxplot() + facet_grid(~Region)


m4 <- glm(Species ~ Microsite + Region/site, family = "poisson", data = env)
summary(m4)
library(glmmTMB)
library(AED)
m5 <- glmmTMB(abun ~ Microsite + (1|Region/site), family = "poisson", data = env)
summary(m5)


m.abun <- env %>% filter(abun <350) %>% glm(abun ~ Microsite + Region/site, family = "quasipoisson", data = .)
summary(m.abun)

m.abun <- env %>% filter(abun <350 & Microsite != "larrea") %>% glmmTMB(abun ~ Microsite + (1|Region/site), family = "nbinom2", data = .)
summary(m.abun)


#let's do some ephedra only
eph <- filter(env, Microsite != "larrea")
rdm <- filter(rdm, Microsite != "larrea")
ggplot(eph, aes(Microsite, Species)) + geom_boxplot() + facet_grid(~Region) + theme_Publication() + ylab("Morphospecies Richness")
eph %>% filter(abun <350) %>% ggplot(aes(Microsite, abun)) + geom_boxplot() + facet_grid(~Region) + theme_Publication() + ylab("Morphospecies Abundance")
ggplot(rdm, aes(Microsite, RDM)) + geom_boxplot() + facet_grid(~Region) + theme_Publication() + ylab("Residual Dry Matter (g)")


m.abun <- eph %>% filter(abun <350) %>% glmmTMB(abun ~ Microsite + RDM + (1|Region/site), family = "nbinom2", data = .)
summary(m.abun)

m2 <- glmmTMB(Species ~ Microsite + RDM + (1|Region/site), family = "poisson", data = eph)
summary(m2)

shapiro.test(rdm$RDM)
m3 <- glmmTMB(RDM ~ Microsite + (1|Region/site), family = "gaussian", data = rdm)
shapiro.test(residuals(m3))
summary(m3)

ggplot(rdm, aes(Microsite, bur.drip)) + geom_boxplot() + facet_grid(~Region) + theme_Publication() + ylab("Number of burrows")

s1 <- glmmTMB(bur.drip ~  RDM + (1|Region/site), family = "poisson", data = eph)
summary(s1)

#interesting stuff is happening, let's check out some ordinations
row.names(comm) <- comm$X
comm <- comm[-1]
row.names(env) <- env$X
env <- env[-1]
all.equal(row.names(comm), row.names(env))
#yay

remove <- c(as.character(factor(eph$X)))
row.names(comm) <- comm$X
keep <- which(row.names(comm) %in% remove)
keep
comm <- comm[keep,]
all.equal(row.names(comm), row.names(eph))


row.names(comm) <- comm$X
comm <- comm[-1]
row.names(eph) <- eph$X
eph <- eph[-1]
eph <- eph[match(row.names(comm), row.names(eph)),]
all.equal(row.names(comm), row.names(eph))



dist <- vegdist(comm)
d1 <- adonis(dist ~ RDM + site + Region + Microsite, data = env)
summary(d1)
d1
# calculate Bray-Curtis distance among samples
comm.bc.dist <- vegdist(comm, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
# plot cluster diagram
plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")
comm.bc.mds <- metaMDS(comm, dist = "bray")
stressplot(comm.bc.mds)
ordiplot(comm.bc.mds, display = "sites", type = "text")
ordipointlabel(comm.bc.mds)

mds.fig <- ordiplot(comm.bc.mds, type = "none")
# plot just the samples, colour by habitat, pch=19 means plot a circle
points(mds.fig, "sites", pch = 19, col = "green", select = eph$Microsite == 
         "open")
points(mds.fig, "sites", pch = 19, col = "blue", select = eph$Microsite == 
         "ephedra")
points(mds.fig, "sites", pch = 19, col = "red", select = eph$Microsite == 
         "open" & env$Region == "Mojave")
points(mds.fig, "sites", pch = 19, col = "orange", select = env$Microsite == 
         "ephedra" & env$Region == "Mojave")

points(mds.fig, "sites", pch = 19, col = "purple", select = env$Microsite == 
         "open" & env$Region == "Cuyama")
points(mds.fig, "sites", pch = 19, col = "black", select = env$Microsite == 
         "ephedra" & env$Region == "Cuyama")



points(mds.fig, "sites", pch = 19, col = "red", select = env$Microsite == 
         "larrea")
# add confidence ellipses around habitat types
ordiellipse(comm.bc.mds, metadata$habitat, conf = 0.95, label = TRUE)
# overlay the cluster results we calculated earlier
ordicluster(comm.bc.mds, comm.bc.clust, col = "gray")

example_NMDS=metaMDS(comm, # Our community-by-species matrix
                     k=2) # The number of reduced dimensions
example_NMDS=metaMDS(comm,k=2,trymax=100)
stressplot(example_NMDS)


plot(example_NMDS)
colors=c(rep("red",5),rep("blue",5))
ordiplot(example_NMDS,type="n")
#Plot convex hulls with colors baesd on treatment
for(i in unique(env$Microsite)) {
  ordihull(example_NMDS$point[grep(i,env$Microsite),],draw="polygon",
           groups=env$Microsite[env$Microsite==i],col=colors[grep(i,env$Microsite)],label=F) } 
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),
                                            rep("blue",5)),air=0.01,cex=1.25)

#I want an ordination of just ephedra

b <- betadisper(dist, eph$RDM)
summary(b)
b
anova(b)
plot(b)
boxplot(b)
TukeyHSD(b)
permutest(b)
str(eph)




library(betapart)
betapart.core(comm)
comm <- decostand(comm, method = "hellinger")
r1 <- cca(comm ~ Microsite + site, data = eph)
summary(r1)
anova(r1)
anova.cca(r1, by = "margin")

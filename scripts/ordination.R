#ordinations for community composition

library(vegan)
library(ggfortify)
library(ggvegan)
comm <- read.csv("clean_data/comm_eph.csv")
row.names(comm) <- comm$X
comm <- comm[,-1]
env <- read.csv("clean_data/cov_eph.csv")


#bray-curtis

comm.bc.dist <- vegdist(comm, method = "bray")


comm.bc.mds <- metaMDS(comm, trymax = 50)
mds2 <- metaMDS(comm, dist = "bray", trymax = 30)
stressplot(comm.bc.mds)

plot(comm.bc.mds)
comm.bc.mds



gof <- goodness

a1 <- adonis(comm.bc.dist ~ Microsite * RDM + rdm.cov, data = env)
a1


#comm <- decostand(comm, method = "hellinger")
#not transforming for CCA
r1 <- cca(comm ~ Microsite + rdm.cov + RDM + site, data = env)
summary(r1)
r1
goodness(r1)
a1 <- anova(r1, by = "margin")
anova(r1, by = "axis")
a1


autoplot(r1)

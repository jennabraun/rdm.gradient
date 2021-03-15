#ordinations for community composition

library(vegan)
library(ggfortify)
library(ggvegan)

comm <- read.csv("clean_data/comm.csv")
row.names(comm) <- comm$X
comm <- comm[,-1]
env <- read.csv("clean_data/cov.csv")


env <- left_join(env, effect, by = "site")
env <- distinct(env)
#bray-curtis

comm.bc.dist <- vegdist(comm, method = "bray")


comm.bc.mds <- metaMDS(comm, trymax = 50)
mds2 <- metaMDS(comm, dist = "bray", trymax = 30)
stressplot(comm.bc.mds)

plot(comm.bc.mds)
comm.bc.mds



a1 <- adonis(comm.bc.dist ~ Microsite * RDM + rdm.cov, data = env)
a1


#comm <- decostand(comm, method = "hellinger")
#not transforming for CCA
r1 <- cca(comm ~ Microsite  + rdm.cov + RDM +  site,  data = env)
summary(r1)
r1
alias(r1, names = TRUE)
goodness(r1)
a1 <- anova(r1, by = "margin")
a1
a2 <- anova(r1, by = "terms")
a2

anova(r1, by = "axis")

test <- autoplot(r1, layers = "species", data = env) + theme_bw()

test + aes(r1, colour = "Microsite")

pca_fort <- fortify(r1, display = "sites") %>%
  bind_cols(env)

pca_fort

ggplot(pca_fort, aes(x = CCA1, y = CCA2, colour = Microsite, shape = Region, size = Region)) + geom_point() +
  scale_color_viridis_d() +
  scale_size(range = 2) +
  coord_equal() + theme_bw() + scale_size_manual(values = c(1, 2, 3)) + xlab("CCA1 (23%)") + ylab("CC2 (17%)")




#test for dissimilarity
a1 <- anosim(comm, env$Microsite, permutations = 999, distance = "bray")
summary(a1)
a2 <- anosim(comm, env$Region, permutations = 999, distance = "bray")
summary(a2)

ggbiplot(r1)

#calculating species overlap and comparing distances

open <- filter(eph, Microsite == "open")
shrub <- filter(eph, Microsite == "ephedra")
comm.open <- comm[match(open$uniID, rownames(comm)),]
comm.shrub <- comm[match(shrub$uniID, rownames(comm)),]
all.equal(rownames(comm.open), as.character(open$uniID))
all.equal(rownames(comm.shrub), as.character(shrub$uniID))

open.dist <- vegdist(comm.open, method = "bray")
shrub.dist <- vegdist(comm.shrub, method = "bray")
open.dist <- as.matrix(open.dist)

#indicator species analysis
library(indicspecies)
factor(env$Microsite, levels = c("ephedra", "open"))
env$group <- paste(env$Region, env$Microsite)
indval <- multipatt(comm, env$group, control = how(nperm=999))
summary(indval)

adonis(comm ~ Microsite  + rdm.cov + RDM + arid + ESI, data = env)

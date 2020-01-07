## Bootstrap RII
## Alex Filazzola
## August 24 2016
## Description
# RII is an effect estimate proposed by Armas et al 2004 to compare a plants performance with neighbours and without.
# It was generated from a similar index that compares a subjects success between control and treatment plots (t-c)/(t+c)
# Heterogeneity stats are typically available for comparisons, but for simpler, within study examples, permutations tests can be used with parametric statistics. 
# Michalet et al. 2014 conducted these analyses in a 2014 Functional Ecology paper between cushion and open microsites. 

##Citations
#Michalet, R., Schöb, C., Lortie, C. J., Brooker, R. W., & Callaway, R. M. (2014). Partitioning net interactions among plants along altitudinal gradients to study community responses to climate change. Functional Ecology, 28(1), 75-86.
#Armas, C., Ordiales, R., & Pugnaire, F. I. (2004). Measuring plant interactions: a new comparative index. Ecology, 85(10), 2682-2686.

## This function calculates RII by subsetting the data and randomizing the treatment and control plots. It will work with unequal samples.
# x = dataframe
# treatment = column with treatment and control variables
# treatment.var = treatment specified in quotations
# control.var = control specified in quotations
# variable = response variable of interest
# number of times to commit bootstrapping. 
se <- function(x) sd(x)/sqrt(length(x)) ## SE

boot.rii <- function(x, treatment,variable,iteration){
  s1 <- filter(x, Microsite == "ephedra") ## subset the treatment group
  o1 <- filter(x, Microsite == "open")  ## subset the control group
  min.samp <- min(length(s1[,1]),length(o1[,1]))  ## minimum number of samples
  rii.avg.total <- c() ## set up blank mean vector
  rii.se.total <- c() ## set up blank mean vector
  for (i in 1:iteration){ ##loop the sampling of treatment and control groups and calculate RII
    set.seed(i) ## control randomization to return same values
    treat.samp<- sample(s1[,variable],min.samp)
    control.samp<-sample(o1[,variable],min.samp)
    return1 <- (treat.samp - control.samp) / (treat.samp+control.samp)
    return1[is.na(return1)] <- 0
    rii.avg <- mean(return1)
    rii.se <- se(return1)
    rii.avg.total <- c(rii.avg.total,rii.avg) ## bind all the means together
    rii.se.total <- c(rii.se.total,rii.se) ## bind all the confidence intervals together
  }
  rii.avg <- mean(rii.avg.total)
  rii.se <- mean(rii.se.total)
  treat <- c(treatment)
  rii.results <- data.frame(factor=treat,average=rii.avg,error=rii.se)
  #print(rii.avg.total)
}


richness.effect <- eph %>% split(.$site) %>% map(boot.rii, treatment="Microsite", variable="Species", 9999)
richness.effect <- do.call(rbind, richness.effect)
richness.effect$site <- row.names(richness.effect)

abun.effect <- eph %>% split(.$site) %>% map(boot.rii, treatment="Microsite", variable="abun", 9999)
abun.effect <- do.call(rbind, abun.effect)
abun.effect$site <- row.names(abun.effect)

eph.rdm <- filter(rdm, Microsite != "larrea")
rdm.effect <- eph.rdm %>% split(.$site) %>% map(boot.rii, treatment="Microsite", variable="RDM", 9999)
rdm.effect <- do.call(rbind, rdm.effect)
rdm.effect$site <- row.names(rdm.effect)

richness.effect$response <- c("richness")
abun.effect$response <- c("abun")
rdm.effect$response <- c("rdm")


effect <- rbind(richness.effect, abun.effect, rdm.effect)


sites <- read.csv("documents/regional_sites.csv")
sites$site <- sites$site.name
effect <- left_join(effect, sites, by = "site")
effect$arid <- effect$aridity.demartonne.annual

write.csv(effect, "clean_data/rii.csv")

#calculate basic rii for rdm, that way have a value for individual shrubs
#this is just for epehdra
rdm <- filter(rdm, Microsite != "larrea")
data <- rdm
y <- "RDM"
rii <- function(data, y)
{
  # Parse out shrub and open
  s1 <- subset(data, Microsite == "ephedra", select=RDM)
  o1 <- subset(data, Microsite == "open", select=RDM)
  # Rii formula
  return1 <- (s1 - o1) / (s1+o1)
}

rdm$uniID <- paste(rdm$site, rdm$ID)
#first make sure datasets matchup nicely
s <- filter(rdm, Microsite == "ephedra")
o <- filter(rdm, Microsite == "open")            
all.equal(s$uniID, o$uniID)
#wicked

rii.rdm <- rii(rdm, "RDM")
rii.rdm <- rename(rii.rdm, rii = RDM)
x1 <- rdm[seq(1, nrow(rdm), by = 2),]
rii <- cbind(rii.rdm, x1)

write.csv(rii, "clean_data/rii_individual.csv")

z <- "rdm.cov"
rii.cov <- function(data, z)
{
  # Parse out shrub and open
  s1 <- subset(data, Microsite == "ephedra", select="rdm.cov")
  o1 <- subset(data, Microsite == "open", select="rdm.cov")
  # Rii formula
  return1 <- (s1 - o1) / (s1+o1)
}


#one of the covers is NA


rdm <- dplyr::filter(rdm, X != 359 & X!=360)
mean(rdm$rdm.cov)

rii.rdm.cov <- rii.cov(rdm, "rdm.cov")
rii.rdm.cov <- rename(rii.rdm.cov, rii = rdm.cov)
x1 <- rdm[seq(1, nrow(rdm), by = 2),]
rii <- cbind(rii.rdm.cov, x1)

write.csv(rii, "clean_data/rii_individual_rdmcover.csv")

#sem modelling

library(tidyverse)
library(lavaan)
env <- read.csv("clean_data/cov.csv")
eph <- filter(env, Microsite != "larrea")
eph$RDM <- as.numeric(as.character(eph$RDM))
Y <- eph$Species
X <- eph$Microsite
X <- factor(X)
#X <- ordered(X)
M <- eph$RDM
Data <- data.frame(Y, X, M)
model <- ' # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M
           # indirect effect (a*b)
             ab := a*b
           # total effect
             total := c + (a*b)
         '
fit <- sem(model, data = Data)
summary(fit, fit.measures = T)

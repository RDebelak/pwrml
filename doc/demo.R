## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(eval=TRUE)

## ----message=FALSE------------------------------------------------------------
library(pwrml)
library(mirt)
library(dplyr)

# To install the most recent mirt version (currently, some matrices may otherwise not be calculated correctly):
# library(devtools)
# install_github('philchalmers/mirt')

## -----------------------------------------------------------------------------
dat <- expand.table(LSAT7)
mirtfit <- mirt(dat,1,verbose = FALSE)

## -----------------------------------------------------------------------------
hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = mirtfit)

## -----------------------------------------------------------------------------
ncps <- calculate_ncps(hyp=hyp)

## -----------------------------------------------------------------------------
power(hyp=hyp,ncp=ncps,alpha=.05,ssize=500)
ssize(hyp=hyp,ncp=ncps,alpha=.05,power=.80)

## -----------------------------------------------------------------------------
altpars <- list(
        a = rlnorm(5,sdlog = .4),
        d = rnorm(5)
        )
altpars
hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)


## -----------------------------------------------------------------------------
group1 = group2 <- list(
        a = rlnorm(5,sdlog = .2),
        d = rnorm(5)
        )

group2$a[1] = (group2$a[1])^2
group2$d[1] = group2$d[1] + .5

altpars <- list(group1,group2)

altpars

hyp <- setup_hypothesis(type = "DIF2PL", altpars = altpars)


## -----------------------------------------------------------------------------
ncps <- calculate_ncps(hyp=hyp)

## -----------------------------------------------------------------------------
ncps <- calculate_ncps(hyp=hyp,sampling=TRUE,sampling.npers = 10^4,approx.npers=10^4)

## -----------------------------------------------------------------------------
n = seq(100,2000)
pow = power(hyp=hyp,ncp=ncps["Gradient"],alpha=.05,ssize=n)
plot(n,pow)

## -----------------------------------------------------------------------------
altpars <- list(
        a = rlnorm(5,sdlog = .4),
        d = rnorm(5)
        )

hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)

data <- setup.data(hyp=hyp,n=500)

## -----------------------------------------------------------------------------
fitted <- mml.fit(data = data,hyp = hyp)

## -----------------------------------------------------------------------------
fitted <- mml.fit(data = data,hyp = hyp,infmat.unres = "ApproxFisher",infmat.res="ApproxFisher",approx.npers = 10^4)

## -----------------------------------------------------------------------------
stats_obs <- stat_obs(fitted)
stats_obs
pvals <- pchisq(stats_obs,df=nrow(hyp$resmod$Amat),ncp=0,lower.tail=FALSE)
pvals

## -----------------------------------------------------------------------------
altpars <- list(
        a = rlnorm(5,sdlog = .4),
        d = rnorm(5)
        )

hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)

calctime(hyp,n.items=7)

## -----------------------------------------------------------------------------
altpars <- list(
        a = rlnorm(7,sdlog = .4),
        d = rnorm(7)
        )

hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)

calctime(hyp,n.items=7)



## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(eval=TRUE)

## ----message=FALSE------------------------------------------------------------
library(pwrml)
library(mirt)
library(dplyr)

# To install the most recent mirt version (results may otherwise be flawed):
# library(devtools)
# install_github('philchalmers/mirt')

## -----------------------------------------------------------------------------
  res = function(altpars,nullpars = NULL) {

    n.items = length(altpars[[1]])

    re = list(
      n.items = n.items,
      itemtype = "2PL",
      Amat = c(0,1,rep(0,(n.items-1)*2)) %>% matrix(.,ncol=n.items*2,byrow=TRUE),
      cvec = 0,
      model = mirt::mirt.model(paste('F = 1-',n.items,'
                           FIXED = (1, d)
                           START = (1,d,0)'))
    )
    return(re)
  }


## -----------------------------------------------------------------------------
unres = function(altpars) {

    re = list(
      parsets = altpars,
      model = 1,
      itemtype = "2PL",
      longpars = pars.long(pars = altpars,itemtype="2PL")
    )

    return(re)
  }

## -----------------------------------------------------------------------------
  maximizeL = function(hyp) {
    # Hypothesis-specific algorithm to find the maximum likelihood restricted parameter set


    maxlpreload= function(pars) {
      # returns the density for each response pattern under the model parameters pars

      patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))

      pre = c()
      for (i in 1:nrow(patterns)) {
        pre[i] = g(patterns[i,],pars)
      }

      return(pre)
    }


    maxl = function(x,pars,pre) {
      # calculates the likelihood of parameters x given model "pars"
      patterns = as.matrix(expand.grid(lapply(1:length(pars$a),function(x) c(0,1))))

      x = list(
        a=c(x,pars$a[2:length(pars$a)]),
        d=c(0,pars$d[2:length(pars$d)])
        )

      res  = c()
      for (i in 1:nrow(patterns)) {
        px = pre[i]
        qx = g(patterns[i,],x)
        res[i] =  {px*log(qx)}
      }
      re = -sum(res)
    }
    resmod = hyp$resmod
    unresmod = hyp$unresmod

    pars = unresmod$parsets
    load.functions(unresmod$itemtype)

    startval = pars$a[1]

    maxlpre = maxlpreload(pars)

    optpar = optim(startval,function(x) {maxl(x,pars,maxlpre)},method = "BFGS")
    re = pars
    re$a = c(optpar$par[1],pars$a[2:length(pars$a)])
    re$d = c(0,pars$d[2:length(pars$d)])

    return(re)
  }


## -----------------------------------------------------------------------------
h_basic = list(
  res = res,
  unres = unres,
  maximizeL = maximizeL
)

## -----------------------------------------------------------------------------
altpars <- list(
        a = rlnorm(5,sdlog = .4),
        d = rnorm(5)
        )

altpars$d[1]=.2

hyp <- setup_hypothesis(type = h_basic, altpars = altpars)

ncps <- calculate_ncps(hyp=hyp)

n = seq(100,2000)
pow = power(hyp=hyp,ncp=ncps["Gradient"],alpha=.05,ssize=n)
plot(n,pow)



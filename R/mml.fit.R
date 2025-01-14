

#' Fit restricted and unrestriced model by ML using mirt
#'
#' @param hyp Hypothesis Object created from the setup_hypothesis function
#' @param infmat.unres String, type of the information matrix for the unrestricted model, "Fisher" or "ApproxFisher" are currently implemented
#' @param infmat.res String, type of the information matrix for the restricted model, "Fisher" or "ApproxFisher" are currently implemented
#' @param free_mean boolean, option to estimate free means between groups
#' @param approx.npers integer, sample size for approximating the Fisher expected information matrix
#' @param data dataframe with data to be fitted
#'
#' @return
#' @export
#'
#' @examples
#'
#' altpars <- list(
#' a = rlnorm(5,sdlog = .4),
#' d = rnorm(5))
#' hyp <- setup_hypothesis(type = "1PLvs2PL", altpars = altpars)
#' data <- setup.data(hyp=hyp,n=500)
#' fitted <- mml.fit(data = data,hyp = hyp)
#'
mml.fit = function(hyp,data,infmat.unres = "Fisher",infmat.res="Fisher",free_mean = FALSE,approx.npers=10^6) {

  if (!is.data.frame(data)) {
    group = data$group
    data = data$data
  }

  multigroup.unres = isTRUE(hyp$unresmod$multigroup)
  multigroup.res = isTRUE(hyp$resmod$multigroup)

  if (isTRUE(free_mean)) {invariance = "free_mean"} else {invariance = ""}

  if (multigroup.unres) {
    unres = mirt::multipleGroup(data,model = hyp$unresmod$model,itemtype = hyp$unresmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance=invariance)
  }

  if (!multigroup.unres) {
    unres =  mirt::mirt(data,model = hyp$unresmod$model,itemtype = hyp$unresmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
  }

  if (multigroup.res) {
    res = mirt::multipleGroup(data,model = hyp$resmod$model,itemtype = hyp$resmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE,group=group,invariance=invariance)
  }

  if (!multigroup.res) {
    res =  mirt::mirt(data,model = hyp$resmod$model,itemtype = hyp$resmod$itemtype,technical = list(NCYCLES = 5000),verbose = FALSE)
    }

  # Caluclating Infmats

  parsr = coef_short(res,itemtype = hyp$resmod$itemtype)
  pars = coef_short(unres,itemtype = hyp$unresmod$itemtype)

    unres@vcov = (infmat(pars,method=infmat.unres,approx.npers = approx.npers,multigroup = multigroup.unres,model = hyp$unresmod$model) * nrow(data) )%>% solve()

    res@vcov = (infmat(parsr,method=infmat.unres,approx.npers = approx.npers,multigroup = multigroup.unres,model = hyp$unresmod$model) * nrow(data) )%>% solve()


    re = list(unres=unres,res=res,hyp=hyp)
  return(re)
}


#' Calculate statistics from fitted mirt models
#'
#' @param fitted Object created by mml.fit function.
#'
#' @return
#' @export
#'
#' @examples
stat_obs = function(fitted,stat=c("Wald","LR","Score","Gradient")) {

  a = rep(NA,4) # Prepare Result Vector

  lx = NULL # Prepare empty lx Vector for Gradient stat

  if ("Wald" %in% stat) {a[1] = wald_obs(fitted)}
  if ("LR" %in% stat) {a[2] = lr_obs(fitted)}
  if ("Score" %in% stat) {
    sctemp = score_obs(fitted,export.lx = TRUE)
    a[3] =sctemp$stat
    lx = sctemp$lx
  }
  if ("Gradient" %in% stat) {a[4] = grad_obs(fitted,lx=lx)}

  re = a[!is.na(a)]
  names(re) = stat

  return(re)
}


wald_obs = function(fitted) {
  resmod = fitted$hyp$resmod
  pars=pars.long(pars = fitted$unres,itemtype = resmod$itemtype,from.mirt=TRUE)
  A = resmod$Amat

  dif = A%*%pars-resmod$cvec
  sigma=mirt::vcov(fitted$unres)
  re = t(dif) %*% solve(A%*% sigma %*% t(A)) %*% dif %>% as.numeric()
  return(re)
}

lr_obs = function(fitted) {
  re = 2*(mirt::logLik(fitted$unres)-mirt::logLik(fitted$res))
  return(re)
}


score_obs = function(fitted,export.lx=FALSE) {

  resmod = fitted$hyp$resmod
  unresmod = fitted$hyp$unresmod
  parsr = coef_short(fitted$res,itemtype = resmod$itemtype)
  load.functions(resmod$itemtype)
  patterns = extract.mirt(fitted$res, "data")
  hashs = apply(patterns,1,digest::digest) %>% factor(.,levels=unique(.))
  rownames(patterns)=hashs
  up = unique(patterns)
  freq = table(hashs)

  if (isTRUE(resmod$multigroup))  { # multigroup model

    freq1 = table(hashs[1:(nrow(patterns)/2)])
    freq2 = table(hashs[(nrow(patterns)/2+1):nrow(patterns)])

    parsr = parsr[[1]]

    ly = lapply(rownames(up),function(i) {l = ldot(up[i,],parsr); list(freq[i] * l,freq1[i] * l,freq2[i] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) %>% array(.,dim=c(length(.),1))


  } else { # single group model

    ly = lapply(rownames(up),function(i) ldot(up[i,],parsr)*freq[i] )
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))
  }

  sigma = mirt::vcov(fitted$res)

  re = t(lx) %*% sigma %*% lx %>% as.numeric()

  if (isTRUE(export.lx)) {
    re = list(stat=re,lx=lx)
  }

  return(re)
}

grad_obs = function(fitted,lx=NULL) {

  resmod = fitted$hyp$resmod
  unresmod = fitted$hyp$unresmod
  parsr = coef_short(fitted$res,itemtype = resmod$itemtype)
  load.functions(resmod$itemtype)
  patterns = extract.mirt(fitted$res, "data")
  hashs = apply(patterns,1,digest::digest) %>% factor(.,levels=unique(.))
  rownames(patterns)=hashs
  up = unique(patterns)
  freq = table(hashs)

  multigroup = isTRUE(unresmod$multigroup)

  if (multigroup & is.null(lx)) { # Multigroup Model

    freq1 = table(hashs[1:(nrow(patterns)/2)])
    freq2 = table(hashs[(nrow(patterns)/2+1):nrow(patterns)])

    parsr = parsr[[1]]

    ly = lapply(rownames(up),function(i) {l = ldot(up[i,],parsr); list(freq[i] * l,freq1[i] * l,freq2[i] * l)} )
    lx = lapply(ly,function(x) x[[1]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx1 = lapply(ly,function(x) x[[2]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))
    lx2 = lapply(ly,function(x) x[[3]]) %>% do.call(rbind,.) %>% colSums() %>% array(.,dim=c(length(.),1))

    lx = c(lx1[1:2,],lx[3:length(lx),],lx2[1:2,]) %>% array(.,dim=c(length(.),1))

  }

  if (!multigroup & is.null(lx)) { # Single Group Model
    ly = lapply(rownames(up),function(i) ldot(up[i,],parsr)*freq[i] )
    lx = do.call(rbind,ly) %>% colSums() %>% array(.,dim=c(length(.),1))
  }

  pars=pars.long(pars = fitted$unres,itemtype = resmod$itemtype,from.mirt=TRUE)
  A = resmod$Amat
  dif = (A%*%pars-resmod$cvec)

  lambda = findlambda(lx,A)

  re = t(lambda) %*% dif %>% as.numeric() %>% abs()

  return(re)
}



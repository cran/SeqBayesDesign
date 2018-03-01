#library(MASS)

### data input ###
#' @export
data.read <- function(filename, exp.type = c("ALT", "CFT.EC"), ...) {
  ### the file is coded that c(stress levels, failure time, censored, weight)
  if(class(filename) == "character"){
    res <- read.csv(file = filename, ...)
  }else if(class(filename) == "data.frame"){
    res <- filename
  }
  if(exp.type == "ALT") {
    colName <- colnames(res)
    res[,1] <- as.numeric(res[,1])
    res[,2] <- as.numeric(res[,2])
    res[,3] <- as.numeric(res[,3])
    res[,4] <- as.numeric(res[,4])
    colnames(res) <- colName
    rownames(res) <- 1:dim(res)[1]
  } else if(exp.type == "CFT.EC") {
    colName <- colnames(res)
    res$Angle <- res$Angle*pi/180
    res <- cbind(res, as.numeric(res[,2] == 1))
    colnames(res) <- c(colName, "UTS")
    rownames(res) <- 1:dim(res)[1]
  }
  return(res)
}
### set data ###
#' @export
dat.setting <- function(dat, exp.type = c("ALT", "CFT.EC"),
                        use.level, max.level = NULL, Cen.time,
                        mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C"),
                        show=TRUE) {
  res <- list()
  colName <- colnames(dat)
  if (exp.type == "ALT") {
    res$data <- dat
    rownames(res$data) <- 1:length(res$data[, 1])
    colnames(res$data) <- c("x", "Y", colName[3:length(colName)])
    res$x <- as.numeric(dat[, 1])
    res$Y <- as.numeric(dat[, 2])
    res$Censored <- dat[, 3]
    res$wts <- dat[, 4]
    res$Ncen <- sum(dat$Censored*dat[, 4])
    res$Nobs <- sum(dat[is.na(dat[, 4]) == 0, 4])
    res$max.level <- as.numeric(max.level)
    res$Cen.time <- as.numeric(Cen.time)
  } else if(exp.type == "CFT.EC") {
    if (sum(dat$UTS == 1) != 0) {
      sigmau <- as.numeric(tapply(X = dat[dat$UTS == 1, 1],
                                  INDEX = dat$Angle[dat$UTS == 1], FUN= mean))
    } else {
      sigmau <- max.level
    }
    dat <- dat[dat$UTS != 1,]
    res$data <- dat
    rownames(res$data) <- 1:length(res$data[, 1])
    colnames(res$data) <- c("x", "Y", colName[3:length(colName)])
    res$x <- dat[, 1]
    res$Y <- dat[, 2]
    res$Censored <- dat[, 3]
    res$wts <- dat[, 4]
    res$sigmau <- sigmau
    res$Rval <- as.numeric(as.character(dat$Rval))[1]
    b <- table(dat$Frequency)
    res$Freq.mod <- as.numeric(names(which.max(b)))
    res$Angle <- dat$Angle[1]
    res$Freq <- dat$Frequency
    res$Ncen <- sum(res$Censored)
    res$Nobs <- sum((dat$UTS != 1)*dat[, 4])
    res$max.level <- sigmau
    res$Cen.time <- Cen.time
  }
  if (length(use.level) == 1) {
    res$use.level$use <- use.level
    res$use.level$wts <- 1
  } else {
    res$use.level$use <- use.level[[1]]
    res$use.level$wts <- use.level[[2]]
  }
  res$mu.fun <- mu.fun
  res$exp.type <- exp.type
  if (show == TRUE) {
    print(list(data = res$data, max.level = res$max.level))}
  invisible(res)
}
################################################################################
#' @export
std.level <- function(use.level, max.level, test.level,
                        mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C")) {
  if(mu.fun == "linear")
  {
    res <- (test.level-use.level$use)/(max.level-use.level$use)
  }else if(mu.fun == "Arrhenius")
  {
    test.level <- test.level + 273.15
    use.level$use <- use.level$use + 273.15
    max.level <- max.level + 273.15
    res <- (1/test.level-1/use.level$use)/(1/max.level-1/use.level$use)
  }else if(mu.fun == "Inv.P")
  {
    res <- (log(test.level)-log(use.level$use))/(log(max.level)-log(use.level$use))
  }else if(mu.fun == "E-C")
  {
    res <- test.level/max.level
  }
  return(res)
}
################################################################################
ori.level <- function(use.level, max.level, st.level,
                        mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C")) {
  if(mu.fun == "linear")
  {
    res <- st.level*(max.level-use.level$use) + use.level$use
  }else if(mu.fun == "Arrhenius")
  {
    use.level$use <- use.level$use + 273.15
    max.level <- max.level + 273.15
    test.level <- st.level*(1/max.level-1/use.level$use)+1/use.level$use
    test.level <- 1/test.level
    res <- test.level - 273.15
  }else if(mu.fun == "Inv.P")
  {
    res <- exp(st.level*(log(max.level)-log(use.level$use))+log(use.level$use))
  }else if(mu.fun == "E-C")
  {
    res <- st.level*max.level
  }
  return(res)
}
################################################################################
std.stress.fun <- function(A, B, stlevels,
                             mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C"),
                             dat) {
  if(mu.fun == "linear" | mu.fun == "Arrhenius" | mu.fun =="Inv.P")
  {
    mu <- A+ B*stlevels
  } else if(mu.fun == "E-C") {
    Rval <- as.numeric(as.character(dat$Rval))
    phi.R <- ifelse(Rval<=1, Rval, 1/Rval)
    gamma.alpha <- 1.6-phi.R*abs(sin(dat$Angle))
    ss <- 1/stlevels
    cc <- (ss-1)*ss^(gamma.alpha-1)*(1-phi.R)^(-gamma.alpha)
    DD <- B*dat$Freq^B*cc+A
    mu <- (log(DD)-log(A))/B
  }
  return(mu)
}
################################################################################
ExpT <- function(pars, stlevels, Cen.time, model = c("lnor","wei"),
                 mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C"), data) {
  A <- pars[1]
  B <- pars[2]
  nu <- pars[3]
  mu <- std.stress.fun(A, B, stlevels, mu.fun, data)
  expTime <- rep(NA, length(stlevels))
  for (i in 1:length(stlevels)) {
    if (model=="lnor") {
      pdf.int <- function(y) {
        normpdf <- dlnorm(y, mu[i], nu)
        res <- y*normpdf
        return(res)
      }
      Time <- integrate(pdf.int, lower=0, upper=Cen.time,
                        stop.on.error = FALSE)[[1]]
      Time <- Time + Cen.time * (1 - plnorm(Cen.time, mu[i], nu))
      expTime[i] = Time
    }
  }
  return(expTime)
}
################################################################################
mu.ders <- function(A, B, stlevels,
                    mu.fun = c("linear", "Arrhenius", "Inv.P", "E-C"), data) {
  if(mu.fun == "linear" | mu.fun == "Arrhenius" | mu.fun == "Inv.P") {
    mu1 <- 1
    mu2 <- stlevels
    mu11 <- 0
    mu22 <- 0
    mu12 <- 0
  } else if(mu.fun == "E-C") {
    phi.R <- ifelse(data$Rval<=1, data$Rval, 1/data$Rval)
    gamma.alpha <- 1.6-phi.R*sin(data$Angle)
    ss <- 1/stlevels
    cc <- (ss-1)*ss^(gamma.alpha-1)*(1-phi.R)^(-gamma.alpha)
    ff <- data$Freq
    DD <- B*ff^B*cc+A
    EE <- B*ff^B*log(ff)*cc+ff^B*cc
    FF <- B*ff^B*(log(ff))^2*cc+2*ff^B*log(ff)*cc

    mu1 <- 1/(B*DD)-1/(A*B)
    mu2 <- log(A)/(B^2)-log(DD)/(B^2)+EE/(B*DD)
    mu11 <- 1/(A^2*B)-1/(B*DD^2)
    mu12 <- 1/(A*B^2)-1/(B^2*DD)-EE/(B*DD^2)
    mu22 <- -2*log(A)/B^3+2*log(DD)/B^3-2*EE/(B^2*DD)-EE^2/(B*DD^2)+FF/(B*DD)
  }
  res <- list(mu1=mu1, mu2=mu2, mu11=mu11, mu12=mu12, mu22=mu22)
  return(res)
}
################################################################################
mnius.loglikelihood <- function(dat, pars, model = c("lnor","wei"),
                                mu.fun = c("linear", "Arrhenius", "Inv.P",
                                           "E-C")) {
  if(mu.fun == "linear"| mu.fun == "Arrhenius" | mu.fun == "Inv.P") {
    A <- pars[1]
    B <- -exp(pars[2])
    nu <- exp(pars[3])
    stlevel <- std.level(dat$use.level, dat$max.level, dat$data[, 1], mu.fun)
    mu <- std.stress.fun(A, B, stlevel, mu.fun, dat)
    wts <- dat$data[, 4]
  } else if(mu.fun == "E-C") {
    A <- exp(pars[1])
    B <- exp(pars[2])
    nu <- exp(pars[3])
    stlevel <- std.level(dat$use.level, dat$max.level, dat$data[, 1], mu.fun)
    mu <- std.stress.fun(A, B, stlevel, mu.fun, dat)
    wts <- rep(1, dat$Nobs)
  }
  cen <- dat$Censored
  y <- dat$data[, 2]
  yy <- log(y)
  if (model=="lnor") {
    normcdf <- ifelse(pnorm(yy, mu, nu) >= 0.99999999, 0.99999999,
                      pnorm(yy, mu, nu))
    normcdf <- ifelse(normcdf <= 1-0.99999999, 1-0.99999999, normcdf)
    normpdf <- ifelse(dnorm(yy, mu, nu)/(y) <= 1-0.99999999, 1-0.99999999,
                      dnorm(yy, mu, nu)/(y))
    ll <- wts*((log(normpdf))*(1-cen)+(log(1-normcdf))*cen)
    res <- (-1)*sum(ll)
  } else if(model=="wei") {
    zz <- (yy-mu)/nu
    ll <- wts*((log(exp(zz-exp(zz)))-log(nu)-yy)*(1-cen) +
                 (log(exp(-exp(zz))))*cen)
    res <- (-1)*sum(ll)
  }
  return(res)
}
################################################################################
#' @export
par.MLE <- function(dat, model, mu.fun, initial, method = "BFGS",
                    hessian = TRUE,...) {
  #call=match.call()
  f <- function(p) {
    mnius.loglikelihood(dat, p, model, mu.fun)
  }
  if (mu.fun == "linear" | mu.fun == "Arrhenius" | mu.fun ==  "Inv.P") {
    init.par <- c(initial[1], log(-initial[2]), log(initial[3]))
    oout <- optim(init.par, f, method = method, hessian = hessian,...)
    coef <- oout$par
    Realcoef <- c(coef[1], -exp(coef[2]), exp(coef[3]))
  } else if (mu.fun == "E-C") {
    init.par <- log(initial)
    oout <- optim(init.par, f, method = method, hessian = hessian,...)
    coef <- oout$par
    Realcoef=exp(coef)
  }
  if (hessian) {
    vcov <- ginv(oout$hessian)
  } else {
    vcov <- NULL
  }
  min <- oout$value
  return(list(est.par = Realcoef, log.likelihood = -min, vcov = vcov))
}
########## Choose responese madel ############
#' @export
lifetime.model <- function(dat, initials) {
  res <- list()
  lnor <- par.MLE(dat, "lnor", mu.fun = dat$mu.fun, initials)
  wei <- par.MLE(dat, "wei", mu.fun = dat$mu.fun, initials)
  if(lnor$log.likelihood > wei$log.likelihood){
    res$lifetime.model <- "Log-normal"
    res$est.par <- lnor$est.par
    res$log.likelihood <- lnor$log.likelihood
  }else{
    res$lifetime.model <- "Weibull"
    res$est.par <- wei$est.par
    res$log.likelihood <- wei$log.likelihood
  }
  return(res)
}
################################################################################
ind.FIM <- function(pars, Cen.time=2000000, dat, stlevel, model, mu.fun) {
  if(model == "lnor") {
    Phi <- pnorm
    phi <- dnorm
    phi1Phi <- function(z) {
      phi(z)/(1-Phi(z))
    }
    dphi1Phi <- function(z) {
      (-z*phi(z)*(1-Phi(z))+phi(z)*phi(z))/(1-Phi(z))^2
    }
    dphiphi <- function(z) {
      return(-z)
    }
    ddphiphi <- function(z) {
      return(-1)
    }
  }else if(model == "wei") {
    Phi <- function(z) {
      return(1-exp(-exp(z)))
    }
    phi <- function(z) {
      return(exp(z-exp(z)))
    }
    phi1Phi <- function(z) {
      return(exp(z))
    }
    dphi1Phi <- function(z) {
      return(exp(z))
    }
    dphiphi <- function(z) {
      return(1-exp(z))
    }
    ddphiphi <- function(z) {
      return(-exp(z))
    }
  }
  #############
  A <- pars[1]
  B <- pars[2]
  nu <- pars[3]
  mu <- std.stress.fun(A, B, stlevel, mu.fun, dat)
  mu.tmp <- mu.ders(A, B, stlevel, mu.fun, dat)
  mu1 <- mu.tmp$mu1
  mu2 <- mu.tmp$mu2
  mu11 <- mu.tmp$mu11
  mu12 <- mu.tmp$mu12
  mu22 <- mu.tmp$mu22
  ######
  zm <- (log(Cen.time) - mu)/nu
  #print(zm)
  Phizm1 <- try(1 - Phi(zm))
  if (class(Phizm1) == "try-error" | Phizm1<1e-8) {
    f11b <- 0
    f22b <- 0
    f33b <- 0
    f12b <- 0
    f13b <- 0
    f23b <- 0
  } else {
    f11b <- mu1^2*dphi1Phi(zm)*Phizm1-nu*mu11*phi(zm)
    f22b <- mu2^2*dphi1Phi(zm)*Phizm1-nu*mu22*phi(zm)
    f33b <- (2*zm*phi1Phi(zm)+zm^2*dphi1Phi(zm))*Phizm1
    f12b <- mu1*mu2*dphi1Phi(zm)*Phizm1-nu*mu12*phi(zm)
    f13b <- mu1*(phi1Phi(zm)+dphi1Phi(zm)*zm)*Phizm1
    f23b <- mu2*(phi1Phi(zm)+dphi1Phi(zm)*zm)*Phizm1
  }
  f11a.int <- function(z) {
    res <- (-1)*mu1^2*ddphiphi(z)*phi(z)+nu*mu11*dphiphi(z)*phi(z)
    return(res)
  }
  f22a.int <- function(z) {
    res <- (-1)*mu2^2*ddphiphi(z)*phi(z)+nu*mu22*dphiphi(z)*phi(z)
    return(res)
  }
  f33a.int <- function(z) {
    res <- (-1)*(1+2*z*dphiphi(z)+z^2*ddphiphi(z))*phi(z)
    return(res)
  }
  f12a.int <- function(z) {
    res <- (-1)*mu1*mu2*ddphiphi(z)*phi(z)+nu*mu12*dphiphi(z)*phi(z)
    return(res)
  }
  f13a.int <- function(z) {
    res <- (-1)*mu1*(dphiphi(z)+z*ddphiphi(z))*phi(z)
    return(res)
  }
  f23a.int <- function(z) {
    res <- (-1)*mu2*(dphiphi(z)+z*ddphiphi(z))*phi(z)
    return(res)
  }

  f11a <- integrate(f11a.int, lower=-Inf, upper=zm,
                    rel.tol = .Machine$double.eps^0.5)[[1]]
  f22a <- integrate(f22a.int, lower=-Inf, upper=zm,
                    rel.tol = .Machine$double.eps^0.5)[[1]]
  f33a <- integrate(f33a.int, lower=-Inf, upper=zm,
                    rel.tol = .Machine$double.eps^0.5)[[1]]
  f12a <- integrate(f12a.int, lower=-Inf, upper=zm,
                    rel.tol = .Machine$double.eps^0.5)[[1]]
  f13a <- integrate(f13a.int, lower=-Inf, upper=zm,
                    rel.tol = .Machine$double.eps^0.5,
                    stop.on.error = FALSE)[[1]]
  f23a <- integrate(f23a.int, lower=-Inf, upper=zm,
                    rel.tol = .Machine$double.eps^0.5)[[1]]

  f11 <- f11a+f11b
  f22 <- f22a+f22b
  f33 <- f33a+f33b
  f12 <- f12a+f12b
  f13 <- f13a+f13b
  f23 <- f23a+f23b

  res <- c(f11, f22, f33, f12, f13, f23)
  return(res)
}

################################################################################
FImatrix <- function(pars, Cen.time, dat, stlevel, model, mu.fun, freq, n.new){
  nu <- pars[3]
  N <- dat$Nobs + n.new  ### Sequential needed
  pi.vec <- c(dat$wts, 1)/N ### wts
  Imat <- matrix(0, 3, 3)
  dat1 <- dat
  for(i in 1:length(stlevel)) {
    pi.val <- pi.vec[i]
    dat1$Freq <- freq[i]
    xres <- ind.FIM(pars, dat1$Cen.time, dat1, stlevel[i], model, dat1$mu.fun)
    #xres=lsinf.CD(A=A, B=B, nu=nu, Nm=Cen.time, sigmai=sigmai, sigmau=sigmau,
    #              Rval=Rval, ff=ff[i], angle=angle)
    tmp.mat <- (1/nu^2)*matrix(c(xres[1], xres[4], xres[5],
                                 xres[4], xres[2], xres[6],
                                 xres[5], xres[6], xres[3]), 3, 3)
    Imat <- Imat+pi.val*tmp.mat
  }
  res <- N*Imat
  return(res)
}
###################avar.fun###################
Eval.seq.criteria <- function(candidate, par.draw, hist.stlevel, dat,
                              use.pattern, quantile, Sinv, model, mu.fun) {
  size <- length(par.draw)/3
  par.draw <- matrix(par.draw, ncol=size)
  if(length(candidate) == 0) {
    DETer <- matrix(rep(NA, 1*size), 1, size)
    TAVar <- matrix(rep(NA, 1*size), 1, size)
  } else {
    DETer <- matrix(rep(NA,length(candidate)*size), length(candidate), size)
    TAVar <- matrix(rep(NA,length(candidate)*size), length(candidate), size)
  }

  for (n.draw in 1:size) {
    A <- par.draw[1, n.draw]
    B <- par.draw[2, n.draw]
    if(length(candidate) == 0){
      seq.Fmat <- FImatrix(par.draw[, n.draw], dat$Cen.time, dat, hist.stlevel,
                           model, mu.fun, dat$Freq, 0)
      DETer[1, n.draw] <- log(det(seq.Fmat))
      dat1 <- dat
      dat1$Freq <- dat1$Freq.mod
      der.mu <- mu.ders(A, B, use.pattern$use, mu.fun, dat1)
      cvec <- matrix(c(der.mu$mu1,der.mu$mu2,
                       rep(qnorm(quantile, mean = 0, sd = 1),
                           length(use.pattern$use))), length(use.pattern$use), 3)
      app.var <- ginv(seq.Fmat)
      Vi <- c()
      for(row in 1:length(use.pattern$use)) {
        Vi <- c(Vi, sum(cvec[row,] %*% app.var*cvec[row,]))
      }
      wavar <- sum(use.pattern$wts*Vi)
      TAVar[1, n.draw] <- ifelse(wavar>0, wavar, NA)
    }else {
      if(is.na(dat$x[1]) == 1) {
        F1=Sinv
      } else {
        F1 <- FImatrix(par.draw[, n.draw], dat$Cen.time, dat, hist.stlevel,
                       model, mu.fun, dat$Freq, 1)
        #F1=Seq.Fmat.CD(A, B, nu, Nm=2000000, dat1.NEW, q.vec=dat1.NEW$x/dat1.NEW$sigmau, dat1.NEW$Freq)
      }
      for (n.new in 1:length(candidate)) {
        seq.Fmat <- c()
        F2 <- FImatrix(par.draw[, n.draw], dat$Cen.time, dat, candidate[n.new],
                       model, mu.fun, dat$Freq.mod, 1)
        #F2=Seq.Fmat.CD(A, B, nu, Nm=2000000, dat1.NEW, candidate[n.new], dat1.NEW$Freq.mod)

        #seq.Fmat <- F1+F2+ginv(Sinv)
        seq.Fmat=F1+F2+Sinv
        DETer[n.new, n.draw] <- log(det(seq.Fmat))
        dat1 <- dat
        dat1$Freq <- dat1$Freq.mod
        der.mu <- mu.ders(A, B, use.pattern$use, mu.fun, dat1)
        cvec <- matrix(c(der.mu$mu1,der.mu$mu2,
                         rep(qnorm(quantile, mean = 0, sd = 1),
                             length(use.pattern$use))),
                         length(use.pattern$use), 3)
        app.var <- ginv(seq.Fmat)
        Vi <- c()
        for(row in 1:length(use.pattern$use)) {
          Vi <- c(Vi, sum(cvec[row,] %*% app.var*cvec[row,]))
        }
        wavar <- sum(use.pattern$wts*Vi)
        TAVar[n.new, n.draw] <- ifelse(wavar>0, wavar, NA)
      }
    }
  }
  res.D <- rowMeans(DETer, na.rm=T)
  res.AVar <- rowMeans(TAVar, na.rm=T)
  return(list(res.D = res.D, res.AVar = res.AVar, TAVar = TAVar))
}
### next point #####
next.point <- function(candidate, par.sample, dat, model, quantile = 0.1, lambda, prior, Sinv){
  #Sinv=diag(c(12/(prior[2]-prior[1])^2,12/(prior[4]-prior[3])^2,(prior[5]-1)^2*(prior[5]-2)/prior[6]^2))
  hist.stlevel = as.vector(t(std.level(dat$use.level, dat$max.level,
                                         dat$data[1], dat$mu.fun)))
  Eval=Eval.seq.criteria(candidate, par.sample, hist.stlevel, dat, dat$use.level, quantile,
                         Sinv, model, dat$mu.fun)

  RES <- matrix(rep(NA, 3*length(candidate)), 3, length(candidate))
  RES[1,] <- Eval$res.D
  RES[2,] <- Eval$res.AVar
  row.names(RES) <- c("D-optimality",  "C-optimality", "Dual-optimality")
  colnames(RES) <- as.factor(candidate)
  opt.point.D = candidate[which(Eval$res.D == max(Eval$res.D))]
  opt.point.C = candidate[which(Eval$res.AVar == min(Eval$res.AVar))]

  opt.cri.D = Eval$res.D[which(Eval$res.D == max(Eval$res.D))]
  opt.cri.C = Eval$res.AVar[which(Eval$res.AVar == min(Eval$res.AVar))]

  res.Dual <- lambda*(Eval$res.D)/opt.cri.D+(1-lambda)*opt.cri.C/(Eval$res.AVar)
  RES[3,] <- res.Dual
  opt.point.Dual = candidate[which(res.Dual == max(res.Dual))]
  opt.cri.Dual = res.Dual[which(res.Dual == max(res.Dual))]

  if (lambda == 1) {
    next.point = data.frame(next.point = opt.point.D, log.det = opt.cri.D,
                            avar =  Eval$res.AVar[which(Eval$res.D == max(Eval$res.D))],
                            row.names = "D-optimality")
  }else if(lambda == 0) {
    next.point = data.frame(next.point = opt.point.C,
                            log.det = Eval$res.D[which(Eval$res.AVar == min(Eval$res.AVar))],
                            avar =  opt.cri.C,
                            row.names = "C-optimality")
  }else {
    next.point = data.frame(next.point = opt.point.Dual,
                            log.det = Eval$res.D[which(res.Dual == max(res.Dual))],
                            avar =  Eval$res.AVar[which(res.Dual == max(res.Dual))],
                            row.names = "Dual-optimality")
  }
  return(list(next.point = next.point, eval = RES))
}

### Generate data #####
gen.data.EC <- function(partrue, st.next.design, hist.data, model, mu.fun, exp.type)
{
  if(exp.type == "CFT.EC"){
    hist.data$Freq <- hist.data$Freq.mod
    stress.fun.NEW <- std.stress.fun(partrue[1],partrue[2], st.next.design,
                                       mu.fun, hist.data)
    if(model == "lnor") {
      next.cycle.NEW <- rlnorm(1, stress.fun.NEW, partrue[3])
    }else if(model == "wei") {
      next.cycle.NEW <- rweibull(1, 1/partrue[3], exp(stress.fun.NEW))
    }
    next.cen.NEW <- ifelse(next.cycle.NEW > hist.data$Cen.time, 1, 0)
    next.cycle.NEW=ifelse(next.cycle.NEW > hist.data$Cen.time,
                           hist.data$Cen.time, next.cycle.NEW)
    data.NEW <- rbind(hist.data$data, c(st.next.design*hist.data$sigmau, next.cycle.NEW, next.cen.NEW, 1, hist.data$Rval, hist.data$Freq.mod, hist.data$Angle, 0))
    data.NEW <- data.NEW[!is.na(data.NEW[, 1]),]
    dat1.NEW <- dat.setting(data.NEW, exp.type, hist.data$use.level, hist.data$max.level, hist.data$Cen.time, mu.fun, show = FALSE)
    dat1.NEW$new.data <- data.NEW[length(hist.data$data[,1]) + 1,]
    cat("stlevel is ", st.next.design, "\n")
    cat("Stress level is ",  st.next.design*hist.data$sigmau, "\n")
  }else if(exp.type == "ALT") {
    stress.fun.NEW <- std.stress.fun(partrue[1],partrue[2], st.next.design, mu.fun, hist.data)
    if(model == "lnor") {
      next.cycle.NEW <- rlnorm(1, stress.fun.NEW, partrue[3])
    }else if(model == "wei") {
      next.cycle.NEW <- rweibull(1, 1/partrue[3], exp(stress.fun.NEW))
    }
    next.cen.NEW <- ifelse(next.cycle.NEW > hist.data$Cen.time, 1, 0)
    next.cycle.NEW=ifelse(next.cycle.NEW > hist.data$Cen.time,
                          hist.data$Cen.time, next.cycle.NEW)
    ori <- ori.level(hist.data$use.level, hist.data$max.level, st.next.design, mu.fun)
    data.NEW <- rbind(hist.data$data, c(ori, next.cycle.NEW, next.cen.NEW, 1))
    data.NEW <- data.NEW[!is.na(data.NEW[, 1]),]
    dat1.NEW <- dat.setting(data.NEW, exp.type, hist.data$use.level,
                            hist.data$max.level, hist.data$Cen.time,
                            mu.fun, show = FALSE)
    dat1.NEW$new.data <- data.NEW[length(hist.data$data[,1]) + 1,]
    cat("stlevel is ", st.next.design, "\n")
    cat("Stress level is ", ori, "\n")
  }
  return(dat1.NEW)
}

###################
MCMC.draw.EC <- function(dat, partrue, n.int, model, mu.fun, prior, priorDis,
                         L.lag=200, transp=0.5, show = TRUE) {
  initials <- c(partrue[1], partrue[2], partrue[3])
  print(initials)
  stlevel <- std.level(dat$use.level, dat$max.level, dat$data[,1], mu.fun)
  if(mu.fun == "E-C"){
    FF <- FImatrix(initials, dat$Cen.time, dat, stlevel, model, mu.fun, dat$Freq, 0)
    trans <- matrix(c(-1/initials[1], 0, 1/initials[2], transp*initials[2]^(transp-1)),2,2)
    COV <- ginv(FF)[1:2,1:2]
    Cov <- trans %*% COV %*% t(trans)
    draw1 <- TMCMC(dat, n.int, initials, model, mu.fun, Cov, prior, priorDis, transp)
    draws <- draw1$par
    acp <- c(mean(draw1$acp[1,]), mean(draw1$acp[2,]))
  }else {
    FF <- FImatrix(initials, dat$Cen.time, dat, stlevel, model, mu.fun,
                   rep(0, length(stlevel)), n.new = 0)
    Cov <- ginv(FF)
    draw1 <- TMCMCALT(dat, n.int, initials, model, mu.fun, Cov, prior, priorDis, 0.1)
    draws <- draw1$par
    acp <- c(mean(draw1$acp[1,]), mean(draw1$acp[2,]), mean(draw1$acp[3,]))
  }
  if(show == TRUE){
    par(mfrow=c(1,3))
  }
    yy1 <- acf(draws[1, (n.int/10):n.int], L.lag, plot=show, main="A")
    yy2 <- acf(draws[2, (n.int/10):n.int], L.lag, plot=show, main="B")
    yy3 <- acf(draws[3, (n.int/10):n.int], L.lag, plot=show, main="nu")

  #ci=qnorm((1 + 0.95)/2)/sqrt(yy1$n.used)
  ci <- 0.01
  aa <- c(1:(L.lag+1))
  lag <- max(c(min(c(L.lag, aa[yy1$acf < ci])), min(c(L.lag, aa[abs(yy2$acf) < ci])),
               min(c(L.lag, aa[abs(yy3$acf) < ci]))))
  pts <- seq(n.int/10, n.int, lag)
  par.draw1 <- draws[, pts]
  if(show == TRUE){
    par(mfrow=c(1, 3))
  }
    yy1 <- acf(par.draw1[1, ], L.lag/2, plot=show, main="A")
    yy2 <- acf(par.draw1[2, ], L.lag/2, plot=show, main="B")
    yy3 <- acf(par.draw1[3, ], L.lag/2, plot=show, main="nu")
  cov.matrix <- cov(t(par.draw1))
  invisible(list( lag = lag, ori.draws = draws, samples = par.draw1,
                  n.pts = length(pts), acp=acp, cov = cov.matrix))
}
###################
### Seq_Design ###
#' @export
SBD.sim <- function(N.design, N.data, dat, prior, priorDis, n.int=100000,
                    n.sample=1000, partrue, model, mu.fun=dat$mu.fun,
                    candidate=seq(0.35, 0.75, 0.05), quantile=0.1, lambda=0,
                    L.lag=100, transp=0.5, show.acf = FALSE){
  use.pattern <- dat$use.level
  lag <- c()
  covInf <- list()
  if(priorDis == "uniform"){
    Sinv <- diag(c((prior[2]-prior[1])^2/12, (prior[4]-prior[3])^2/12,
                  prior[6]^2/((prior[5]-1)^2*(prior[5]-2))))
  }else if(priorDis == "normal"){
    Sinv <- diag(c(prior[2]^2, prior[4]^2, prior[6]^2/((prior[5]-1)^2*(prior[5]-2))))
  }
  Bayesian.NEW <- matrix(rep(NA,3*(N.design+1)),3,N.design+1) ### collecting Bayesian estimators
  datNEW <- dat
  hist.stlevel <- as.vector(t(std.level(datNEW$use.level, datNEW$max.level,
                                        datNEW$data[1], datNEW$mu.fun)))

  if (N.data != 0){
    draws <- MCMC.draw.EC(datNEW, partrue, n.int, model, mu.fun, prior, priorDis,
                          L.lag, transp, show = show.acf)
    lag <- draws$lag
    covInf[[1]] <- draws$cov
    #print(lag)
    sample.draw <- sample(seq(1,length(draws$samples[1,])), n.sample)
    par.draw1 <- draws$samples[,sample.draw]
    hist.inf <- Eval.seq.criteria(c(), par.draw1, hist.stlevel, datNEW,
                                 use.pattern, quantile, Sinv, model, mu.fun)
    hist.D <- hist.inf$res.D
    hist.avar <- hist.inf$res.AVar
    Bayesian.NEW[,1] <- c(mean(par.draw1[1,]),mean(par.draw1[2,]),mean(par.draw1[3,]))
    print(hist.avar)
  }else{
    par.draw1 <- matrix(rep(NA, 3*n.int), 3, n.int)
    par.draw1[1,] <- runif(n.int, prior[1], prior[2])
    par.draw1[2,] <- runif(n.int, prior[3], prior[4])
    par.draw1[3,] <- sqrt(1/rgamma(n.int, shape=prior[5], scale=prior[6]))
    Bayesian.NEW[,1] <- c(mean(par.draw1[1,]),mean(par.draw1[2,]),mean(par.draw1[3,]))
    draws <- par.draw1
    covInf[[1]]  <-  cov(t(par.draw1))
  }
  if(N.design == 0){
    invisible(list( hist.avar= hist.avar, ori.draw = draws, samples = par.draw1,
                    Bayes = Bayesian.NEW[,1], cov = covInf))
  }else {
    n.cand <- length(candidate)
    OPT.Summary <- matrix(rep(NA, (N.design+1)*3), N.design+1, 3)
    for (n.design in 1:N.design){
      cat('Design number: ', n.design, "\n")
      sample.draw <- sample(seq(1,length(par.draw1[1,])), n.sample)
      par.draw <- par.draw1[,sample.draw]
      if (N.data != 0) {
        nextpoint <- next.point(candidate, par.draw, datNEW, model, quantile,
                                lambda, prior, Sinv)
      }else {
        if(dat$exp.type == "CFT.EC"){
          data.NEW <- matrix(c(NA, 0, datNEW$Rval, datNEW$Freq.mod, datNEW$Angle, NA, 0),1, 7)
          colnames(data.NEW) <- c("Stress", "Cycles", "Rval", "Freq", "Angle", "Censored" ,"UTS")
          data.NEW <- data.frame(data.NEW)
          dat1.NEW <- dat.setting(data.NEW, data.NEW$exp.type, datNEW$use.level,
                               datNEW$sigmau, datNEW$Cen.time, datNEW$mu.fun)
          #Eval <- Eval_ObjFun(candidate, par.draw, datNEW, use.pattern, quantile, Sinv)
          #Eval <- Eval.seq.criteria(candidate, par.sample, hist.stlevel, dat1.NEW, dat$use.level, quantile, Sinv, dat$model, dat$mu.fun)
        }else{
          data.NEW <- matrix(c(NA, 0, datNEW$Rval, datNEW$Freq.mod, datNEW$Angle, NA, 0),1, 7)
        }
      }
      OPT.Summary[(n.design+1), 1:3] <- as.matrix(nextpoint$next.point)
      nextp <- as.numeric(nextpoint$next.point[1])
      cat(rownames(nextpoint$next.point), "\n")
      ########
      datNEW <- gen.data.EC(partrue, nextp, datNEW, model, datNEW$mu.fun, datNEW$exp.type)
      N.data <- N.data + 1
      ########
      draws <- MCMC.draw.EC(datNEW, partrue, n.int, model, mu.fun, prior, priorDis,
                            L.lag, transp, show = show.acf)
      covInf[[n.design+1]] = draws$cov
      lag <- c(lag, draws$lag)
      #print(draws$lag)
      par.draw1 <- draws$samples
      Bayesian.NEW[, n.design+1] <- c(mean(par.draw1[1,]), mean(par.draw1[2,]),
                                      mean(par.draw1[3,]))
      cat("Bayesian estimates:", round(Bayesian.NEW[, n.design+1], 4), "\n", "\n")
    }
    OPT.Summary[1,] <- c(NA, hist.D, hist.avar)
    row.names(OPT.Summary) <- 0:N.design
    colnames(OPT.Summary) <- c("next.point", "log.det", "avar")
    print(OPT.Summary)
    #cat("\n", "Data set:", "\n")
    #print(dat1.NEW$data)
    invisible(list(Bayes = Bayesian.NEW, Final.opt = OPT.Summary,
                   data = datNEW$data, hist.inf= c(hist.D, hist.avar)))
  }
}
### Experiment-type Seq_Design (no generate data)###
#' @export
SBD.next.pt <- function(N.data, dat, prior, priorDis, n.int = 100000, n.sample = 1000,
                    initial, model, mu.fun=dat$mu.fun, candidate = seq( 0.35, 0.75, 0.05),
                    quantile=0.1, lambda=0, L.lag=100, transp=0.5, show.acf = FALSE) {
  use.pattern <- dat$use.level
  lag <- c()
  Bayesian.NEW <- matrix(rep(NA,3*1), 3, 1) ### collecting Bayesian estimators
  datNEW <- dat
  hist.stlevel <- as.vector(t(std.level(datNEW$use.level, datNEW$max.level,
                                        datNEW$data[1], datNEW$mu.fun)))

  if (N.data != 0){
    draws <- MCMC.draw.EC(datNEW, initial, n.int, model, mu.fun, prior, priorDis,
                          L.lag, transp, show = show.acf)
    lag <- draws$lag
    par.draw <- draws$samples
    Bayesian.NEW <- c(mean(par.draw[1,]), mean(par.draw[2,]), mean(par.draw[3,]))
    hist.inf <- Eval.criteria(Bayesian.NEW, hist.stlevel, datNEW, use.pattern,
                              quantile, model, mu.fun)
    hist.D = hist.inf$res.D
    hist.avar = hist.inf$res.AVar
  }else{
    par.draw <- matrix(rep(NA, 3*n.int), 3, n.int)
    par.draw[1,] <- runif(n.int, prior[1], prior[2])
    par.draw[2,] <- runif(n.int, prior[3], prior[4])
    par.draw[3,] <- sqrt(1/rgamma(n.int, shape=prior[5], scale=prior[6]))
    Bayesian.NEW <- c(mean(par.draw[1,]), mean(par.draw[2,]), mean(par.draw[3,]))
    draws <- par.draw
  }
  OPT.Summary <- matrix(rep(NA, 2*3), 2, 3)
  cat('Next Design:', "\n")
  sample.draw <- sample(seq(1,length(par.draw[1,])), n.sample)
  par.draw <- par.draw[, sample.draw]
  if (N.data != 0) {
    nextpoint <- next.point(candidate, par.draw, datNEW, model, quantile, lambda,
                            prior, draws$cov)
  }else{
    data.NEW <- matrix(c(NA, 0, datNEW$Rval, datNEW$Freq.mod, datNEW$Angle, NA, 0), 1, 7)
    colnames(data.NEW) <- c("Stress", "Cycles", "Rval", "Freq", "Angle", "Censored" ,"UTS")
    data.NEW <- data.frame(data.NEW)
    datNEW <- dat.setting(data.NEW, data.NEW$exp.type, datNEW$use.level, datNEW$sigmau,
                          datNEW$Cen.time, datNEW$mu.fun)
    #Eval <- Eval_ObjFun(candidate, par.draw, datNEW, use.pattern, quantile, Sinv)
  }
  OPT.Summary[2, 1:3] <- as.matrix(nextpoint$next.point)
  nextp <- as.numeric(nextpoint$next.point[1])
  cat(rownames(nextpoint$next.point), "\n")
  data.NEW <- rbind(datNEW$data, c(nextp*dat$sigmau, NA, NA, 1, dat$Rval,
                                   dat$Freq.mod, dat$Angle, 0))
  data.NEW <- data.NEW[!is.na(data.NEW[, 1]),]
  datNEW <- dat.setting(data.NEW, dat$exp.type, dat$use.level, dat$max.level,
                          dat$Cen.time, dat$mu.fun, show = FALSE)
  OPT.Summary[1,] <- c(NA, hist.D, hist.avar)
  row.names(OPT.Summary) <- 0:1
  colnames(OPT.Summary) <- c("next.point", "log.det", "avar")
  print(OPT.Summary)
  #cat("\n", "Data set:", "\n")
  #print(dat1.NEW$data)
  invisible(list( Bayes = Bayesian.NEW, Final.opt = OPT.Summary, data = datNEW$data))
}
################################################################################
use.profile <- function(quseL=0.05, quseU=0.25, pts=20, pattern=1)
{
  xx=seq(0, 1, pts)
  quse=quseL+(quseU-quseL)*xx

  if(pattern==1)
  {
    quse.wts=dbeta(xx, 2, 5)
  }

  if(pattern==2)
  {
    quse.wts=dbeta(xx, 5, 2)
  }

  if(pattern==3)
  {
    quse.wts=dbeta(xx, 3, 3)
  }

  if(pattern==4)
  {
    quse.wts=dbeta(xx, 2,5)+dbeta(xx, 10, 5)
  }

  quse.wts=quse.wts/sum(quse.wts)

  par(mai=c(0.2,0.2,0.2,0.2))
  #barplot(quse.wts, width=(quseU-quseL)/pts, axes=F)

  res=list(use=quse, wts=quse.wts)
  return(res)
}

### Plots ###
#' @export
SBD.plot <- function(dat, Obj, cri = c("C-optimality", "D-optimality"),
                     y.at = seq(0.3, 1, 0.1)) {
  hist.stlevel = std.level(dat$use.level, dat$max.level, dat$x, dat$mu.fun)
  #par(mar=c(8, 5, 5, 1), oma=c(2, 0, 0, 0))
  layout(matrix(c(1,1,1,1,2,2), 2, 3))

  plot(length(dat$data[,1]):length(Obj$data[,1]), Obj$Final.opt[, 1], type='b', col=2,
       pch=19, ylim= c(0, 1), cex=1.5,
       xlim=c(0, length(Obj$data[,1]) + 1), cex.axis=2, xlab='', ylab='',
       xaxt = "n", yaxt = "n")
  axis(side=1, at = 1:length(Obj$data[,1]), cex.axis=1.2)
  axis(side=2, at = y.at, cex.axis=1.5)
  mtext('Number of samples in current data', side=1, line=2, cex=1)
  mtext('Standardized level', side=2, line=2, cex=1.5)
  lines(c(dat$Nobs+0.5, dat$Nobs+0.5), c(-1,2), lty=5)
  #aa=sort(hist.stlevel)
  points(1:length(dat$data[,1]), hist.stlevel, pch=8, cex=1.5)
  legend("bottomleft", c("historical stress level", "sequential design"),
         pch=c(8, 19), cex=1.2, col=c(1,2))
  mtext(cri, side=3, cex=1.5, line = 1)

  if(cri == "C-optimality"){y <- Obj$Final.opt[, 3]}else{y <- Obj$Final.opt[, 2]}
  Nnew <- length(Obj$data[, 1])
  plot(length(dat$data[,1]):length(Obj$data[,1]), y, type='b', col=1, pch=19,
       xlim=c(length(dat$data[,1])-1, Nnew+1), cex=1.5, cex.axis=1.2,
       xlab='', ylab='', main="")
  mtext('Avar', side=2, line=2, cex=1.5)
  mtext(cri, side=3, cex=1.5, line = 1)
  mtext('Sample size', side=1, line=2, cex=1)
}

### MLE non-sequential design ###
### C-optimality ###
### MLE non-sequential design C-opt ###
################################################################################
Eval.criteria <- function(pars, stlevel, dat, use.pattern, quantile, model, mu.fun, Sinv=c()) {
  A <- pars[1]
  B <- pars[2]
  F1 <- FImatrix(pars, dat$Cen.time, dat, stlevel, model, mu.fun, dat$Freq, 0)
  #F1=Fmat.CD(A, B, nu, Nm=2000000, dat, q.vec=dat$x/dat$sigmau, pi.vec=rep(1/length(dat$Y), length(dat$Y)), kk=length(dat$Y), ff=dat$Freq)
  res.D <- log(det(F1))
  dat1 <- dat
  dat1$Freq <- dat1$Freq.mod
  der.mu <- mu.ders(A, B, dat$use.level$use, mu.fun, dat1)
  cvec <- matrix(c(der.mu$mu1, der.mu$mu2,
                     rep(qnorm(quantile, mean = 0, sd = 1), length(dat$use.level$use))),
                   length(dat$use.level$use), 3)
  if(length(Sinv) == 0){
    app.var <- ginv(F1)
    } else {
      app.var <- ginv(F1+Sinv)
    }
  Vi <- c()
  for (row in 1:length(dat$use.level$use)) {
    Vi <- c(Vi, sum(cvec[row,] %*% app.var*cvec[row,]))
  }
  wavar <- sum(dat$use.level$wts*Vi)
  res.AVar <- ifelse(wavar>0, wavar, NA)
  return(list( res.D = res.D, res.AVar = res.AVar))
}

########################################################
trad.design <- function(dat, partrue, N.design = 12, design.type = c("TOD", "EQD"), mu.fun,
                          can.stress = seq(0.35, 0.7, 0.05), stressH = 0.75, N.middle = 0, quantile) {
  design.C <- list()
  design.D <- list()
  design.EQD <- list()
  dq <- matrix(rep(NA, length(can.stress), N.design-N.middle-1), length(can.stress), N.design-N.middle-1)
  cq <- matrix(rep(NA, length(can.stress), N.design-N.middle-1), length(can.stress), N.design-N.middle-1)
  mdp <- c()
  mcp <- c()
  mmdp <- c()
  mmcp <- c()
  if (design.type == "EQD") {
    design.EQD$opt <- c(can.stress, stressH)
    design.EQD$pts <- sort(rep(design.EQD$opt, N.design/length(design.EQD$opt)))
    design.EQD$n <- rep(N.design/length(design.EQD$opt), length(design.EQD$opt))
    par.est = NA
  } else if (design.type == "TOD") {
    par.est <- par.MLE(dat, mnius.loglikelihood, "lnor", dat$mu.fun, starts=log(partrue))$est.par
    for (j in 1:length(can.stress)) {
      for (i in 1:(N.design-N.middle-1)) {
        lower <- can.stress[j]
        stlevel <- c(rep(lower, i), rep((lower + stressH)/2, N.middle), rep(stressH, N.design-i-N.middle))
        dat$Freq <- rep(dat$Freq.mod, N.design)
        dat$Nobs <- N.design
        eval <- Eval.criteria(partrue, stlevel, dat, dat$use.level, quantile, dat$mu.fun)
        cq[j,i] <- eval$res.AVar
        dq[j,i] <- eval$res.D
      }
      position.D <- which(dq[j,] == max(dq[j,]))
      position.C <- which(cq[j,] == min(cq[j,]))

      if (length(position.D) == 0) {dqq = NA} else {dqq = position.D}
      mdp <- c(mdp, dqq)
      mmdp <- c(mmdp, dq[j, mdp[j]])
      if (length(position.C) == 0) {cqq = NA} else {cqq = position.C}
      mcp <- c(mcp, cqq)
      mmcp <- c(mmcp, cq[j, mcp[j]])
    }
    pos.D <- min(which(mmdp == max(mmdp[mmdp>0], na.rm=T)))
    pos.C <- min(which(mmcp == min(mmcp[mmcp>0], na.rm=T)))

    if (is.finite(pos.D) == FALSE) {
      design.D$opt <- rep(NA, 2)
      design.D$pts <- rep(NA, N.design)
      design.D$n <- rep(NA, 2)
    } else {
      nLow <- mdp[pos.D]
      nHigh <- N.design-N.middle-nLow
      stressL <- can.stress[pos.D]
      stressM <- (stressL + stressH)/2
      design.D$pts <- c(rep(stressL, nLow), rep(stressM, N.middle), rep(stressH, nHigh))
      if (N.middle == 0) {
        design.D$opt <- c(stressL, stressH)
        design.D$n <- c(nLow, nHigh)
      } else {
        design.D$opt <- c(stressL, stressM, stressH)
        design.D$n <- c(nLow, N.middle, nHigh)
      }
    }

    if (is.finite(pos.C) == FALSE) {
      design.C$opt <- rep(NA, 2)
      design.C$pts <- rep(NA, N.design)
      design.C$n <- rep(NA, 2)
    } else {
      nLow <- mcp[pos.C]
      nHigh <- N.design-N.middle-nLow
      stressL <- can.stress[pos.C]
      stressM <- (stressL + stressH)/2
      design.C$pts <- c(rep(stressL, nLow), rep(stressM, N.middle), rep(stressH, nHigh))
      if (N.middle == 0) {
        design.C$opt <- c(stressL, stressH)
        design.C$n <- c(nLow, nHigh)
      } else {
        design.C$opt <- c(stressL, stressM, stressH)
        design.C$n <- c(nLow, N.middle, nHigh)
      }
    }
  }
  return(list(C_opt = design.C, D_opt = design.D, EQD = design.EQD, planningValue = par.est))
}


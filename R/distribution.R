.fixInf <- function(data) {
  if (any(is.infinite(data))) {
    data <- apply(data, 2, function(x) {
      if (any(is.infinite(x))) {
        x[ind <- which(is.infinite(x))]=NA
        x[ind]=max(x, na.rm = TRUE) + 1
      }
      x
    })
  }
  data
}


# cor a symmetric correlation matrix
# sds standard deviations of the resulting covariance.
cor2cov <- function(cor, sds) {
    if (length(sds) != length(diag(cor))) stop("inputs are of mismatched dimension")
    cor * sds * rep(sds, each=nrow(cor))
}


#Sigma is correlation matrix
rmvnorm <- function (n, mu, Sigma, tol = 1e-06, empirical = TRUE){
  p=length(mu)
  if (!all(dim(Sigma) == c(p, p))) 
    stop("incompatible arguments")
  eS=eigen(Sigma, symmetric = TRUE)
  ev=eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) 
    stop("'Sigma' is not positive definite")
  X=matrix(rnorm(p * n), n)
  if (empirical) {
    X=scale(X, TRUE, FALSE)
    X=X %*% svd(X, nu = 0, nv = length(mu))$v
    X=scale(X, FALSE, TRUE)
  }
  X=drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  t(X)
}

# negative binomial (NB)
# Sigma is covariance matrix
rmvnbinom <- function(n, mu, Sigma, phi) {
  if(max(abs(Sigma))>1){
  	Cor=cov2cor(Sigma)
  }else{
 	Cor=Sigma
  }
  normd=rmvnorm(n, rep(0, nrow(Cor)), Sigma = Cor)
  unif=pnorm(normd)
  dat=t(qnbinom(t(unif), mu = mu, size = 1/phi))
  dat= .fixInf(dat)
  dat
}


# zero-inflated negative binomial (ZINB)
dzinbinom <- function(x, mu, theta, size, pi, log = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta=size
  rval=log(1 - pi) + dnbinom(x, mu = mu, size = theta, log = TRUE)
  if(any(x0 <- (x == 0L))) rval[x0]=log(exp(rval) + pi)[x0]
  if(log) rval else exp(rval)
}

pzinbinom <- function(q, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta=size
  rval=log(1 - pi) + pnbinom(q, mu = mu, size = theta, lower.tail = lower.tail, log.p = TRUE)
  if(any(q0 <- (is.finite(rval) & (lower.tail | q < 0)))) rval[q0]=log(exp(rval) + pi)[q0]
  if(log.p) rval else exp(rval)
}

qzinbinom <- function(p, mu, theta, size, pi, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta=size
  if(log.p) p=exp(p)
  if(!lower.tail) p=1 - p
  p=pmax(0, (p - pi)/(1 - pi))
  rval=qnbinom(p, mu = mu, size = theta, lower.tail = TRUE, log.p = FALSE)
  rval
}


rzinbinom <- function(n, mu, theta, size, pi) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(!missing(theta) & !missing(size)) stop("only 'theta' or 'size' may be specified")
  if(!missing(size)) theta=size
  rval=rnbinom(n, mu = mu, size = theta)
  rval[runif(n) < pi]=0
  rval
}

# Sigma is covariance matrix
rmvzinbinom = function(n, mu, Sigma, phi,p0) {
   if(max(abs(Sigma))>1){
  	Cor=cov2cor(Sigma)
  }else{
 	Cor=Sigma
  }
  normd=rmvnorm(n, rep(0, nrow(Cor)), Sigma = Cor)  #The normal-to-anything framework
  unif=pnorm(normd)
  dat=qzinbinom(p=t(unif), mu=mu, theta=1/phi,pi=p0)
  dat=t(matrix(dat,ncol=nrow(unif)))
  dat= .fixInf(dat)
  dat
}

# zero-inflated poisson (ZIP)
dzipois <- function(x, lambda, pi, log = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  rval=log(1 - pi) + dpois(x, lambda = lambda, log = TRUE)
  if(any(x0 <- (x == 0L))) rval[x0]=log(exp(rval) + pi)[x0]
  if(log) rval else exp(rval)
}

pzipois <- function(q, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  rval=log(1 - pi) + ppois(q, lambda = lambda, lower.tail = lower.tail, log.p = TRUE)
  if(any(q0 <- (is.finite(rval) & (lower.tail | q < 0)))) rval[q0]=log(exp(rval) + pi)[q0]
  if(log.p) rval else exp(rval)
}

qzipois <- function(p, lambda, pi, lower.tail = TRUE, log.p = FALSE) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  if(log.p) p=exp(p)
  if(!lower.tail) p=1 - p
  p=pmax(0, (p - pi)/(1 - pi))
  rval=qpois(p, lambda = lambda, lower.tail = TRUE, log.p = FALSE)
  rval
}

rzipois <- function(n, lambda, pi) {
  if(any(pi < 0) | any(pi > 1))  warning("'pi' must be in [0, 1]")
  rval=rpois(n, lambda = lambda)
  rval[runif(n) < pi]=0
  rval
}


rmvZIP = function(n, mu, Sigma,p0) {
  if(max(abs(Sigma))>1){
  	Cor=cov2cor(Sigma)
  }else{
 	Cor=Sigma
  }
  normd=rmvnorm(n, rep(0, nrow(Cor)), Sigma = Cor)  #The normal-to-anything framework
  unif=pnorm(normd)
  dat=qzipois(p=t(unif), lambda=mu,pi=p0)
  dat=t(matrix(dat,ncol=nrow(unif)))
  dat= .fixInf(dat)
  dat
}








MAGMA<-function(x,distrib=c('NB','ZINB','ZIP')){ #row is sample, col is taxon
  distrib=match.arg(distrib)
  networks=magma(data=x, distrib=distrib,criterion.select="stars")
  lambda=networks$opt.lambda
  network=magma(data=x, distrib=distrib, lambda = lambda) #with fixed lambda
  icov=network$icov
  icov=as.data.frame(icov)
  Sigma=prec2cov(icov)
  Sigma=cov2cor(Sigma)
  Sigma
}


###################################################################################
### the code is obtained  and adjusted from https://github.com/huayingfang/CCLasso

cclasso <- function(x, counts = FALSE, pseudo = 0.5, k_cv = 3, 
                    lam_int = c(1e-4, 1), k_max = 20, n_boot = 20) {
  n <- nrow(x);
  p <- ncol(x);

  if(counts) {
    x <- x + pseudo;
    x <- x / rowSums(x);
  }
  x <- log(x);
  vx2 <- stats::var(x);

  # Diagonal weight for loss function
  rmean_vx2 <- rowMeans(vx2);
  wd <- 1/diag(vx2 - rmean_vx2 - rep(rmean_vx2, each = p) + mean(rmean_vx2));
  wd2 <- sqrt(wd);
  
  # Some global parameters for optimization with single lambda
  rho <- 1;
  u_f <- eigen(diag(p) - 1/p)$vectors;
  wd_u <- (t(u_f) %*% (wd * u_f))[-p, -p];
  wd_u_eig <- eigen(wd_u);
  d0_wd <- 1 / ( (rep(wd_u_eig$values, each = p-1) + wd_u_eig$values) / 
    (2 * rho) + 1 );
  u0_wd <- wd_u_eig$vectors;

  # Golden section method for the selection of lambda (log10 scale)
  sigma <- vx2;
  lam_int2 <- log10(range(lam_int));
  a1 <- lam_int2[1]; 
  b1 <- lam_int2[2];
  # Store lambda and corresponding cross validation's loss
  lams <- NULL; 
  fvals <- NULL;
  # Two trial points in first 
  a2 <- a1 + 0.382 * (b1 - a1); 
  b2 <- a1 + 0.618 * (b1 - a1);
  fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
    sigma = sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
    wd2 = wd2);
  lams <- c(lams, b2); 
  fvals <- c(fvals, fb2$cv_loss);
  fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
    sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
    wd2 = wd2);
  lams <- c(lams, a2); 
  fvals <- c(fvals, fa2$cv_loss);
  # Error tolerance for convergence
  err_lam2 <- 1e-1 * max(1, lam_int2);
  err_fval <- 1e-4;
    
  err <- b1 - a1;
  k <- 0;
  while(err > err_lam2 && k < k_max) {
    fval_max <- max(fa2$cv_loss, fb2$cv_loss);

    if(fa2$cv_loss > fb2$cv_loss) {
      a1 <- a2;      
      a2 <- b2;
      fa2 <- fb2;
      b2 <- a1 + 0.618 * (b1 - a1);
      fb2 <- cv_loss_cclasso(lambda2 = 10^b2 / rho, x = x, k_cv = k_cv, 
        sigma = fa2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
        wd2 = wd2);

      lams <- c(lams, b2); 
      fvals <- c(fvals, fb2$cv_loss);
    } else {
      b1 <- b2;
      b2 <- a2;
      fb2 <- fa2;
      a2 <- a1 + 0.382 * (b1 - a1);
      fa2 <- cv_loss_cclasso(lambda2 = 10^a2 / rho, x = x, k_cv = k_cv, 
        sigma = fb2$sigma, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd, 
        wd2 = wd2);

      lams <- c(lams, a2); 
      fvals <- c(fvals, fa2$cv_loss);
    }
    fval_min <- min(fa2$cv_loss, fb2$cv_loss);      

    k <- k + 1;
    err <- b1 - a1;
    if(abs(fval_max - fval_min) / (1 + fval_min) <= err_fval) {
      break;
    }
  }
  info_cv <- list(lams = lams, fvals = fvals, k = k + 2, 
    lam_int = 10^c(a1, b1)); 
  if(a1 == lam_int2[1] || b1 == lam_int2[2]) {
    cat("WARNING:\n", "\tOptimal lambda is near boundary! ([", 10^a1, ",", 
      10^b1, "])\n", sep = "");
  }

  lambda <- 10^((a2 + b2)/2);
  # Bootstrap for cclasso
  lambda2 <- lambda / rho;
  info_boot <- boot_cclasso(x = x, sigma = fb2$sigma, lambda2 = lambda2, 
    n_boot = n_boot, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);

  return(list(var_w = info_boot$var_w, cor_w = info_boot$cor_w, 
    p_vals = info_boot$p_vals, lambda = lambda, info_cv = info_cv));
}
#-------------------------------------------------------------------------------
# Bootstrap for cclasso
boot_cclasso <- function(x, sigma, lambda2, n_boot = 20, 
                         wd, u_f, u0_wd, d0_wd) {
  n <- nrow(x);
  p <- ncol(x);

  # Store the result of bootstrap
  cors_boot <- matrix(0, nrow = p * (p - 1)/2, ncol = n_boot + 1);
  vars_boot <- matrix(0, nrow = p, ncol = n_boot + 1);
  cors_mat <- matrix(0, p, p);
  ind_low <- lower.tri(cors_mat);

  # Bootstrap procedure
  sam_boot <- matrix(sample(1:n, size = n * n_boot, replace = T), 
    ncol = n_boot);
  for(k in 1:n_boot) {
    ind_samp <- sam_boot[, k];
    sigma2 <- cclasso_sub(sigma = sigma, vx = var(x[ind_samp, ]), 
      lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
        
    vars_boot[, k] <- diag(sigma2);
    Is <- 1 / sqrt(vars_boot[, k]);
    cors_mat <- Is * sigma2 * rep(Is, each = p);
    cors_boot[, k] <- cors_mat[ind_low];
  }
  Is <- 1 / sqrt(diag(sigma));
  cors_mat <- sigma * Is * rep(Is, each = p);
  cors_boot[, n_boot + 1] <- cors_mat[ind_low];
  vars_boot[, n_boot + 1] <- diag(sigma);

  #----------------------------------------  
  # Variance estimation via bootstrap
  vars2 <- rowMeans(vars_boot);
  #----------------------------------------
  # Correlations' relationship for artificial null sample
  tol_cor <- 1e-3;
  sam_art0 <- matrix(rnorm(n * p), nrow = n) * rep(sqrt(vars2), each = n);
  cors_art0 <- cor(sam_art0)[ind_low];
  sam_art <- sam_art0 - log(rowSums(exp(sam_art0)));
  sigma_art <- cclasso_sub(sigma = sigma, vx = var(sam_art), 
    lambda2 = lambda2, wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
  Is <- 1 / sqrt(diag(sigma_art));
  cors_mat <- Is * sigma_art * rep(Is, each = p);
  cors_art2 <- cors_mat[ind_low];
  # Bias of estimation between absolute data and cclasso of compositional data
  cors0m2 <- log( ((1 + cors_art0) * (1 - cors_art2)) / ((1 + cors_art2) * 
    (1 - cors_art0)) );
  tmp <- abs(cors_art2) >= tol_cor;
  bias02 <- ifelse(sum(tmp), median(abs(cors0m2)[tmp]), 0);
  # Modification of estimation for cclasso
  cors2 <- log( (1 + cors_boot) / (1 - cors_boot) );    
  cors2mod <- (cors_boot >= tol_cor) * (cors2 + bias02) + 
    (cors_boot <= -tol_cor) * (cors2 - bias02);
  cors2mod <- 1 - rowMeans(2 / (exp(cors2mod) + 1));
  cors2_mat <- diag(p);
  cors2_mat[ind_low] <- cors2mod;
  cors2_mat <- t(cors2_mat);
  cors2_mat[ind_low] <- cors2mod;
  # P-values with null distribution of correlation estimations of absolute data
  p_vals <- pt(cors2mod * sqrt((n - 2) / (1 - cors2mod^2)), df = n - 2);
  p_vals <- ifelse(p_vals <= 0.5, p_vals, 1 - p_vals);
  pval_mat <- diag(p);
  pval_mat[ind_low] <- p_vals;
  pval_mat <- t(pval_mat);
  pval_mat[ind_low] <- p_vals;
  #----------------------------------------

  return(list(var_w = vars2, cor_w = cors2_mat, p_vals = pval_mat));    
}
#-------------------------------------------------------------------------------
# cross validation's loss of cclasso for single lambda
cv_loss_cclasso <- function(lambda2, x, k_cv, sigma, 
                            wd, u_f, u0_wd, d0_wd, wd2) {
  n <- nrow(x);
  p <- ncol(x);

  n_b <- floor(n / k_cv);
  cv_loss <- 0;
  for(k in 1:k_cv) {
    itest <- (n_b * (k - 1) + 1):(n_b * k);
    vxk <- stats::var(x[itest, ]);
    vx2k <- stats::var(x[-itest, ]);

    sigma <- cclasso_sub(sigma = sigma, vx = vx2k, lambda2 = lambda2, 
      wd = wd, u_f = u_f, u0_wd = u0_wd, d0_wd = d0_wd);
    
    dsig <- sigma - vxk;
    tmp <- rowMeans(dsig);
    dsig <- dsig - tmp - rep(tmp, each = p) + mean(tmp);
    cv_loss <- cv_loss + base::norm(wd2 * dsig, "F")^2;
  }

  return(list(cv_loss = cv_loss, sigma = sigma));
}
#-------------------------------------------------------------------------------
# cclasso for single lambda
cclasso_sub <- function(sigma, vx, lambda2, 
                        wd, u_f, u0_wd, d0_wd, 
                        k_max = 200, x_tol = 1e-4) {
  p <- ncol(sigma);
  sigma2 <- sigma;
  LAMBDA <- matrix(0, p, p);
  lambda2 <- matrix(lambda2, p, p);
  diag(lambda2) <- 0;

  k <- 0;
  err <- 1;
  while(err > x_tol && k < k_max) {
    # Update sigma
    x_sigma <- t(u_f) %*% ((sigma2 - vx) - LAMBDA) %*% u_f;
    x_sigma[-p,-p] <- u0_wd %*% ((t(u0_wd) %*% x_sigma[-p, -p] %*% u0_wd) * 
      d0_wd) %*% t(u0_wd);
    sigma_new <- vx + u_f %*% x_sigma %*% t(u_f);
    # Update sigma2
    A <- LAMBDA + sigma_new;
    sigma2_new <- (A > lambda2) * (A - lambda2) + (A < -lambda2) * 
      (A + lambda2);
    # Update Lambda
    LAMBDA <- LAMBDA + (sigma_new - sigma2_new);

    err <- max( abs(sigma_new - sigma)/(abs(sigma) + 1), 
      abs(sigma2_new - sigma2)/(abs(sigma2) + 1) ); 
    k <- k + 1;
    sigma <- sigma_new;
    sigma2 <- sigma2_new;
  }

  if(k >= k_max) {
    cat("WARNING of cclasso_sub:\n", "\tMaximum Iteration:", k_max, 
      "&& Relative error:", err, "!\n");
  }
  
  return(sigma);
}

CCLasso<-function(x){
  networks=cclasso(x,counts = T)
  Sigma=networks$cor_w
  Sigma
}


SparCC<-function(x){
  #library(SpiecEasi)
  networks=sparcc(x)
  Sigma=networks$Cov
  Sigma=cov2cor(Sigma)
  Sigma
}


SPIECEASI <-function(x){
  #library(SpiecEasi)
  networks=spiec.easi(x)
  Sigma=getOptCov(networks)
  Sigma=as.matrix(Sigma)
  Sigma=cov2cor(Sigma)
  Sigma
}




CCREPE<-function(x){
  x.rowsum=apply(x,1,sum)
  x.norm=x/x.rowsum
  networks=ccrepe(x=x.norm)
  Sigma=networks$sim.score
  diag(Sigma)=1
  Sigma
}



estSigma<-function(x,method=c('MAGMA','CCLasso','SparCC','SPIECEASI','CCREPE'),distrib=c('NB','ZINB','ZIP')){ #row is sample, col is OTU
  method=match.arg(method)
  if(method=='MAGMA') distrib=match.arg(distrib)
  est.method=switch(method,MAGMA=MAGMA, CCLasso=CCLasso, SparCC=SparCC,SPIECEASI=SPIECEASI,CCREPE=CCREPE)
  if(method=='MAGMA') Sigma=est.method(x,distrib=distrib)
  else Sigma=est.method(x)
  Sigma
}






<<<<<<< HEAD
=======


>>>>>>> Initial commit
fit.NB<-function(x){
  x0=as.data.frame(x)
  x0=as.matrix(x0)
  notu=ncol(x0)
  sfs=apply(x0,1,mean); sfs=sfs/median(sfs)
  params=matrix(NA,notu,ncol=2)
  colnames(params)=c('mu','phi')
<<<<<<< HEAD
  for (i in 1:notu) {
    if(sum(x0[,i])>2){
      fit=try(glm.nb(x0[,i]~1+offset(log(sfs))))
      if(methods::is(fit,"try-error")) {
        fit=try(glm(x0[,i]~1+offset(log(sfs)),family = poisson))
=======
  fits=list()
  for (i in 1:notu) {
    if(sum(x0[,i])>2){
      fit=try(glm.nb(x0[,i]~1+offset(log(sfs))),silent=T)
      fits[[i]]=fit
      if(methods::is(fit,"try-error")) {
        fit=try(glm(x0[,i]~1+offset(log(sfs)),family = poisson),silent=T)
>>>>>>> Initial commit
        mu=exp(coef(fit))
        params[i,]=c(mu,0)
      }else{
        phi=1/fit$theta
        mu=exp(coef(fit))
        params[i,]=c(mu,phi)
      }
    }else{
      fit=try(glm(x0[,i]~1+offset(log(sfs)),family = poisson))
      mu=exp(coef(fit))
      params[i,]=c(mu,0)
    }
  }
<<<<<<< HEAD
  list(otu.mat=x0,mu=params[,'mu'],phi=params[,'phi'])
}



=======
  list(mu=params[,'mu'],phi=params[,'phi'],fits=fits)
  
}


>>>>>>> Initial commit
fit.ZINB<-function(x){
  x0=as.data.frame(x)
  x0=as.matrix(x0)
  notu=ncol(x0)
  sfs=apply(x0,1,mean); sfs=sfs/median(sfs)
  params=matrix(NA,notu,ncol=3)
  colnames(params)=c('mu','phi','p0')
<<<<<<< HEAD
=======
  fits=list()
>>>>>>> Initial commit
  for (i in 1:notu) {
    if(sum(x0[,i])>5){
      ctrl=pscl::zeroinfl.control(method = "L-BFGS-B")
      ctrl$reltol=NULL
      ctrl$factr=1e-3/.Machine$double.eps
      fit=try(pscl::zeroinfl(x0[,i]~1+offset(log(sfs))|1, dist = "negbin", link = "logit" , control=ctrl), silent=T)
<<<<<<< HEAD
=======
      fits[[i]]=fit
>>>>>>> Initial commit
      if(methods::is(fit,"try-error")) {
        fit=try(glm.nb(x0[,i]~1+offset(log(sfs))))
        phi=1/fit$theta
        mu=exp(coef(fit))
        params[i,]=c(mu,phi,0)
      }else{
        mu=exp(fit$coefficients$count) 
        theta=fit$theta
        p0=plogis(fit$coefficients$zero)
        params[i,]=c(mu,1/theta,p0)
      }
    }else{
      fit=try(glm.nb(x0[,i]~1+offset(sfs)))
      phi=1/fit$theta
      mu=exp(coef(fit))
      params[i,]=c(mu,phi,0)
    }
  }
<<<<<<< HEAD
  list(otu.mat=x0,mu=params[,'mu'],phi=params[,'phi'],p0=params[,'p0'])
=======
  list(mu=params[,'mu'],phi=params[,'phi'],p0=params[,'p0'],fits=fits)
>>>>>>> Initial commit
}



fit.ZIP<-function(x){
  x0=as.data.frame(x)
  x0=as.matrix(x0)
  notu=ncol(x0)
  sfs=apply(x0,1,mean); sfs=sfs/median(sfs)
  params=matrix(NA,notu,ncol=2)
  colnames(params)=c('mu','p0')
<<<<<<< HEAD
=======
  logliks=aics=rep(NA,notu)
  fits=list()
  
>>>>>>> Initial commit
  for (i in 1:notu) {
    if(sum(x0[,i])>5){
      ctrl =pscl::zeroinfl.control(method = "L-BFGS-B")
      ctrl$reltol= NULL
      ctrl$factr=1e-3/.Machine$double.eps
      fit=try(pscl::zeroinfl(x0[,i]~1+offset(log(sfs))|1, dist = "poisson", link = "logit" , control=ctrl), silent=T)
<<<<<<< HEAD
=======
      fits[[i]]=fit
>>>>>>> Initial commit
      if(methods::is(fit,"try-error")) {
        fit=try(glm(x0[,i]~1+offset(log(sfs)),family = quasipoisson))
        mu=exp(coef(fit))
        params[i,]=c(mu,0)
      }else{
        mu=exp(fit$coefficients$count) 
        p0=plogis(fit$coefficients$zero)
        params[i,]=c(mu,p0)
      }
    }else{
      fit=try(glm(x0[,i]~1+offset(log(sfs)),family = quasipoisson))
      mu=exp(coef(fit))
      params[i,]=c(mu,0)
    }
  }
<<<<<<< HEAD
  list(otu.mat=x0,mu=params[,'mu'],p0=params[,'p0'])
}



=======
  list(mu=params[,'mu'],p0=params[,'p0'],fits=fits)
}


>>>>>>> Initial commit
fit.DM<-function(x){
  params=list()
  res=dirmult(x)
  params$mu0=res$pi
  params$phi0=res$theta
<<<<<<< HEAD
=======
  #params$otu.mat=x
>>>>>>> Initial commit
  params
}


<<<<<<< HEAD
=======

>>>>>>> Initial commit
estParams<-function(x,Sigma=NULL,method=c('MAGMA','CCLasso','SparCC','SPIECEASI','CCREPE'),distrib=c('NB','ZINB','ZIP','DM')){
  distrib=match.arg(distrib)
  method=match.arg(method)
  if(is.null(Sigma) & distrib!='DM') Sigma=estSigma(x,method,distrib)
  params=list()
  params$Sigma=Sigma
  fit.distrib=switch(distrib,ZIP=fit.ZIP, NB=fit.NB, ZINB=fit.ZINB,DM=fit.DM)
  out=fit.distrib(x)
  params$otu.mat=out$otu.mat
  params$mu=out$mu
  params$phi=out$phi
  params$p0=out$p0
  params
}

















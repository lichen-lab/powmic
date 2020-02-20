
setParams.NB <- function(nTaxa=1000,p.DA=0.05,Sigma,lmu0, lphi0,lfc.mu,depth=1,
                         sim.seed=123456) {
  
  if(missing(sim.seed))
    sim.seed = 1
  set.seed(sim.seed)
  
  if(!is.null(Sigma)){
    nTaxa=nrow(Sigma)
  }
  if(length(lmu0)!=length(lphi0)){
    stop('length of lmu0 and lphi0 should be the same')
  }
  if(length(lmu0)==1) { 
    lmu0=rep(lmu0, nTaxa)
    lphi0=rep(lphi0, nTaxa)
  }else if (length(lmu0)!=nTaxa) {
    id=sample(1:length(lmu0),nTaxa, replace=TRUE)
    lmu0=lmu0[id]
    lphi0=lphi0[id]
  } 
  
  list(nTaxa=nTaxa, p.DA=p.DA, lmu0=lmu0, lphi0=lphi0,lfc.mu=lfc.mu,Sigma=Sigma,
       depth=depth,
       sim.seed=sim.seed)
}


simCounts.NB<-function(sim.params, n1, n2){
  
  designs=c(rep(0, n1), rep(1, n2))
  lmu0=sim.params$lmu0
  lphi0=sim.params$lphi0
  lfc.mu=sim.params$lfc.mu
  nlfc.mu=length(unique(lfc.mu))
  Sigma=sim.params$Sigma
  depth=sim.params$depth
  nTaxa=sim.params$nTaxa
  p.DA=sim.params$p.DA
  nDA=round(nTaxa*p.DA)
  lfc.mu1 = lfc.mu2= numeric(nTaxa)
  
  if(length(lmu0)!=nTaxa){
    id=sample(1:length(lmu0),nTaxa,replace=T)
    lmu0=est$lmu[id]
    lphi0=est$lphi[id]
  }
  
  #if(comp=='no'){
    DAid = sample(nTaxa,nDA)
    DAid1=sample(DAid,length(DAid)/2)
    DAid2=setdiff(DAid,DAid1)
    if(nlfc.mu==1){
     lfc.mu1[DAid1] = rep(lfc.mu,length(DAid1))
     lfc.mu2[DAid2] = rep(lfc.mu,length(DAid2))
    }else{
     lfc.mu1[DAid1] = sample(lfc.mu,length(DAid1),replace=T)
      lfc.mu2[DAid2] =sample(lfc.mu,length(DAid2),replace=T)
    }
   
  # }else if(comp=='yes'){
  #   #DAid = sample(nTaxa,nDA)
  #   #lfc.mu2[DAid] = rep(lfc.mu,nDA)
  #   DAid=sort(lmu0,index.return=T,decreasing=T)$ix
  # 
  #   DAid=sample(DAid[1:(length(lmu0)/2)],nDA)
  #   
  #   #DAid=DAid[1:nDA]
  #   lfc.mu2[DAid] = rep(lfc.mu,nDA)
  # }
  
  if(is.null(Sigma)){
    mu1=matrix(exp(lmu0+rep(lfc.mu1*depth,n1)),ncol=n1)
    mu2=matrix(exp(lmu0+rep(lfc.mu2*depth,n2)),ncol=n2)
    mu=cbind(mu1,mu2)
    x=rnbinom(n=length(mu), mu=mu, size=1/exp(lphi0))
    x=matrix(x,ncol=n1+n2)
  }else{
    mu1=exp(lmu0+lfc.mu1)
    mu2=exp(lmu0+lfc.mu2)
    x1=rmvnbinom(n1, mu1*depth, Sigma, exp(lphi0))
    x2=rmvnbinom(n2, mu2*depth, Sigma, exp(lphi0))
    x=cbind(t(x1),t(x2))
  }
  
  list(counts=x, designs=designs, DAid=DAid, sim.params=sim.params)
  
}



setParams.ZIP <- function(nTaxa=1000,p.DA=0.05,Sigma,lmu0,lp0,lfc.mu,
                          depth=1,
                          sim.seed=123456) {
  
  if(missing(sim.seed))
    sim.seed = 1
  set.seed(sim.seed)
  
  if(!is.null(Sigma)){
    nTaxa=nrow(Sigma)
  }
  if(length(lmu0)!=length(lp0)){
    stop('length of lmu0 and lp0 should be the same')
  }
  if(length(lmu0)==1) { 
    lmu0=rep(lmu0, nTaxa)
    lp0=rep(lp0,nTaxa)
  }else if (length(lmu0)!=nTaxa) {
    id=sample(1:length(lmu0),nTaxa, replace=TRUE)
    lmu0=lmu0[id]
    lp0=lp0[id]
    
  } 
  
  list(nTaxa=nTaxa, p.DA=p.DA, lmu0=lmu0,lp0=lp0,lfc.mu=lfc.mu,Sigma=Sigma,
       depth=depth,
       sim.seed=sim.seed)
}


simCounts.ZIP<-function(sim.params, n1, n2){
  
  designs=c(rep(0, n1), rep(1, n2))
  lmu0=sim.params$lmu0
  lp0=sim.params$lp0
  lfc.mu=sim.params$lfc.mu
  nlfc.mu=length(unique(lfc.mu))
  Sigma=sim.params$Sigma
  depth=sim.params$depth
  nTaxa=sim.params$nTaxa
  p.DA=sim.params$p.DA
  nDA=round(nTaxa*p.DA)
  lfc.mu1 = lfc.mu2= numeric(nTaxa)
  
  if(length(lmu0)!=nTaxa){
    id=sample(1:length(lmu0),nTaxa,replace=T)
    lmu0=est$lmu[id]
    lp0=est$lp0[id]
  }
  
  
  #if(comp=='no'){
    DAid = sample(nTaxa,nDA)
    DAid1=sample(DAid,length(DAid)/2)
    DAid2=setdiff(DAid,DAid1)
    if(nlfc.mu==1){
     lfc.mu1[DAid1] = rep(lfc.mu,length(DAid1))
     lfc.mu2[DAid2] = rep(lfc.mu,length(DAid2))
    }else{
     lfc.mu1[DAid1] = sample(lfc.mu,length(DAid1),replace=T)
     lfc.mu2[DAid2] = sample(lfc.mu,length(DAid2),replace=T)
    }
  # }else if(comp=='yes'){
  #   #DAid = sample(nTaxa,nDA)
  #   #lfc.mu2[DAid] = rep(lfc.mu,nDA)
  #   DAid=sort(lmu0,index.return=T)$ix
  #   DAid=sample(DAid[1:(length(lmu0)/3)],nDA)
  #   lfc.mu2[DAid] = rep(lfc.mu,nDA)
  # }
  # 
  if(is.null(Sigma)){
    mu1=matrix(exp(lmu0+rep(lfc.mu1*depth,n1)),ncol=n1)
    mu2=matrix(exp(lmu0+rep(lfc.mu2*depth,n2)),ncol=n2)
    mu=cbind(mu1,mu2)
    x=rzipois(n=length(mu),lambda=mu,pi=exp(lp0))
    x=matrix(x,ncol=n1+n2)
  }else{
    mu1=exp(lmu0+lfc.mu1)
    mu2=exp(lmu0+lfc.mu2)
    x1=rmvZIP(n1, mu1*depth, Sigma, exp(lp0))
    x2=rmvZIP(n2, mu2*depth, Sigma, exp(lp0))
    x=cbind(t(x1),t(x2))
  }
  
  list(counts=x, designs=designs, DAid=DAid, sim.params=sim.params)
  
}




setParams.ZINB <- function(nTaxa=1000,p.DA=0.05,Sigma,lmu0, lphi0,lp0,lfc.mu,
                           depth=1,
                           sim.seed=123456) {
  
  if(missing(sim.seed))
    sim.seed = 1
  set.seed(sim.seed)
  
  if(!is.null(Sigma)){
    nTaxa=nrow(Sigma)
  }
  if(length(lmu0)!=length(lphi0) & length(lphi0)!=length(lp0) ){
    stop('length of lmu0, lphi0 and lp0 should be the same')
  }
  if(length(lmu0)==1) { 
    lmu0=rep(lmu0, nTaxa)
    lphi0=rep(lphi0, nTaxa)
    lp0=rep(lp0,nTaxa)
  }else if (length(lmu0)!=nTaxa) {
    id=sample(1:length(lmu0),nTaxa, replace=TRUE)
    lmu0=lmu0[id]
    lphi0=lphi0[id]
    lp0=lp0[id]
    
  } 
  
  list(nTaxa=nTaxa, p.DA=p.DA, lmu0=lmu0, lphi0=lphi0,lp0=lp0,lfc.mu=lfc.mu,Sigma=Sigma,
       depth=depth,
       sim.seed=sim.seed)
}


simCounts.ZINB<-function(sim.params, n1, n2){
  
  
  designs=c(rep(0, n1), rep(1, n2))
  lmu0=sim.params$lmu0
  lphi0=sim.params$lphi0
  lp0=sim.params$lp0
  lfc.mu=sim.params$lfc.mu
  nlfc.mu=length(lfc.mu)
  Sigma=sim.params$Sigma
  depth=sim.params$depth
  nTaxa=sim.params$nTaxa
  p.DA=sim.params$p.DA
  nDA=round(nTaxa*p.DA)
  lfc.mu1 = lfc.mu2= numeric(nTaxa)
  
  if(length(lmu0)!=nTaxa){
    id=sample(1:length(lmu0),nTaxa,replace=T)
    lmu0=est$lmu[id]
    lphi0=est$lphi[id]
    lp0=est$lp0[id]
  }
  
  
  #if(comp=='no'){
    DAid = sample(nTaxa,nDA)
    DAid1=sample(DAid,length(DAid)/2)
    DAid2=setdiff(DAid,DAid1)
    if(nlfc.mu==1){
     lfc.mu1[DAid1] = rep(lfc.mu,length(DAid1))
     lfc.mu2[DAid2] = rep(lfc.mu,length(DAid2))
    }else{
     lfc.mu1[DAid1] = sample(lfc.mu,length(DAid1),replace=T)
      lfc.mu2[DAid2] =sample(lfc.mu,length(DAid2),replace=T)
    }
  # }else if(comp=='yes'){
  #   #DAid = sample(nTaxa,nDA)
  #   #lfc.mu2[DAid] = rep(lfc.mu,nDA)
  #   DAid=sort(lmu0,index.return=T)$ix
  #   DAid=sample(DAid[1:(length(lmu0)/3)],nDA)
  #   lfc.mu2[DAid] = rep(lfc.mu,nDA)
  # }
  # 
  if(is.null(Sigma)){
    mu1=matrix(exp(lmu0+rep(lfc.mu1*depth,n1)),ncol=n1)
    mu2=matrix(exp(lmu0+rep(lfc.mu2*depth,n2)),ncol=n2)
    mu=cbind(mu1,mu2)
    x=rzinbinom(n=length(mu),mu=mu,theta=1/exp(lphi0), pi=exp(lp0))
    x=matrix(x,ncol=n1+n2)
  }else{
    mu1=exp(lmu0+lfc.mu1)
    mu2=exp(lmu0+lfc.mu2)
    x1=rmvzinbinom(n1, mu1*depth, Sigma, exp(lphi0),exp(lp0))
    x2=rmvzinbinom(n2, mu2*depth, Sigma, exp(lphi0),exp(lp0))
    x=cbind(t(x1),t(x2))
  }
  
  list(counts=x, designs=designs, DAid=DAid, sim.params=sim.params)
  
}




setParams.DM <- function(nTaxa=1000, p.DA=0.05,lmu0, lphi0, lfc.mu,
                         libsize=libsize,depth=depth,sim.seed=123456) {
  
  if(missing(sim.seed))
    sim.seed = 1
  set.seed(sim.seed)
  
  if(length(lmu0)==1) { 
    lmu0=rep(lmu0, nTaxa)
  }else if (length(lmu0)!=nTaxa) {
    id=sample(1:length(lmu0),nTaxa, replace=TRUE)
    lmu0=lmu0[id]
  } 
  
  list(nTaxa=nTaxa, p.DA=p.DA, lmu0=lmu0, lphi0=lphi0,lfc.mu=lfc.mu,
       libsize=libsize,depth=depth,
       sim.seed=sim.seed)
}



simCounts.DM<-function(sim.params, n1, n2){
  
  designs=c(rep(0, n1), rep(1, n2))
  nTaxa = sim.params$nTaxa
  p.DA = sim.params$p.DA
  nDA=round(nTaxa*p.DA)
  lfc.mu=sim.params$lfc.mu
  nlfc.mu=length(lfc.mu)
  lfc.mu1 = lfc.mu2= numeric(nTaxa)
  mu0=exp(sim.params$lmu0)
  mu0=mu0/sum(mu0)
  phi0=exp(sim.params$lphi0)
  libsize=sim.params$libsize
  depth=sim.params$depth
  
  #if(comp=='no'){
    DAid = sample(nTaxa,nDA)
    DAid1=sample(DAid,length(DAid)/2)
    DAid2=setdiff(DAid,DAid1)
    if(nlfc.mu==1){
     lfc.mu1[DAid1] = rep(lfc.mu,length(DAid1))
     lfc.mu2[DAid2] = rep(lfc.mu,length(DAid2))
    }else{
     lfc.mu1[DAid1] = sample(lfc.mu,length(DAid1),replace=T)
      lfc.mu2[DAid2] =sample(lfc.mu,length(DAid2),replace=T)
    }
  # }else if(comp=='yes'){
  #   #DAid = sample(nTaxa,nDA)
  #   #lfc.mu2[DAid] = rep(lfc.mu,nDA)
  #   DAid=sort(lmu0,index.return=T)$ix
  #   DAid=sample(DAid[1:(length(lmu0)/3)],nDA)
  #   lfc.mu2[DAid] = rep(lfc.mu,nDA)
  # }
  
  g0=(1 - phi0) / phi0
  x1 = matrix(NA,nTaxa, n1)
  mu1=mu0*exp(lfc.mu1)
  mu1 = mu0*g0
  rmu1 = rdirichlet(n1, mu1)
  for (i in 1:n1) {
    x1[,i] = rmultinom(1,libsize, prob = rmu1[i,])[, 1]
  }
  
  x2 = matrix(NA,nTaxa, n2)
  mu2=mu0*exp(lfc.mu2)
  mu2=mu2/sum(mu2)
  mu2=mu2*g0
  rmu2 = rdirichlet(n2, mu2)
  for (i in 1:n2) {
    x2[,i] = rmultinom(1,libsize*depth, prob = rmu2[i,])[, 1]
  }
  
  x=cbind(x1,x2)
  
  list(counts=x, designs=designs, DAid=DAid, sim.params=sim.params)
  
}



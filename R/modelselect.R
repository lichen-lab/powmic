
modelselect.AIC<-function(x){
  
  nOTU=ncol(x)
  fit.NB=fit.NB(x)
  fit.ZINB=fit.ZINB(x)
  fit.ZIP=fit.ZIP(x)

  AIC.NB=AIC.ZINB=AIC.ZIP=rep(NA,nOTU)
  
  for(i in 1:nOTU){
    if(class(fit.NB$fits[[i]])!="try-error" & class(fit.NB$fits[[i]])!="NULL"){
      AIC.NB[i]=stats::AIC(fit.NB$fits[[i]])
    }
    if(class(fit.ZINB$fits[[i]])!="try-error" & class(fit.ZINB$fits[[i]])!="NULL"){
      AIC.ZINB[i]=stats::AIC(fit.ZINB$fits[[i]])
    }
    if(class(fit.ZIP$fits[[i]])!="try-error" & class(fit.ZIP$fits[[i]])!="NULL"){
      AIC.ZIP[i]=stats::AIC(fit.ZIP$fits[[i]])
    }
  }
  
  AIC.NB[is.na(AIC.NB)]=10000
  AIC.ZINB[is.na(AIC.ZINB)]=10000
  AIC.ZIP[is.na(AIC.ZIP)]=10000
  
  
  AICS=cbind(AIC.NB,AIC.ZINB,AIC.ZIP); colnames(AICS)=c('NB','ZINB','ZIP')
  tops=apply(AICS,1,which.min)
  percentage=table(tops)/sum(table(tops))
  
  df=data.frame(method=c('NB','ZINB','ZIP'),percentage=percentage)
  p=ggplot(df, aes(x=method, y=percentage, fill=method)) +
    geom_bar(stat="identity")+theme_minimal()+
    geom_text(aes(label=round(percentage,3)), vjust=1.6, color="white",
              position = position_dodge(0.9), size=3.5)+
    labs(title="AIC for model selection",x="Method", y = "Percentage")
  p
  
  
}



modelselect.vuong<-function(x){
  
  nOTU=ncol(x)
  fit.NB=fit.NB(x)
  fit.ZINB=fit.ZINB(x)
  fit.ZIP=fit.ZIP(x)

  ps1=ps2=ps4=rep(1,nOTU)
  for(i in 1:nOTU){
    vuong1 = try(nonnest2::vuongtest(fit.NB$fits[[i]],fit.ZINB$fits[[i]],nested=T), silent=T)
    if(!inherits(vuong1, 'try-error')) {
      ps1[i]=vuong1$p_LRT$A
    }
    vuong2 = try(nonnest2::vuongtest(fit.NB$fits[[i]],fit.ZIP$fits[[i]],nested=F), silent=T)
    if(!inherits(vuong2, 'try-error')) {
      ps2[i]=vuong2$p_LRT$A
    }
    vuong4 = try(nonnest2::vuongtest(fit.ZIP$fits[[i]],fit.NB$fits[[i]],nested=F), silent=T)
    if(!inherits(vuong4, 'try-error')) {
      ps4[i]=vuong4$p_LRT$A
    }
  }
  
  percentage=c(mean(ps1<0.05),mean(ps2<0.05),mean(ps4<0.05))
  df=data.frame(method=c('NB>ZINB','NB>ZIP','ZIP>NB'),percentage=percentage)
  p=ggplot(df, aes(x=method, y=percentage, fill=method)) +
    geom_bar(stat="identity")+theme_minimal()+
    geom_text(aes(label=round(percentage,3)), vjust=1.6, color="white",
              position = position_dodge(0.9), size=3.5)+
    labs(title="vuongtest",x="Method", y = "Percentage")
  p
  
  
  
}









